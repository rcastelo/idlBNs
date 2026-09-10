#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <string.h>
#include <math.h>
#include "cliquepick.h"

/* R_alloc, R_unif_index and unif_rand come from the R headers above */

/*
 * CLIQUE-PICKING
 *
 * Counting and uniform sampling of the acyclic moral orientations of a chordal
 * graph, after Wienoebst, Bannach & Liskiewicz, "Polynomial-Time Algorithms for
 * Counting and Sampling Markov Equivalent DAGs with Applications", JMLR 24
 * (2023). Replaces the enumerate-and-deduplicate routine in imec.c, which is
 * Theta(m! * #AMO) and becomes unusable at about m = 8.
 *
 * The formula is Proposition 25: for a rooted clique tree (T, r, iota) of a
 * UCCG G,
 *
 *     #AMO(G) = sum over nodes v of T of
 *                   phi(iota(v), FP(v, T)) * product over H in C_G(iota(v))
 *                                                of #AMO(H),
 *
 * where C_G(K) are the undirected connected components left after fixing the
 * clique K as source (Corollary 18), and FP(v, T) (Definition 23) collects the
 * intersections iota(x_i) ^ iota(x_{i+1}) along the root-to-v path that are
 * contained in iota(v). Those form a chain X_1 ( X_2 ( ... ( X_l (Lemma 24),
 * and phi -- the number of permutations of the clique that do not begin with
 * any X_i -- follows Lemma 22:
 *
 *     phi_i = (|S| - |X_i|)! - sum over j > i of (|X_j| - |X_i|)! * phi_j
 *
 * with X_0 = {} prepended, the answer being phi_0.
 *
 * Sampling (Section 5) reuses the same decomposition: draw a clique node with
 * probability proportional to its term in the sum, draw a permutation of that
 * clique uniformly among those with no forbidden prefix, emit it, and recurse
 * into the subproblems. The result is a vertex ORDER; orienting every edge from
 * the earlier to the later endpoint gives the AMO.
 *
 * SCOPE. Everything here works on one chain component at a time and represents
 * vertex sets as 64-bit masks, so components of up to 64 vertices are handled;
 * the caller checks that bound. The measured components at p = 500 reach 12.
 * Counts are doubles: exact to 2^53, which covers any component whose class the
 * caller could enumerate anyway, and for sampling only the ratios matter.
 *
 * Validated in tests/test_cliquepick.R against the brute-force enumerator for
 * every chordal graph up to 7 vertices and random ones up to 9, on the count,
 * on phi (against counting allowed permutations directly), and on uniformity.
 */

#define CP_MAXV 64

static inline int cp_popcount(cp_vset x) {
    int n = 0;
    while (x) { x &= x - 1; n++; }
    return n;
}
static inline int cp_first(cp_vset x) {           /* index of the lowest bit */
    int i = 0;
    while (!(x & 1u)) { x >>= 1; i++; }
    return i;
}

/* ------------------------------------------------------- maximum cardinality
 * search over `uni`, visiting the vertices of `prio` first. Fills order[] and
 * returns its length. O(m^2), which is ample at these sizes.
 */
static int
cp_mcs(const cp_vset *adj, int m, cp_vset uni, cp_vset prio, int *order) {
    int wt[CP_MAXV];
    cp_vset left = uni;
    for (int v = 0; v < m; v++) wt[v] = 0;
    int n = 0;
    while (left) {
        cp_vset cand = (prio & left) ? (prio & left) : left;
        int best = -1;
        for (cp_vset c = cand; c; c &= c - 1) {
            int v = cp_first(c);
            if (best < 0 || wt[v] > wt[best]) best = v;
        }
        order[n++] = best;
        left &= ~((cp_vset) 1 << best);
        for (cp_vset c = adj[best] & left; c; c &= c - 1) wt[cp_first(c)]++;
    }
    return n;
}

/* ------------------------------------------------------------- clique tree
 * Maximal cliques and parent pointers, built from an MCS order exactly as in
 * Berry & Simonet's Algorithm 4: walking the order, the clique of x is its
 * already-visited neighbourhood plus x, and a new clique node opens whenever
 * that neighbourhood differs from the current one. Returns the node count.
 */
static int
cp_clique_tree(const cp_vset *adj, int m, cp_vset uni, const int *order, int n,
               cp_vset *K, int *parent) {
    int inv[CP_MAXV], cliq[CP_MAXV];
    for (int i = 0; i < n; i++) inv[order[i]] = i;
    int s = -1;
    cp_vset visited = 0;
    for (int i = 0; i < n; i++) {
        int x = order[i];
        cp_vset S = adj[x] & visited & uni;
        if (s < 0 || K[s] != S) {
            s++;
            K[s] = S;
            parent[s] = -1;
            if (S) {                       /* attach under the clique of the
                                              latest-visited vertex of S */
                int k = -1;
                for (cp_vset c = S; c; c &= c - 1) {
                    int w = cp_first(c);
                    if (k < 0 || inv[w] > inv[k]) k = w;
                }
                parent[s] = cliq[k];
            }
        }
        K[s] |= (cp_vset) 1 << x;
        cliq[x] = s;
        visited |= (cp_vset) 1 << x;
    }
    return s + 1;
}

/* ------------------------------------------------------------------ C_G(K)
 * The undirected connected components of G_K[V \ K]: orient every edge out of
 * the clique K, close under Meek's rules, and take what stays undirected
 * (Corollary 18). dir[u] holds the heads of the arcs out of u.
 */
static int
cp_subproblems(const cp_vset *adj, int m, cp_vset uni, cp_vset K, cp_vset *out) {
    cp_vset dir[CP_MAXV];
    for (int v = 0; v < m; v++) dir[v] = 0;
    for (cp_vset c = K; c; c &= c - 1) {
        int u = cp_first(c);
        dir[u] |= adj[u] & uni & ~K;                    /* K is the source */
    }
    /* Meek closure. R1 and R2 suffice for a chordal graph oriented from a
       source clique, but R3 and R4 are cheap here and cost nothing to keep. */
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int a = 0; a < m; a++) {
            if (!((uni >> a) & 1)) continue;
            for (int b = 0; b < m; b++) {
                if (a == b || !((uni >> b) & 1)) continue;
                int undir_ab = ((adj[a] >> b) & 1) &&
                               !((dir[a] >> b) & 1) && !((dir[b] >> a) & 1);
                if (!undir_ab) continue;
                int orient = 0;
                for (int c = 0; c < m && !orient; c++) {
                    if (c == a || c == b || !((uni >> c) & 1)) continue;
                    /* R1: c -> a, c not adjacent to b */
                    if (((dir[c] >> a) & 1) && !((adj[c] >> b) & 1)) orient = 1;
                    /* R2: a -> c -> b */
                    if (!orient && ((dir[a] >> c) & 1) && ((dir[c] >> b) & 1)) orient = 1;
                }
                if (!orient) {
                    /* R3: c -> b, d -> b, a - c, a - d, c and d non-adjacent */
                    for (int c = 0; c < m && !orient; c++) {
                        if (c == a || c == b || !((uni >> c) & 1)) continue;
                        if (!(((dir[c] >> b) & 1) && ((adj[a] >> c) & 1) &&
                              !((dir[a] >> c) & 1) && !((dir[c] >> a) & 1))) continue;
                        for (int d = c + 1; d < m && !orient; d++) {
                            if (d == a || d == b || !((uni >> d) & 1)) continue;
                            if (((dir[d] >> b) & 1) && ((adj[a] >> d) & 1) &&
                                !((dir[a] >> d) & 1) && !((dir[d] >> a) & 1) &&
                                !((adj[c] >> d) & 1)) orient = 1;
                        }
                    }
                }
                if (orient) { dir[a] |= (cp_vset) 1 << b; changed = 1; }
            }
        }
    }
    /* connected components of what is still undirected, outside K */
    cp_vset rest = uni & ~K, seen = 0;
    int nc = 0;
    while (rest & ~seen) {
        cp_vset start = (rest & ~seen);
        int v = cp_first(start);
        cp_vset comp = (cp_vset) 1 << v, frontier = comp;
        while (frontier) {
            int w = cp_first(frontier);
            frontier &= frontier - 1;
            for (cp_vset c = adj[w] & rest & ~comp; c; c &= c - 1) {
                int z = cp_first(c);
                if (!((dir[w] >> z) & 1) && !((dir[z] >> w) & 1)) {
                    comp |= (cp_vset) 1 << z;
                    frontier |= (cp_vset) 1 << z;
                }
            }
        }
        seen |= comp;
        if (cp_popcount(comp) >= 1) out[nc++] = comp;
    }

    /*
     * The components must come out in an order consistent with the arcs
     * BETWEEN them. By Corollary 18.2, two adjacent vertices lying in different
     * undirected components of G_K[V \ K] are already oriented, from the
     * earlier to the later in the order; so emitting the components in an
     * arbitrary order would contradict those arcs and produce an orientation
     * that is not an AMO at all. Counting never noticed -- it is a product over
     * components -- but sampling does. Topologically sort them.
     */
    int indeg[CP_MAXV];
    for (int a = 0; a < nc; a++) indeg[a] = 0;
    for (int a = 0; a < nc; a++)
        for (int b = 0; b < nc; b++) {
            if (a == b) continue;
            int arc = 0;
            for (cp_vset c = out[a]; c && !arc; c &= c - 1) {
                int x = cp_first(c);
                if (dir[x] & out[b]) arc = 1;
            }
            if (arc) indeg[b]++;
        }
    cp_vset tmp[CP_MAXV];
    char placed[CP_MAXV];
    for (int a = 0; a < nc; a++) placed[a] = 0;
    for (int k = 0; k < nc; k++) {
        int pick = -1;
        for (int a = 0; a < nc; a++)
            if (!placed[a] && indeg[a] == 0) { pick = a; break; }
        if (pick < 0) {                    /* cannot happen for a valid G_K */
            for (int a = 0; a < nc; a++) if (!placed[a]) { pick = a; break; }
        }
        placed[pick] = 1;
        tmp[k] = out[pick];
        for (int b = 0; b < nc; b++) {
            if (placed[b]) continue;
            int arc = 0;
            for (cp_vset c = out[pick]; c && !arc; c &= c - 1) {
                int x = cp_first(c);
                if (dir[x] & out[b]) arc = 1;
            }
            if (arc) indeg[b]--;
        }
    }
    for (int a = 0; a < nc; a++) out[a] = tmp[a];
    return nc;
}

/* ------------------------------------------------------------------- phi */
static double
cp_fac(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

/*
 * Permutations of a clique of size csize that begin with none of the nested
 * sets of sizes fp[1..nfp-1]; fp[0] = 0. Lemma 22, memoised over the chain.
 */
static double
cp_phi(int csize, const int *fp, int nfp) {
    double memo[CP_MAXV + 2];
    for (int i = nfp - 1; i >= 0; i--) {
        double s = cp_fac(csize - fp[i]);
        for (int j = i + 1; j < nfp; j++)
            s -= cp_fac(fp[j] - fp[i]) * memo[j];
        memo[i] = s;
    }
    return memo[0];
}

/* --------------------------------------------------------- the memo table */

static int
cp_lookup(const cp_ctx *ctx, cp_vset key) {
    for (int i = 0; i < ctx->nnode; i++)
        if (ctx->node[i].key == key) return i;
    return -1;
}

static int
cp_new_node(cp_ctx *ctx, cp_vset key) {
    if (ctx->nnode == ctx->cap) {
        int nc = ctx->cap ? ctx->cap * 2 : 32;
        cp_node *nn = (cp_node *) R_alloc((size_t) nc, sizeof(cp_node));
        if (ctx->node) memcpy(nn, ctx->node, (size_t) ctx->nnode * sizeof(cp_node));
        ctx->node = nn; ctx->cap = nc;
    }
    cp_node *nd = &ctx->node[ctx->nnode];
    memset(nd, 0, sizeof(*nd));
    nd->key = key;
    return ctx->nnode++;
}

/*
 * cp_build -- Proposition 25, memoised over vertex subsets. Returns the node
 * id of `uni`, or -1 if anything went out of range.
 */
int
cp_build(cp_ctx *ctx, const cp_vset *adj, int m, cp_vset uni) {
    if (!ctx->ok) return -1;
    int hit = cp_lookup(ctx, uni);
    if (hit >= 0) return hit;

    ctx->adj = adj; ctx->m = m;
    int id = cp_new_node(ctx, uni);

    int order[CP_MAXV], parent[CP_MAXV];
    cp_vset K[CP_MAXV];
    int n = cp_mcs(adj, m, uni, 0, order);
    int ncl = cp_clique_tree(adj, m, uni, order, n, K, parent);

    cp_node *nd = &ctx->node[id];
    nd->ncl = ncl;
    nd->cl  = (cp_vset *) R_alloc((size_t) ncl, sizeof(cp_vset));
    nd->w   = (double *)  R_alloc((size_t) ncl, sizeof(double));
    nd->nsp = (int *)     R_alloc((size_t) ncl, sizeof(int));
    nd->sp  = (int **)    R_alloc((size_t) ncl, sizeof(int *));
    nd->nfp = (int *)     R_alloc((size_t) ncl, sizeof(int));
    nd->fp  = (cp_vset **)R_alloc((size_t) ncl, sizeof(cp_vset *));

    double total = 0.0;
    for (int v = 0; v < ncl; v++) {
        /* the subproblems C_G(K[v]) and their counts */
        cp_vset comps[CP_MAXV];
        int nc = cp_subproblems(adj, m, uni, K[v], comps);
        int *ids = (int *) R_alloc((size_t) (nc > 0 ? nc : 1), sizeof(int));
        double prod = 1.0;
        for (int c = 0; c < nc; c++) {
            int cid = cp_build(ctx, adj, m, comps[c]);
            if (cid < 0) { ctx->ok = 0; return -1; }
            ids[c] = cid;
            prod *= ctx->node[cid].total;
            /* re-fetch: the node array can move when the recursion grows it */
            nd = &ctx->node[id];
        }
        nd->nsp[v] = nc; nd->sp[v] = ids;

        /* FP(v, T): the intersections along the root-to-v path that are
           contained in K[v] (Definition 23); nested, so sorting by size and
           dropping duplicates gives the chain X_1 ( ... ( X_l */
        cp_vset chain[CP_MAXV];
        int nch = 0;
        for (int cur = v; parent[cur] >= 0; cur = parent[cur]) {
            cp_vset S = K[cur] & K[parent[cur]];
            if (S && (S & ~K[v]) == 0) {
                int dup = 0;
                for (int q = 0; q < nch; q++) if (chain[q] == S) { dup = 1; break; }
                if (!dup) chain[nch++] = S;
            }
        }
        for (int a = 1; a < nch; a++) {          /* insertion sort by size */
            cp_vset x = chain[a];
            int b = a - 1;
            while (b >= 0 && cp_popcount(chain[b]) > cp_popcount(x)) {
                chain[b + 1] = chain[b]; b--;
            }
            chain[b + 1] = x;
        }
        nd->nfp[v] = nch;
        nd->fp[v] = (cp_vset *) R_alloc((size_t) (nch > 0 ? nch : 1), sizeof(cp_vset));
        int fps[CP_MAXV + 2];
        fps[0] = 0;
        for (int a = 0; a < nch; a++) {
            nd->fp[v][a] = chain[a];
            fps[a + 1] = cp_popcount(chain[a]);
        }
        nd->cl[v] = K[v];
        nd->w[v]  = prod * cp_phi(cp_popcount(K[v]), fps, nch + 1);
        total += nd->w[v];
    }
    nd = &ctx->node[id];
    nd->total = total;
    if (!R_FINITE(total) || total < 0.0) { ctx->ok = 0; return -1; }
    return id;
}

double
cp_count(cp_ctx *ctx, int id) {
    return (id < 0 || !ctx->ok) ? NA_REAL : ctx->node[id].total;
}

/* ---------------------------------------------------------------- sampling */

/*
 * o[s], the paper's first-occurrence table: for a vertex first appearing in the
 * prefix set X_k, o[s] = |X_k|; for a vertex in none of them, |K| + 1. A
 * permutation is forbidden exactly when, scanning it, the running maximum of o
 * equals the position -- that is when its first i elements are precisely X_k.
 */
static int
cp_allowed(const int *perm, int len, const int *o) {
    int mx = 0;
    for (int i = 0; i < len; i++) {
        int ov = o[perm[i]];
        if (ov > mx) mx = ov;
        if (mx == i + 1) return 0;
        if (mx > len) return 1;
    }
    return 1;
}

void
cp_sample(cp_ctx *ctx, int id, int *ord, int *tick) {
    const cp_node *nd = &ctx->node[id];

    /* a clique node, with probability proportional to its term in the sum. A
       one-clique subproblem needs no draw, which keeps singleton components
       free of RNG contact and the two engines' streams aligned. */
    int j = 0;
    if (nd->ncl > 1) {
        double u = unif_rand() * nd->total, acc = 0.0;
        j = nd->ncl - 1;
        for (int i = 0; i < nd->ncl; i++) {
            acc += nd->w[i];
            if (u <= acc) { j = i; break; }
        }
    }

    int K[CP_MAXV], len = 0;
    for (cp_vset c = nd->cl[j]; c; c &= c - 1) K[len++] = cp_first(c);

    int o[CP_MAXV];
    for (int i = 0; i < len; i++) o[K[i]] = len + 1;
    for (int a = nd->nfp[j] - 1; a >= 0; a--) {      /* smallest set wins */
        int sz = cp_popcount(nd->fp[j][a]);
        for (cp_vset c = nd->fp[j][a]; c; c &= c - 1) o[cp_first(c)] = sz;
    }

    /* rejection sampling; the paper bounds the expected number of trials by 2 */
    int perm[CP_MAXV];
    for (int trial = 0; ; trial++) {
        for (int i = 0; i < len; i++) perm[i] = K[i];
        for (int i = len - 1; i > 0; i--) {          /* Fisher-Yates */
            int k = (int) R_unif_index((double) (i + 1));
            int t = perm[i]; perm[i] = perm[k]; perm[k] = t;
        }
        if (cp_allowed(perm, len, o)) break;
        if (trial > 100000)
            error("cp_sample: no allowed permutation after 100000 trials");
    }

    for (int i = 0; i < len; i++) ord[(*tick)++] = perm[i];
    for (int c = 0; c < nd->nsp[j]; c++)
        cp_sample(ctx, nd->sp[j][c], ord, tick);
}

/* ------------------------------------------------------- R entry points */

/*
 * The R entry points work component by component with LOCAL vertex indices, so
 * only a chain component has to fit in 64 bits, never the whole graph. They
 * mirror what C_dag_imec_sample() does on a DAG state, and walk the components
 * in the same order, so the two engines consume the random stream identically.
 */
static int
cp_read_components(SEXP A_R, int *pp, int **comp_of, int ***vsets, int **msz,
                   cp_vset **adjbuf, int *pnc) {
    SEXP dim = getAttrib(A_R, R_DimSymbol);
    if (dim == R_NilValue || LENGTH(dim) != 2 || INTEGER(dim)[0] != INTEGER(dim)[1])
        error("cliquepick: expected a square matrix");
    int p = INTEGER(dim)[0];
    const int *a = LOGICAL(A_R);
    *pp = p;
    char *seen = (char *) R_alloc((size_t) p, sizeof(char));
    int *stack = (int *) R_alloc((size_t) p, sizeof(int));
    int *loc = (int *) R_alloc((size_t) p, sizeof(int));
    memset(seen, 0, (size_t) p);
    int **vs = (int **) R_alloc((size_t) (p > 0 ? p : 1), sizeof(int *));
    int *ms = (int *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(int));
    cp_vset *ab = (cp_vset *) R_alloc((size_t) (p > 0 ? p : 1) * 64, sizeof(cp_vset));
    int nc = 0;
    for (int v0 = 0; v0 < p; v0++) {
        if (seen[v0]) continue;
        int top = 0, n = 0;
        int *cur = (int *) R_alloc((size_t) p, sizeof(int));
        seen[v0] = 1; stack[top++] = v0; cur[n++] = v0;
        while (top > 0) {
            int w = stack[--top];
            for (int z = 0; z < p; z++)
                if (!seen[z] && (a[w + (size_t) p * z] || a[z + (size_t) p * w])) {
                    seen[z] = 1; stack[top++] = z; cur[n++] = z;
                }
        }
        if (n < 2) continue;
        if (n > CP_MAXV) return 0;
        for (int i = 1; i < n; i++) {           /* ascending */
            int x = cur[i], j = i - 1;
            while (j >= 0 && cur[j] > x) { cur[j + 1] = cur[j]; j--; }
            cur[j + 1] = x;
        }
        for (int i = 0; i < n; i++) loc[cur[i]] = i;
        cp_vset *adj = ab + (size_t) nc * 64;
        for (int i = 0; i < n; i++) adj[i] = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j && (a[cur[i] + (size_t) p * cur[j]] ||
                               a[cur[j] + (size_t) p * cur[i]]))
                    adj[i] |= (cp_vset) 1 << j;
        vs[nc] = cur; ms[nc] = n; nc++;
    }
    *vsets = vs; *msz = ms; *adjbuf = ab; *pnc = nc; (void) comp_of;
    return 1;
}

/* C_cp_amo_count -- #AMO of an undirected chordal graph, NA if out of range */
SEXP
C_cp_amo_count(SEXP A_R) {
    void *vmax = vmaxget();
    int p, nc, **vs, *ms; cp_vset *ab;
    if (!cp_read_components(A_R, &p, NULL, &vs, &ms, &ab, &nc)) {
        vmaxset(vmax); return ScalarReal(NA_REAL);
    }
    double tot = 1.0;
    for (int c = 0; c < nc; c++) {
        cp_ctx ctx; memset(&ctx, 0, sizeof(ctx)); ctx.ok = 1;
        cp_vset uni = (ms[c] == 64) ? ~(cp_vset) 0 : (((cp_vset) 1 << ms[c]) - 1);
        int id = cp_build(&ctx, ab + (size_t) c * 64, ms[c], uni);
        if (id < 0) { vmaxset(vmax); return ScalarReal(NA_REAL); }
        tot *= ctx.node[id].total;
    }
    vmaxset(vmax);
    return ScalarReal(tot);
}

/*
 * C_cp_amo_sample -- a uniform AMO as a 1-based order over all p vertices.
 * Orienting each edge from the earlier to the later vertex gives the
 * orientation; vertices in no chain component carry no undirected edge, so
 * where they land never matters.
 */
SEXP
C_cp_amo_sample(SEXP A_R) {
    void *vmax = vmaxget();
    int p, nc, **vs, *ms; cp_vset *ab;
    if (!cp_read_components(A_R, &p, NULL, &vs, &ms, &ab, &nc)) {
        vmaxset(vmax); return R_NilValue;
    }
    int *pos = (int *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(int));
    char *inc = (char *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(char));
    memset(inc, 0, (size_t) p);
    int ord[CP_MAXV];
    if (nc > 0) GetRNGstate();
    for (int c = 0; c < nc; c++) {
        cp_ctx ctx; memset(&ctx, 0, sizeof(ctx)); ctx.ok = 1;
        cp_vset uni = (ms[c] == 64) ? ~(cp_vset) 0 : (((cp_vset) 1 << ms[c]) - 1);
        int id = cp_build(&ctx, ab + (size_t) c * 64, ms[c], uni);
        if (id < 0) { if (nc > 0) PutRNGstate(); vmaxset(vmax); return R_NilValue; }
        int tick = 0;
        cp_sample(&ctx, id, ord, &tick);
        for (int i = 0; i < tick; i++) { pos[vs[c][ord[i]]] = i; inc[vs[c][ord[i]]] = 1; }
    }
    if (nc > 0) PutRNGstate();

    /* emit a full order: each component's vertices in their sampled order,
       then everything else */
    SEXP ans = PROTECT(allocVector(INTSXP, p));
    int *o = INTEGER(ans); int k = 0;
    for (int c = 0; c < nc; c++)
        for (int i = 0; i < ms[c]; i++)
            for (int j = 0; j < ms[c]; j++)
                if (pos[vs[c][j]] == i) { o[k++] = vs[c][j] + 1; break; }
    for (int v = 0; v < p; v++) if (!inc[v]) o[k++] = v + 1;
    UNPROTECT(1);
    vmaxset(vmax);
    return ans;
}
