#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <R_ext/Random.h>
#include <string.h>
#include "dag.h"
#include "dag_R.h"
#include "cliquepick.h"

/*
 * EXACT INTERVENTIONAL-MARKOV-EQUIVALENCE-CLASS OPERATIONS, IN C
 *
 * A port of R/imec.R. Where rcar() explores the I-equivalence class of the
 * current DAG by a random walk of covered-arc reversals -- whose equilibrium
 * is proportional to the number of I-covered arcs of each member, not uniform
 * -- this samples the class exactly, and can enumerate it.
 *
 * The class factorises. By Proposition 4 of Wienoebst, Bannach & Liskiewicz
 * (JMLR 2023), after Hauser & Buehlmann (2012), the undirected components of
 * an I-essential graph are chordal and a DAG lies in the class iff it is
 * obtained by acyclic moral orientations of those components INDEPENDENTLY of
 * each other. So |[D]_I| is a product of per-component AMO counts and one
 * independent uniform draw per component is uniform on the class.
 *
 * FOUR ORDERING RULES, so that this engine and the R one stay bit-identical
 * and tests/test_search_engines.R keeps passing:
 *
 *   1. Arrows are visited in (head, tail) ascending order, matching R's
 *      which(G & !t(G), arr.ind = TRUE), which is column-major. The
 *      protected-arrow labelling updates flags in place, so the order in
 *      which arrows are relabelled within a pass is observable.
 *
 *   2. Chain components come out ordered by their smallest vertex, and each
 *      component's vertex list ascending -- R scans v = 1..p and takes the
 *      component of the first unseen vertex.
 *
 *   3. A component's edges are listed by (higher endpoint, lower endpoint)
 *      ascending, matching which(sub & upper.tri(sub), arr.ind = TRUE), and
 *      its AMOs are enumerated over permutations of the component's vertices
 *      in LEXICOGRAPHIC order with first-appearance de-duplication. The draw
 *      is an index into that list.
 *
 *   4. One R_unif_index() per component, in component order, matching R's
 *      sample.int(length(amos), 1L) inside imec.sample()'s loop, and drawn
 *      only after every component has been enumerated.
 */

#define IDL_AMO_EDGEMAX 32          /* m <= 8 gives at most 28 edges */

/* ---------------------------------------------------------------- targets */

/*
 * tprot[u * p + v] = 1 when some intervention target tells u and v apart,
 * i.e. |I ^ {u, v}| = 1. Those arrows are I-essential and never become
 * undirected. Mirrors .protects() in R/imec.R.
 */
static char *
build_tprot(SEXP tgt, int p) {
    char *tprot = (char *) R_alloc((size_t) p * p, sizeof(char));
    memset(tprot, 0, (size_t) p * p);
    if (tgt == R_NilValue)
        return tprot;
    R_xlen_t K = XLENGTH(tgt);
    char *in = (char *) R_alloc((size_t) p, sizeof(char));
    for (R_xlen_t k = 0; k < K; k++) {
        SEXP Ik = VECTOR_ELT(tgt, k);
        memset(in, 0, (size_t) p);
        R_xlen_t nk = XLENGTH(Ik);
        const int *iv = INTEGER(Ik);
        for (R_xlen_t j = 0; j < nk; j++)
            if (iv[j] != NA_INTEGER && iv[j] >= 1 && iv[j] <= p)
                in[iv[j] - 1] = 1;
        for (int u = 0; u < p; u++)
            for (int v = 0; v < p; v++)
                if (in[u] != in[v])
                    tprot[(size_t) u * p + v] = 1;
    }
    return tprot;
}

/* ------------------------------------------------- the I-essential graph */

#define PROTECTED     1
#define UNDECIDABLE   2
#define NOT_PROTECTED 3

/*
 * E[u * p + v] = 1 for an arc u -> v; an undirected edge sets both. Chickering's
 * protected-arrow labelling with target protection, i.e. pcalg's
 * EssentialGraph::replaceUnprotected() (src/greedy.cpp:981) extended by rule 1
 * above. Validated against pcalg::dag2essgraph in tests/test_imec.R.
 */
static void
iess_graph(const idl_dag *d, const char *tprot, char *E) {
    int p = d->p;
    memset(E, 0, (size_t) p * p);
    for (int u = 0; u < p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++)
            E[(size_t) u * p + chu->v[j]] = 1;
    }

    int na = d->nedges;
    if (na == 0)
        return;
    int *af = (int *) R_alloc((size_t) na, sizeof(int));   /* tail */
    int *at = (int *) R_alloc((size_t) na, sizeof(int));   /* head */
    int *fl = (int *) R_alloc((size_t) na, sizeof(int));
    char *act = (char *) R_alloc((size_t) na, sizeof(char));

    /* rule 1: (head, tail) ascending */
    int n = 0;
    for (int v = 0; v < p; v++)
        for (int u = 0; u < p; u++)
            if (E[(size_t) u * p + v]) {
                af[n] = u; at[n] = v;
                fl[n] = tprot[(size_t) u * p + v] ? PROTECTED : UNDECIDABLE;
                act[n] = 1;
                n++;
            }

    /*
     * Arrows grouped by head. The configuration scans below each look only at
     * the arrows into one vertex, so without this they would be O(arrows) each
     * and the whole labelling O(arrows^2) -- 490 us per call at p = 200. Heads
     * never change, so the index is built once.
     */
    int *hstart = (int *) R_alloc((size_t) p + 1, sizeof(int));
    int *hlist  = (int *) R_alloc((size_t) na, sizeof(int));
    { int *cnt = (int *) R_alloc((size_t) p + 1, sizeof(int));
      memset(cnt, 0, ((size_t) p + 1) * sizeof(int));
      for (int e = 0; e < n; e++) cnt[at[e] + 1]++;
      hstart[0] = 0;
      for (int v = 0; v < p; v++) hstart[v + 1] = hstart[v] + cnt[v + 1];
      memcpy(cnt, hstart, ((size_t) p) * sizeof(int));
      for (int e = 0; e < n; e++) hlist[cnt[at[e]]++] = e; }
#define FOR_HEAD(vv, ee) \
    for (int _q = hstart[vv], ee = (_q < hstart[(vv) + 1] ? hlist[_q] : -1); \
         _q < hstart[(vv) + 1]; _q++, ee = hlist[_q])

    /* an arrow in a v-structure is protected */
    for (int v = 0; v < p; v++)
        FOR_HEAD(v, a) {
            for (int _r = hstart[v]; _r < hstart[v + 1]; _r++) {
                int b = hlist[_r];
                if (b >= a) continue;
                int s1 = af[a], s2 = af[b];
                if (!E[(size_t) s1 * p + s2] && !E[(size_t) s2 * p + s1]) {
                    fl[a] = PROTECTED; fl[b] = PROTECTED;
                }
            }
        }

#define ADJ(u, v)   (E[(size_t)(u) * p + (v)] || E[(size_t)(v) * p + (u)])
#define ISPAR(u, v) (E[(size_t)(u) * p + (v)] && !E[(size_t)(v) * p + (u)])
#define ISNBR(u, v) (E[(size_t)(u) * p + (v)] && E[(size_t)(v) * p + (u)])

    char *wasund = (char *) R_alloc((size_t) na, sizeof(char));
    for (;;) {
        int nund = 0;
        for (int e = 0; e < n; e++) {
            wasund[e] = (char) (act[e] && fl[e] == UNDECIDABLE);
            if (wasund[e]) nund++;
        }
        if (nund == 0)
            break;

        for (int e = 0; e < n; e++) {
            if (!act[e] || fl[e] != UNDECIDABLE) continue;
            int a = af[e], b = at[e], f = NOT_PROTECTED;

            /* (a)  c -> a  with c not adjacent to b */
            FOR_HEAD(a, g) {
                if (f == PROTECTED) break;
                if (!act[g]) continue;
                if (!ADJ(af[g], b))
                    f = (fl[g] == PROTECTED) ? PROTECTED : UNDECIDABLE;
            }
            /* (c)  a -> c -> b */
            FOR_HEAD(b, g) {
                if (f == PROTECTED) break;
                if (!act[g]) continue;
                int c = af[g];
                if (ISPAR(a, c)) {
                    int h = -1;
                    FOR_HEAD(c, q) if (act[q] && af[q] == a) { h = q; break; }
                    f = (fl[g] == PROTECTED && h >= 0 && fl[h] == PROTECTED)
                        ? PROTECTED : UNDECIDABLE;
                }
            }
            /* (d)  c1 -> b, c2 -> b, both neighbours of a, c1 and c2 non-adjacent */
            FOR_HEAD(b, g1) {
                if (f == PROTECTED) break;
                if (!act[g1]) continue;
                for (int _r2 = hstart[b]; _r2 < hstart[b + 1] && f != PROTECTED; _r2++) {
                    int g2 = hlist[_r2];
                    if (g2 >= g1 || !act[g2]) continue;
                    int c1 = af[g1], c2 = af[g2];
                    if (ISNBR(a, c1) && ISNBR(a, c2) && !ADJ(c1, c2))
                        f = (fl[g1] == PROTECTED && fl[g2] == PROTECTED)
                            ? PROTECTED : UNDECIDABLE;
                }
            }
            fl[e] = f;                 /* in place: later arrows see it */
        }

        /* count only the arrows that LEFT the undecidable set this pass; an
           arrow protected from the start must not keep the loop alive */
        int settled = 0;
        for (int e = 0; e < n; e++) {
            if (!act[e]) continue;
            if (fl[e] == NOT_PROTECTED) {
                E[(size_t) at[e] * p + af[e]] = 1;   /* becomes undirected */
                act[e] = 0;
            }
            if (wasund[e] && fl[e] != UNDECIDABLE)
                settled++;
        }
        if (settled == 0)
            error("iess_graph: invalid graph (labelling did not converge)");
    }
}

/* ------------------------------------------------------ chain components */

/*
 * Label the vertices of the undirected part of E. comp[v] is the component
 * index of v or -1 when v has no undirected incident edge. Rule 2: components
 * are numbered in order of their smallest vertex.
 */
static int
uccgs(const char *E, int p, int *comp) {
    for (int v = 0; v < p; v++) comp[v] = -1;
    int nc = 0;
    int *stack = (int *) R_alloc((size_t) p, sizeof(int));
    for (int v = 0; v < p; v++) {
        if (comp[v] >= 0) continue;
        int any = 0;
        for (int w = 0; w < p; w++)
            if (w != v && E[(size_t) v * p + w] && E[(size_t) w * p + v]) { any = 1; break; }
        if (!any) continue;
        int top = 0;
        comp[v] = nc; stack[top++] = v;
        while (top > 0) {
            int x = stack[--top];
            for (int w = 0; w < p; w++)
                if (comp[w] < 0 && E[(size_t) x * p + w] && E[(size_t) w * p + x]) {
                    comp[w] = nc; stack[top++] = w;
                }
        }
        nc++;
    }
    return nc;
}

/* ------------------------------------------------------------------ AMOs */

typedef struct {
    int      m;                       /* vertices in the component          */
    int     *vs;                      /* their indices, ascending           */
    int      ne;                      /* edges                              */
    int     *ea, *eb;                 /* endpoints, local indices, ea < eb   */
    unsigned n;                       /* number of AMOs                     */
    unsigned cap;
    unsigned *mask;                   /* one bitmask per AMO, first-seen order */
} idl_comp;

/* bit e of a mask is 1 when the edge is oriented ea[e] -> eb[e] */
static int
amo_is_moral(const idl_comp *C, unsigned mask) {
    int m = C->m;
    /* parents, per vertex, as bitmasks over local indices */
    unsigned pa[32];
    unsigned adj[32];
    for (int i = 0; i < m; i++) { pa[i] = 0u; adj[i] = 0u; }
    for (int e = 0; e < C->ne; e++) {
        int a = C->ea[e], b = C->eb[e];
        adj[a] |= 1u << b; adj[b] |= 1u << a;
        if (mask & (1u << e)) pa[b] |= 1u << a; else pa[a] |= 1u << b;
    }
    for (int v = 0; v < m; v++)
        for (int i = 0; i < m; i++) {
            if (!(pa[v] & (1u << i))) continue;
            for (int j = 0; j < i; j++) {
                if (!(pa[v] & (1u << j))) continue;
                if (!(adj[i] & (1u << j))) return 0;      /* immorality */
            }
        }
    return 1;
}

/*
 * Storage is R_alloc'd up front at the exact upper bound m!, one AMO per
 * permutation before de-duplication, so nothing here has to be freed on an
 * error unwind.
 */
static unsigned
factorial(int m) {
    unsigned f = 1u;
    for (int i = 2; i <= m; i++) f *= (unsigned) i;
    return f;
}

/*
 * Enumerate the acyclic moral orientations of one component, in the order R
 * produces them (rule 3): over permutations of the component's vertices in
 * lexicographic order, orienting each edge from the earlier to the later
 * vertex, keeping first appearances. Orienting along a vertex order can never
 * produce a cycle, so only morality is tested. Returns 0 when the component is
 * too large to enumerate.
 */
static int
amos_of(idl_comp *C, int max_size) {
    int m = C->m;
    if (m > max_size || C->ne > IDL_AMO_EDGEMAX)
        return 0;
    int *perm = (int *) R_alloc((size_t) m, sizeof(int));
    int *pos  = (int *) R_alloc((size_t) m, sizeof(int));
    for (int i = 0; i < m; i++) perm[i] = i;
    C->cap = factorial(m);
    C->mask = (unsigned *) R_alloc((size_t) C->cap, sizeof(unsigned));
    C->n = 0u;

    /* first-appearance de-duplication; the class is small, so a linear scan
       over the masks kept so far is cheaper than a hash table */
    for (;;) {
        for (int i = 0; i < m; i++) pos[perm[i]] = i;
        unsigned mask = 0u;
        for (int e = 0; e < C->ne; e++)
            if (pos[C->ea[e]] < pos[C->eb[e]]) mask |= 1u << e;
        if (amo_is_moral(C, mask)) {
            int seen = 0;
            for (unsigned k = 0; k < C->n; k++)
                if (C->mask[k] == mask) { seen = 1; break; }
            if (!seen) C->mask[C->n++] = mask;
        }
        /* next permutation, lexicographic */
        int i = m - 2;
        while (i >= 0 && perm[i] >= perm[i + 1]) i--;
        if (i < 0) break;
        int j = m - 1;
        while (perm[j] <= perm[i]) j--;
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
        for (int a = i + 1, b = m - 1; a < b; a++, b--) {
            t = perm[a]; perm[a] = perm[b]; perm[b] = t;
        }
    }
    return 1;
}

/* ------------------------------------------------------- the whole class */

/*
 * Decompose the I-equivalence class of the current DAG. Returns the number of
 * components; *enumerable is 0 when some component is too large, in which case
 * the caller must fall back to the rcar() walk rather than sample a wrong
 * distribution.
 */
static int
decompose(const idl_dag *d, const char *tprot, int max_size,
          char *E, idl_comp **out, int *enumerable) {
    int p = d->p;
    iess_graph(d, tprot, E);
    int *comp = (int *) R_alloc((size_t) p, sizeof(int));
    int nc = uccgs(E, p, comp);
    idl_comp *cs = (idl_comp *) R_alloc((size_t) (nc > 0 ? nc : 1), sizeof(idl_comp));
    *enumerable = 1;
    for (int c = 0; c < nc; c++) {
        idl_comp *C = &cs[c];
        memset(C, 0, sizeof(*C));
        int m = 0;
        for (int v = 0; v < p; v++) if (comp[v] == c) m++;
        C->m = m;
        C->vs = (int *) R_alloc((size_t) m, sizeof(int));
        int *loc = (int *) R_alloc((size_t) p, sizeof(int));
        int k = 0;
        for (int v = 0; v < p; v++) if (comp[v] == c) { C->vs[k] = v; loc[v] = k; k++; }
        /* rule 3: edges by (higher endpoint, lower endpoint) ascending */
        int cap = m * (m - 1) / 2;
        C->ea = (int *) R_alloc((size_t) (cap > 0 ? cap : 1), sizeof(int));
        C->eb = (int *) R_alloc((size_t) (cap > 0 ? cap : 1), sizeof(int));
        C->ne = 0;
        for (int j = 0; j < m; j++)
            for (int i = 0; i < j; i++) {
                int u = C->vs[i], w = C->vs[j];
                if (E[(size_t) u * p + w] && E[(size_t) w * p + u]) {
                    C->ea[C->ne] = loc[u]; C->eb[C->ne] = loc[w]; C->ne++;
                }
            }
        if (!amos_of(C, max_size)) *enumerable = 0;
    }
    *out = cs;
    return nc;
}

/* rebuild the DAG from an edge list: remove everything, then add the new arcs
   in a topological order of the target graph, so no intermediate state can be
   cyclic and idl_dag_add_edge's precondition always holds */
static void
set_edges(idl_dag *d, const int *from, const int *to, int n) {
    int p = d->p;
    int nold = d->nedges;
    int *of = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    int *ot = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    int k = 0;
    for (int u = 0; u < p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++) { of[k] = u; ot[k] = chu->v[j]; k++; }
    }
    for (int e = 0; e < k; e++) idl_dag_remove_edge(d, of[e], ot[e]);

    int *indeg = (int *) R_alloc((size_t) p, sizeof(int));
    memset(indeg, 0, (size_t) p * sizeof(int));
    for (int e = 0; e < n; e++) indeg[to[e]]++;
    int *queue = (int *) R_alloc((size_t) p, sizeof(int));
    int qh = 0, qt = 0;
    for (int v = 0; v < p; v++) if (indeg[v] == 0) queue[qt++] = v;
    char *done = (char *) R_alloc((size_t) p, sizeof(char));
    memset(done, 0, (size_t) p);
    while (qh < qt) {
        int u = queue[qh++];
        done[u] = 1;
        for (int e = 0; e < n; e++) {
            if (from[e] != u) continue;
            idl_dag_add_edge(d, from[e], to[e]);
            if (--indeg[to[e]] == 0) queue[qt++] = to[e];
        }
    }
    for (int v = 0; v < p; v++)
        if (!done[v]) error("set_edges: target edge list is cyclic");
}

/* the member of the class selected by pick[c] in each component */
static void
apply_pick(idl_dag *d, const char *E, const idl_comp *cs, int nc,
           const unsigned *pick) {
    int p = d->p;
    int cap = d->nedges > 0 ? d->nedges : 1;
    int *from = (int *) R_alloc((size_t) cap, sizeof(int));
    int *to   = (int *) R_alloc((size_t) cap, sizeof(int));
    int n = 0;
    for (int u = 0; u < p; u++)
        for (int v = 0; v < p; v++)
            if (E[(size_t) u * p + v] && !E[(size_t) v * p + u]) {
                from[n] = u; to[n] = v; n++;
            }
    for (int c = 0; c < nc; c++) {
        const idl_comp *C = &cs[c];
        unsigned mask = C->mask[pick[c]];
        for (int e = 0; e < C->ne; e++) {
            int a = C->vs[C->ea[e]], b = C->vs[C->eb[e]];
            if (mask & (1u << e)) { from[n] = a; to[n] = b; }
            else                  { from[n] = b; to[n] = a; }
            n++;
        }
    }
    set_edges(d, from, to, n);
}

/* ------------------------------------------------------ R entry points */

/* ===================================================================== *
 * THE SPARSE PATH
 *
 * Sampling runs on every iteration of the search, so it is the one that has to
 * scale. The matrix form below -- a p x p char array, enumerated over all
 * ordered pairs -- costs Theta(p^2) per call whatever the graph's density, and
 * at p = 200 that scaffolding, not Clique-Picking, was most of the time.
 *
 * This representation carries the arcs themselves plus two CSR indexes: the
 * arcs grouped by head, for the configuration scans of the labelling, and the
 * skeleton with each vertex's neighbours ascending, for adjacency queries by
 * binary search. Everything below is O(V + E log deg); nothing iterates pairs.
 *
 * The enumeration used by the exhaustive escape still builds the matrix: it
 * runs only at a local optimum, is capped by the class size, and its cost is
 * dominated by scoring the members anyway.
 * ===================================================================== */

typedef struct {
    int   p, na;
    int  *af, *at;        /* arc tail and head                              */
    char *undir;          /* arc has become undirected                      */
    int  *hstart, *hlist; /* arcs by head                                   */
    int  *sstart, *slist; /* skeleton: neighbours of v, ascending           */
    int  *sarc;           /* ... and the arc joining them                   */
} iess;

/* index of v in u's neighbour list, or -1 */
static inline int
ie_find(const iess *G, int u, int v) {
    int lo = G->sstart[u], hi = G->sstart[u + 1] - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (G->slist[mid] == v) return mid;
        if (G->slist[mid] < v) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}
static inline int ie_adj(const iess *G, int u, int v) { return ie_find(G, u, v) >= 0; }
static inline int ie_ispar(const iess *G, int u, int v) {      /* u -> v */
    int k = ie_find(G, u, v);
    if (k < 0) return 0;
    int e = G->sarc[k];
    return !G->undir[e] && G->af[e] == u;
}
static inline int ie_isnbr(const iess *G, int u, int v) {
    int k = ie_find(G, u, v);
    return k >= 0 && G->undir[G->sarc[k]];
}

/*
 * Per-vertex membership in the target family, as a word array: an arc u -> v is
 * protected by the targets exactly when the two masks differ, which is the
 * |I ^ {u, v}| = 1 condition. Costs O(p * K / 64) instead of the p^2 matrix.
 */
static int
build_tmask(SEXP tgt, int p, uint64_t **maskp) {
    R_xlen_t K = (tgt == R_NilValue) ? 0 : XLENGTH(tgt);
    int nw = (int) ((K + 63) / 64); if (nw < 1) nw = 1;
    uint64_t *mask = (uint64_t *) R_alloc((size_t) p * nw, sizeof(uint64_t));
    memset(mask, 0, (size_t) p * nw * sizeof(uint64_t));
    for (R_xlen_t k = 0; k < K; k++) {
        SEXP Ik = VECTOR_ELT(tgt, k);
        const int *iv = INTEGER(Ik);
        R_xlen_t nk = XLENGTH(Ik);
        for (R_xlen_t j = 0; j < nk; j++)
            if (iv[j] != NA_INTEGER && iv[j] >= 1 && iv[j] <= p)
                mask[(size_t)(iv[j] - 1) * nw + (k >> 6)] |= (uint64_t) 1 << (k & 63);
    }
    *maskp = mask;
    return nw;
}
static inline int
tm_protects(const uint64_t *mask, int nw, int u, int v) {
    for (int w = 0; w < nw; w++)
        if (mask[(size_t) u * nw + w] != mask[(size_t) v * nw + w]) return 1;
    return 0;
}

/*
 * The I-essential graph, on arc lists. Same protected-arrow labelling as the
 * matrix version -- and the same arrow order, (head, tail) ascending -- so the
 * two agree arc for arc; tests/test_imec.R pins the matrix one against
 * pcalg::dag2essgraph and tests/test_c_imec.R pins this against that.
 */
static void
ie_build(const idl_dag *d, SEXP tgt, iess *G) {
    int p = d->p;
    G->p = p;
    int na = d->nedges;
    G->na = na;
    G->af = (int *) R_alloc((size_t) (na > 0 ? na : 1), sizeof(int));
    G->at = (int *) R_alloc((size_t) (na > 0 ? na : 1), sizeof(int));
    G->undir = (char *) R_alloc((size_t) (na > 0 ? na : 1), sizeof(char));
    G->hstart = (int *) R_alloc((size_t) p + 1, sizeof(int));
    G->hlist = (int *) R_alloc((size_t) (na > 0 ? na : 1), sizeof(int));
    G->sstart = (int *) R_alloc((size_t) p + 1, sizeof(int));
    G->slist = (int *) R_alloc((size_t) (2 * na > 0 ? 2 * na : 1), sizeof(int));
    G->sarc  = (int *) R_alloc((size_t) (2 * na > 0 ? 2 * na : 1), sizeof(int));

    /* arcs in (head, tail) ascending order */
    int *cnt = (int *) R_alloc((size_t) p + 1, sizeof(int));
    memset(cnt, 0, ((size_t) p + 1) * sizeof(int));
    for (int u = 0; u < p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++) cnt[chu->v[j] + 1]++;
    }
    G->hstart[0] = 0;
    for (int v = 0; v < p; v++) G->hstart[v + 1] = G->hstart[v] + cnt[v + 1];
    int *fill = (int *) R_alloc((size_t) p, sizeof(int));
    memcpy(fill, G->hstart, (size_t) p * sizeof(int));
    for (int u = 0; u < p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++) {
            int v = chu->v[j], e = fill[v]++;
            G->af[e] = u; G->at[e] = v;
        }
    }
    /* the tails within a head group must be ascending too */
    for (int v = 0; v < p; v++)
        for (int a = G->hstart[v] + 1; a < G->hstart[v + 1]; a++) {
            int x = G->af[a], b = a - 1;
            while (b >= G->hstart[v] && G->af[b] > x) { G->af[b + 1] = G->af[b]; b--; }
            G->af[b + 1] = x;
        }
    for (int e = 0; e < na; e++) { G->undir[e] = 0; G->hlist[e] = e; }

    /* skeleton CSR, neighbours ascending */
    memset(cnt, 0, ((size_t) p + 1) * sizeof(int));
    for (int e = 0; e < na; e++) { cnt[G->af[e] + 1]++; cnt[G->at[e] + 1]++; }
    G->sstart[0] = 0;
    for (int v = 0; v < p; v++) G->sstart[v + 1] = G->sstart[v] + cnt[v + 1];
    memcpy(fill, G->sstart, (size_t) p * sizeof(int));
    for (int e = 0; e < na; e++) {
        int u = G->af[e], v = G->at[e];
        G->slist[fill[u]] = v; G->sarc[fill[u]++] = e;
        G->slist[fill[v]] = u; G->sarc[fill[v]++] = e;
    }
    for (int v = 0; v < p; v++)
        for (int a = G->sstart[v] + 1; a < G->sstart[v + 1]; a++) {
            int xv = G->slist[a], xe = G->sarc[a], b = a - 1;
            while (b >= G->sstart[v] && G->slist[b] > xv) {
                G->slist[b + 1] = G->slist[b]; G->sarc[b + 1] = G->sarc[b]; b--;
            }
            G->slist[b + 1] = xv; G->sarc[b + 1] = xe;
        }

    if (na == 0) return;

    uint64_t *tmask; int nw = build_tmask(tgt, p, &tmask);
    int *fl = (int *) R_alloc((size_t) na, sizeof(int));
    char *wasund = (char *) R_alloc((size_t) na, sizeof(char));
    for (int e = 0; e < na; e++)
        fl[e] = tm_protects(tmask, nw, G->af[e], G->at[e]) ? PROTECTED : UNDECIDABLE;

    /* v-structures */
    for (int v = 0; v < p; v++)
        for (int a = G->hstart[v]; a < G->hstart[v + 1]; a++)
            for (int b = G->hstart[v]; b < a; b++)
                if (!ie_adj(G, G->af[a], G->af[b])) { fl[a] = PROTECTED; fl[b] = PROTECTED; }

    for (;;) {
        int nund = 0;
        for (int e = 0; e < na; e++) {
            wasund[e] = (char) (!G->undir[e] && fl[e] == UNDECIDABLE);
            if (wasund[e]) nund++;
        }
        if (!nund) break;
        for (int e = 0; e < na; e++) {
            if (G->undir[e] || fl[e] != UNDECIDABLE) continue;
            int a = G->af[e], b = G->at[e], f = NOT_PROTECTED;
            for (int g = G->hstart[a]; g < G->hstart[a + 1] && f != PROTECTED; g++) {
                if (G->undir[g]) continue;
                if (!ie_adj(G, G->af[g], b))
                    f = (fl[g] == PROTECTED) ? PROTECTED : UNDECIDABLE;
            }
            for (int g = G->hstart[b]; g < G->hstart[b + 1] && f != PROTECTED; g++) {
                if (G->undir[g]) continue;
                int c = G->af[g];
                if (ie_ispar(G, a, c)) {
                    int k = ie_find(G, a, c), h = (k >= 0) ? G->sarc[k] : -1;
                    f = (fl[g] == PROTECTED && h >= 0 && !G->undir[h] &&
                         fl[h] == PROTECTED) ? PROTECTED : UNDECIDABLE;
                }
            }
            for (int g1 = G->hstart[b]; g1 < G->hstart[b + 1] && f != PROTECTED; g1++) {
                if (G->undir[g1]) continue;
                for (int g2 = G->hstart[b]; g2 < g1 && f != PROTECTED; g2++) {
                    if (G->undir[g2]) continue;
                    int c1 = G->af[g1], c2 = G->af[g2];
                    if (ie_isnbr(G, a, c1) && ie_isnbr(G, a, c2) && !ie_adj(G, c1, c2))
                        f = (fl[g1] == PROTECTED && fl[g2] == PROTECTED)
                            ? PROTECTED : UNDECIDABLE;
                }
            }
            fl[e] = f;
        }
        int settled = 0;
        for (int e = 0; e < na; e++) {
            if (G->undir[e]) continue;
            if (fl[e] == NOT_PROTECTED) G->undir[e] = 1;
            if (wasund[e] && fl[e] != UNDECIDABLE) settled++;
        }
        if (!settled) error("ie_build: invalid graph (labelling did not converge)");
    }
}

/*
 * Chain components of the I-essential graph, over the undirected arcs only.
 * O(V + E), and it visits only the vertices that carry one.
 */
typedef struct { int m; int *v; cp_vset adj[64]; } chaincomp;

static int
ie_components(const iess *G, chaincomp *cc, int *ncomp) {
    int p = G->p;
    char *seen = (char *) R_alloc((size_t) p, sizeof(char));
    int *stack = (int *) R_alloc((size_t) p, sizeof(int));
    int *loc = (int *) R_alloc((size_t) p, sizeof(int));
    memset(seen, 0, (size_t) p);
    int nc = 0;
    for (int v0 = 0; v0 < p; v0++) {
        if (seen[v0]) continue;
        int has = 0;
        for (int k = G->sstart[v0]; k < G->sstart[v0 + 1]; k++)
            if (G->undir[G->sarc[k]]) { has = 1; break; }
        if (!has) { seen[v0] = 1; continue; }
        int top = 0, n = 0;
        int *vs = (int *) R_alloc((size_t) p, sizeof(int));
        seen[v0] = 1; stack[top++] = v0; vs[n++] = v0;
        while (top > 0) {
            int w = stack[--top];
            for (int k = G->sstart[w]; k < G->sstart[w + 1]; k++) {
                if (!G->undir[G->sarc[k]]) continue;
                int z = G->slist[k];
                if (!seen[z]) { seen[z] = 1; stack[top++] = z; vs[n++] = z; }
            }
        }
        if (n < 2) continue;
        if (n > 64) return 0;
        for (int a = 1; a < n; a++) {
            int x = vs[a], b = a - 1;
            while (b >= 0 && vs[b] > x) { vs[b + 1] = vs[b]; b--; }
            vs[b + 1] = x;
        }
        chaincomp *C = &cc[nc];
        C->m = n; C->v = vs;
        for (int a = 0; a < n; a++) loc[vs[a]] = a;
        for (int a = 0; a < n; a++) C->adj[a] = 0;
        for (int a = 0; a < n; a++)
            for (int k = G->sstart[vs[a]]; k < G->sstart[vs[a] + 1]; k++)
                if (G->undir[G->sarc[k]])
                    C->adj[a] |= (cp_vset) 1 << loc[G->slist[k]];
        nc++;
    }
    *ncomp = nc;
    return 1;
}

/*
 * Draw a uniform member of the class and re-orient the DAG into it.
 *
 * pos[] holds a position for every vertex; only the relative order inside a
 * chain component ever decides an orientation, so vertices outside them keep
 * position 0. Only the arcs whose direction actually changes are touched.
 */
static int
sample_and_apply(idl_dag *d, const iess *G, int rng) {
    int p = d->p;
    chaincomp *cc = (chaincomp *) R_alloc((size_t) (p > 0 ? p : 1), sizeof(chaincomp));
    int nc = 0;
    if (!ie_components(G, cc, &nc)) return 0;
    if (nc == 0) return 1;                       /* the class is a singleton */

    int *pos = (int *) R_alloc((size_t) p, sizeof(int));
    memset(pos, 0, (size_t) p * sizeof(int));
    int ord[64];
    if (rng) GetRNGstate();
    for (int c = 0; c < nc; c++) {
        cp_ctx ctx; memset(&ctx, 0, sizeof(ctx)); ctx.ok = 1;
        cp_vset uni = (cc[c].m == 64) ? ~(cp_vset) 0 : (((cp_vset) 1 << cc[c].m) - 1);
        int id = cp_build(&ctx, cc[c].adj, cc[c].m, uni);
        if (id < 0) { if (rng) PutRNGstate(); return 0; }
        int tick = 0;
        cp_sample(&ctx, id, ord, &tick);
        for (int i = 0; i < tick; i++) pos[cc[c].v[ord[i]]] = i;
    }
    if (rng) PutRNGstate();

    int cap = G->na > 0 ? G->na : 1;
    int *cf = (int *) R_alloc((size_t) cap, sizeof(int));
    int *ct = (int *) R_alloc((size_t) cap, sizeof(int));
    int nch = 0;
    for (int e = 0; e < G->na; e++) {
        if (!G->undir[e]) continue;
        int u = G->af[e], v = G->at[e];
        int want_u_v = (pos[u] < pos[v]);
        int have_u_v = idl_dag_has_edge(d, u, v);
        if (want_u_v == have_u_v) continue;
        if (have_u_v) { cf[nch] = u; ct[nch] = v; }
        else          { cf[nch] = v; ct[nch] = u; }
        nch++;
    }
    for (int e = 0; e < nch; e++) idl_dag_remove_edge(d, cf[e], ct[e]);
    int left = nch;
    while (left > 0) {
        int progress = 0;
        for (int e = 0; e < nch; e++) {
            if (cf[e] < 0) continue;
            if (idl_dag_can_add(d, ct[e], cf[e])) {
                idl_dag_add_edge(d, ct[e], cf[e]);
                cf[e] = -1; left--; progress = 1;
            }
        }
        if (!progress) error("sample_and_apply: could not re-orient acyclically");
    }
    return 1;
}

/*
 * C_dag_imec_sample
 *
 * st        externalptr  the DAG state, replaced in place by a draw that is
 *                        uniform on its I-equivalence class
 * tgt_R     VECSXP       the target family, a list of 1-based integer vectors
 * maxsz_R   INTSXP       largest chain component that will be enumerated
 *
 * Returns TRUE when it sampled, FALSE when some component was too large and
 * the caller must fall back to C_dag_rcar. Consumes exactly one R_unif_index()
 * per component, in component order, matching imec.sample() in R/imec.R
 * (rule 4). Nothing that can longjmp happens between Get- and PutRNGstate().
 */
SEXP
C_dag_imec_sample(SEXP st, SEXP tgt_R, SEXP maxsz_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    (void) maxsz_R;                 /* Clique-Picking has no component limit */
    void *vmax = vmaxget();
    iess G;
    ie_build(d, tgt_R, &G);
    int ok = sample_and_apply(d, &G, 1);
    vmaxset(vmax);
    return ScalarLogical(ok ? TRUE : FALSE);
}

/*
 * C_dag_imec_size -- |[D]_I| as a double, or NA when not enumerable. The count
 * is a product of per-component AMO counts and overflows an int quickly, which
 * is why it is returned as a double.
 */
SEXP
C_dag_imec_size(SEXP st, SEXP tgt_R, SEXP maxsz_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    (void) maxsz_R;
    void *vmax = vmaxget();
    iess G;
    ie_build(d, tgt_R, &G);
    chaincomp *cc = (chaincomp *) R_alloc((size_t) (d->p > 0 ? d->p : 1), sizeof(chaincomp));
    int nc = 0;
    if (!ie_components(&G, cc, &nc)) { vmaxset(vmax); return ScalarReal(NA_REAL); }
    double sz = 1.0;
    for (int c = 0; c < nc; c++) {
        cp_ctx ctx; memset(&ctx, 0, sizeof(ctx)); ctx.ok = 1;
        cp_vset uni = (cc[c].m == 64) ? ~(cp_vset) 0 : (((cp_vset) 1 << cc[c].m) - 1);
        int id = cp_build(&ctx, cc[c].adj, cc[c].m, uni);
        if (id < 0) { vmaxset(vmax); return ScalarReal(NA_REAL); }
        sz *= ctx.node[id].total;
    }
    vmaxset(vmax);
    return ScalarReal(sz);
}

/*
 * C_dag_imec_members -- every member of the class, each as a 2 x m integer
 * matrix in C_dag_edgeM() format. Returns NULL when the class is not enumerable
 * (a chain component larger than maxsz) or when it has more than maxmem
 * members. The DAG state is left exactly as it was found.
 */
SEXP
C_dag_imec_members(SEXP st, SEXP tgt_R, SEXP maxsz_R, SEXP maxmem_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int max_size = asInteger(maxsz_R);
    double max_mem = asReal(maxmem_R);
    void *vmax = vmaxget();

    char *tprot = build_tprot(tgt_R, d->p);
    char *E = (char *) R_alloc((size_t) d->p * d->p, sizeof(char));
    idl_comp *cs = NULL;
    int enumerable = 0;
    int nc = decompose(d, tprot, max_size, E, &cs, &enumerable);
    if (!enumerable) { vmaxset(vmax); return R_NilValue; }

    double total = 1.0;
    for (int c = 0; c < nc; c++) total *= (double) cs[c].n;
    /* the exhaustive escape scores one neighbourhood per member, so its cost is
       linear in the class size; refuse rather than materialise a class the
       caller has said is too big to walk */
    if (!ISNA(max_mem) && total > max_mem) { vmaxset(vmax); return R_NilValue; }
    if (total > (double) INT_MAX)
        error("C_dag_imec_members: class of %g members is too large to list", total);
    int ntot = (int) total;

    /* remember the DAG we were handed, to restore it at the end */
    int nold = d->nedges;
    int *of = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    int *ot = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    { int k = 0;
      for (int u = 0; u < d->p; u++) {
          const idl_ivec *chu = &d->ch[u];
          for (int j = 0; j < chu->n; j++) { of[k] = u; ot[k] = chu->v[j]; k++; } } }

    SEXP ans = PROTECT(allocVector(VECSXP, ntot));
    unsigned *pick = (unsigned *) R_alloc((size_t) (nc > 0 ? nc : 1), sizeof(unsigned));
    for (int c = 0; c < nc; c++) pick[c] = 0u;
    for (int t = 0; t < ntot; t++) {
        apply_pick(d, E, cs, nc, pick);
        SEXP em = PROTECT(allocMatrix(INTSXP, 2, d->nedges));
        int *a = INTEGER(em); int k = 0;
        for (int u = 0; u < d->p; u++) {
            const idl_ivec *chu = &d->ch[u];
            for (int j = 0; j < chu->n; j++) {
                a[2 * k] = u + 1; a[2 * k + 1] = chu->v[j] + 1; k++;
            }
        }
        SEXP dn = PROTECT(allocVector(VECSXP, 2));
        SEXP rn = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(rn, 0, mkChar("from"));
        SET_STRING_ELT(rn, 1, mkChar("to"));
        SET_VECTOR_ELT(dn, 0, rn);
        SET_VECTOR_ELT(dn, 1, R_NilValue);
        setAttrib(em, R_DimNamesSymbol, dn);
        SET_VECTOR_ELT(ans, t, em);
        UNPROTECT(3);
        /* odometer over the components, last one moving fastest, which is the
           order R's nested loop in imec.list() produces */
        for (int c = nc - 1; c >= 0; c--) {
            if (++pick[c] < cs[c].n) break;
            pick[c] = 0u;
        }
    }
    set_edges(d, of, ot, nold);
    UNPROTECT(1);
    vmaxset(vmax);
    return ans;
}

/*
 * C_dag_set_edges -- replace the DAG by the given arc set, 1-based, in
 * C_dag_edgeM() format. Used to walk the members of a class during the
 * exhaustive escape.
 */
SEXP
C_dag_set_edges(SEXP st, SEXP em_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    if (TYPEOF(em_R) != INTSXP)
        error("C_dag_set_edges: 'em' must be an integer matrix");
    int n = (int) (XLENGTH(em_R) / 2);
    const int *a = INTEGER(em_R);
    void *vmax = vmaxget();
    int *from = (int *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    int *to   = (int *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    for (int e = 0; e < n; e++) {
        int u = a[2 * e], v = a[2 * e + 1];
        if (u < 1 || u > d->p || v < 1 || v > d->p || u == v)
            error("C_dag_set_edges: arc %d out of range", e + 1);
        from[e] = u - 1; to[e] = v - 1;
    }
    set_edges(d, from, to, n);
    vmaxset(vmax);
    return R_NilValue;
}
