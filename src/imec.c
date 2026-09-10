#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
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
 * TWO ORDERING RULES, so that this engine and the R one stay bit-identical
 * and tests/test_search_engines.R keeps passing:
 *
 *   1. Arrows are visited in (head, tail) ascending order, matching R's
 *      which(G & !t(G), arr.ind = TRUE), which is column-major. The
 *      protected-arrow labelling updates flags in place, so the order in
 *      which arrows are relabelled within a pass is observable.
 *
 *   2. Chain components come out ordered by their smallest vertex, and each
 *      component's vertex list ascending, and are visited in that order by
 *      both the sampler and the enumerator.
 *
 * Counting, sampling and enumeration are all delegated to Clique-Picking in
 * cliquepick.c; nothing here enumerates orientations by hand any more.
 */

#define PROTECTED     1
#define UNDECIDABLE   2
#define NOT_PROTECTED 3

/*
 * Rebuild the DAG from an edge list: remove everything, then add the new arcs
 * in a topological order of the target graph, so no intermediate state can be
 * cyclic and idl_dag_add_edge's precondition always holds.
 *
 * The whole edge list is validated -- no repeated arc, and a full topological
 * order exists -- BEFORE the old graph is touched, so a rejected list leaves
 * the state exactly as it was found. Validating on the way in, as an earlier
 * version did, meant a cyclic list raised only after the old arcs had been
 * removed and the acyclic prefix inserted, leaving the external pointer
 * holding neither the old graph nor the requested one; a repeated arc was not
 * detected at all, and was inserted twice.
 *
 * Arcs are bucketed by tail, which makes both passes O(V + E) and lets the
 * insertion below visit them in exactly the order the single-pass version did:
 * vertices in the order the queue produces them, and within a vertex the order
 * the arcs appear in the input. That order is observable -- it fixes the layout
 * of each ch[] vector, hence the order neighbourhoods are enumerated and the
 * random stream consumed -- so it is preserved deliberately.
 */
static void
set_edges(idl_dag *d, const int *from, const int *to, int n) {
    int p = d->p;

    int *tstart = (int *) R_alloc((size_t) p + 1, sizeof(int));
    int *tlist  = (int *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    int *indeg  = (int *) R_alloc((size_t) p, sizeof(int));
    memset(tstart, 0, ((size_t) p + 1) * sizeof(int));
    memset(indeg, 0, (size_t) p * sizeof(int));
    for (int e = 0; e < n; e++) { tstart[from[e] + 1]++; indeg[to[e]]++; }
    for (int v = 0; v < p; v++) tstart[v + 1] += tstart[v];
    int *fill = (int *) R_alloc((size_t) p, sizeof(int));
    memcpy(fill, tstart, (size_t) p * sizeof(int));
    for (int e = 0; e < n; e++) tlist[fill[from[e]]++] = to[e];

    /* no arc given twice; sorted on a copy, so tlist keeps the input order */
    int *scratch = (int *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    for (int u = 0; u < p; u++) {
        int b = tstart[u], f = tstart[u + 1];
        if (f - b < 2) continue;
        memcpy(scratch, tlist + b, (size_t) (f - b) * sizeof(int));
        R_isort(scratch, f - b);
        for (int j = 1; j < f - b; j++)
            if (scratch[j] == scratch[j - 1])
                error("set_edges: arc %d -> %d given more than once",
                      u + 1, scratch[j] + 1);
    }

    /* acyclic, and the order the arcs will be inserted in */
    int *queue = (int *) R_alloc((size_t) p, sizeof(int));
    int qh = 0, qt = 0;
    for (int v = 0; v < p; v++) if (indeg[v] == 0) queue[qt++] = v;
    while (qh < qt) {
        int u = queue[qh++];
        for (int j = tstart[u]; j < tstart[u + 1]; j++)
            if (--indeg[tlist[j]] == 0) queue[qt++] = tlist[j];
    }
    if (qt < p) error("set_edges: target edge list is cyclic");

    /* ---- validated; only now is the old graph destroyed --------------- */
    int nold = d->nedges;
    int *of = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    int *ot = (int *) R_alloc((size_t) (nold > 0 ? nold : 1), sizeof(int));
    int k = 0;
    for (int u = 0; u < p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++) { of[k] = u; ot[k] = chu->v[j]; k++; }
    }
    for (int e = 0; e < k; e++) idl_dag_remove_edge(d, of[e], ot[e]);

    for (int i = 0; i < p; i++) {
        int u = queue[i];
        for (int j = tstart[u]; j < tstart[u + 1]; j++)
            idl_dag_add_edge(d, u, tlist[j]);
    }

    /* Inserting by topological source order reproduces each vertex's child
       sequence but not its parent insertion order, so a rebuild is not the
       identity on pa[] even when handed the edge list just read out of the
       graph. The original order cannot be recovered from an edge matrix at
       all -- it is the history of the search, while the matrix groups arcs by
       tail -- so make the rebuild canonical instead: ascending, which is what
       the R engine's rebuild from an adjacency matrix produces. A rebuild is
       then idempotent and history-independent, and the exhaustive escape's
       restore leaves the same state in both engines. */
    idl_dag_canonical_order(d);
}

/* the member of the class selected by pick[c] in each component */

/* ===================================================================== *
 * THE SPARSE PATH
 *
 * Sampling runs on every iteration of the search, so it is the one that has to
 * scale. The matrix form this replaced -- a p x p char array, enumerated over
 * all ordered pairs -- cost Theta(p^2) per call whatever the graph's density,
 * and at p = 200 that scaffolding, not Clique-Picking, was most of the time.
 *
 * This representation carries the arcs themselves plus two CSR indexes: the
 * arcs grouped by head, for the configuration scans of the labelling, and the
 * skeleton with each vertex's neighbours ascending, for adjacency queries by
 * binary search. Everything below is O(V + E log deg); nothing iterates pairs.
 *
 * The enumeration used by the exhaustive escape shares this representation:
 * it computes the class size first, refuses an over-budget class before
 * anything is materialised, and then unranks the members one at a time.
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

/*
 * How many of those a graph on p vertices can have. ie_components() keeps
 * only components of at least two vertices, and components are disjoint, so
 * the count is at most floor(p / 2); the + 1 keeps the allocation non-zero at
 * p = 1. At 528 bytes apiece this is worth halving: p entries came to 21 MB
 * at p = 40000. ie_components() rechecks the bound before it writes, so the
 * tighter allocation cannot be overrun if that invariant ever changes.
 */
#define CHAINCOMP_CAP(p) ((size_t) ((p) / 2 + 1))

static int
ie_components(const iess *G, chaincomp *cc, int *ncomp) {
    int p = G->p;
    char *seen = (char *) R_alloc((size_t) p, sizeof(char));
    int *stack = (int *) R_alloc((size_t) p, sizeof(int));
    int *loc = (int *) R_alloc((size_t) p, sizeof(int));
    /* One traversal buffer, reused by every component, with only the n
       vertices actually found copied out and kept. Giving each component its
       own p-element buffer instead cost Theta(p^2): a graph whose chain
       components are all small -- a perfect matching being the extreme, p/2
       components of two vertices -- allocated p/2 buffers of p integers, 800
       MB at p = 20000, to hold p integers' worth of vertices. Scratch is now
       Theta(p) and what is retained is sum(n) <= p. */
    int *buf = (int *) R_alloc((size_t) p, sizeof(int));
    memset(seen, 0, (size_t) p);
    int nc = 0;
    for (int v0 = 0; v0 < p; v0++) {
        if (seen[v0]) continue;
        int has = 0;
        for (int k = G->sstart[v0]; k < G->sstart[v0 + 1]; k++)
            if (G->undir[G->sarc[k]]) { has = 1; break; }
        if (!has) { seen[v0] = 1; continue; }
        int top = 0, n = 0;
        seen[v0] = 1; stack[top++] = v0; buf[n++] = v0;
        while (top > 0) {
            int w = stack[--top];
            for (int k = G->sstart[w]; k < G->sstart[w + 1]; k++) {
                if (!G->undir[G->sarc[k]]) continue;
                int z = G->slist[k];
                if (!seen[z]) { seen[z] = 1; stack[top++] = z; buf[n++] = z; }
            }
        }
        if (n < 2) continue;
        if (n > 64) return 0;
        for (int a = 1; a < n; a++) {
            int x = buf[a], b = a - 1;
            while (b >= 0 && buf[b] > x) { buf[b + 1] = buf[b]; b--; }
            buf[b + 1] = x;
        }
        /* kept: exactly the n vertices of this component, never p */
        int *vs = (int *) R_alloc((size_t) n, sizeof(int));
        memcpy(vs, buf, (size_t) n * sizeof(int));
        if ((size_t) nc >= CHAINCOMP_CAP(p))
            error("ie_components: more than %d chain components on %d "
                  "vertices", (int) CHAINCOMP_CAP(p), p);
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
    chaincomp *cc = (chaincomp *) R_alloc(CHAINCOMP_CAP(p), sizeof(chaincomp));
    int nc = 0;
    if (!ie_components(G, cc, &nc)) return 0;
    if (nc == 0) {                               /* the class is a singleton */
        /* still canonicalise: the R path rebuilds pasets from the sampled
           adjacency on every draw, this one included, so skipping it here
           left the engines in different orders exactly when there was
           nothing to sample */
        idl_dag_canonical_order(d);
        return 1;
    }

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

    /* Only the changed arcs were touched, so pa[] and ch[] are now in an
       order that depends on the search's history. The R engine rebuilds both
       from the sampled adjacency matrix, in ascending order. Match it, or the
       two engines hand the score function the same parent SETS in different
       orders -- the sets agreed in every one of 300 draws tested, the order
       differed in 183 -- and the score function's arithmetic is not
       associative, so the engines' scores drift apart in the last bits. */
    idl_dag_canonical_order(d);
    return 1;
}

/*
 * C_dag_imec_sample
 *
 * st        externalptr  the DAG state, replaced in place by a draw that is
 *                        uniform on its I-equivalence class
 * tgt_R     VECSXP       the target family, a list of 1-based integer vectors
 *
 * Returns TRUE when it sampled, FALSE when some component was too large and
 * the caller must fall back to C_dag_rcar. Consumes exactly one R_unif_index()
 * per component, in component order, matching imec.sample() in R/imec.R
 * (rule 4). Nothing that can longjmp happens between Get- and PutRNGstate().
 */
SEXP
C_dag_imec_sample(SEXP st, SEXP tgt_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
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
C_dag_imec_size(SEXP st, SEXP tgt_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    void *vmax = vmaxget();
    iess G;
    ie_build(d, tgt_R, &G);
    chaincomp *cc = (chaincomp *) R_alloc(CHAINCOMP_CAP(d->p), sizeof(chaincomp));
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
 * matrix in C_dag_edgeM() format. Returns NULL when the class has more than
 * maxmem members. The DAG state is not touched at all -- not saved and put
 * back, but never written in the first place.
 *
 * The class size is computed FIRST, by Clique-Picking and without enumerating
 * anything, and the budget is applied to it before a single member is built.
 * The previous version enumerated the m! vertex orderings of each chain
 * component and de-duplicated them, which allocated an m!-entry table before
 * the budget could be consulted: an 11-vertex path, whose class has eleven
 * members, cost 160 MB and 3.2 seconds to arrive at "too large", and a
 * 13-vertex one overflowed the 32-bit count.
 */
SEXP
C_dag_imec_members(SEXP st, SEXP tgt_R, SEXP maxmem_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    double max_mem = asReal(maxmem_R);
    void *vmax = vmaxget();

    iess G;
    ie_build(d, tgt_R, &G);
    chaincomp *cc = (chaincomp *) R_alloc(CHAINCOMP_CAP(d->p), sizeof(chaincomp));
    int nc = 0;
    if (!ie_components(&G, cc, &nc)) { vmaxset(vmax); return R_NilValue; }

    /* count first: polynomial, and no member is materialised */
    cp_ctx *ctxs = (cp_ctx *) R_alloc((size_t) (nc > 0 ? nc : 1), sizeof(cp_ctx));
    int *ids = (int *) R_alloc((size_t) (nc > 0 ? nc : 1), sizeof(int));
    double total = 1.0;
    for (int c = 0; c < nc; c++) {
        memset(&ctxs[c], 0, sizeof(cp_ctx)); ctxs[c].ok = 1;
        cp_vset uni = (cc[c].m == 64) ? ~(cp_vset) 0 : (((cp_vset) 1 << cc[c].m) - 1);
        ids[c] = cp_build(&ctxs[c], cc[c].adj, cc[c].m, uni);
        if (ids[c] < 0) { vmaxset(vmax); return R_NilValue; }
        total *= ctxs[c].node[ids[c]].total;
    }
    if (!R_FINITE(total) || (!ISNA(max_mem) && total > max_mem) ||
        total > (double) INT_MAX) {
        vmaxset(vmax); return R_NilValue;
    }
    int ntot = (int) total;

    /* Nothing below mutates the DAG: ie_build() reads it, and each member is
       emitted as an edge matrix rather than applied. An earlier version
       enumerated by orienting the state in place and captured the arcs here
       to put back afterwards; that capture-and-restore outlived its reason
       and was not free. set_edges() rebuilds by topological source order,
       which reproduces each vertex's child sequence but NOT its parent
       insertion order, so restoring perturbed pa[] -- the sequence handed to
       the score function -- in 119 of 200 states tested. An exhaustive escape
       that found no improving move still changed the arithmetic of every
       later score. */
    int *pos = (int *) R_alloc((size_t) d->p, sizeof(int));
    int ord[64];
    SEXP ans = PROTECT(allocVector(VECSXP, ntot));
    for (int t = 0; t < ntot; t++) {
        /* decode t into one member index per component */
        memset(pos, 0, (size_t) d->p * sizeof(int));
        long rest = t;
        for (int c = 0; c < nc; c++) {
            long tc = (long) (ctxs[c].node[ids[c]].total + 0.5);
            int tick = 0;
            cp_member(&ctxs[c], ids[c], (double) (rest % tc), ord, &tick);
            rest /= tc;
            if (tick != cc[c].m)
                error("C_dag_imec_members: member covered %d of %d vertices",
                      tick, cc[c].m);
            for (int i = 0; i < tick; i++) pos[cc[c].v[ord[i]]] = i;
        }
        /* the member's arcs: directed ones as they are, undirected by position */
        SEXP em = PROTECT(allocMatrix(INTSXP, 2, G.na));
        int *a = INTEGER(em);
        for (int e = 0; e < G.na; e++) {
            int u = G.af[e], v = G.at[e];
            if (G.undir[e] && pos[u] > pos[v]) { int tmp = u; u = v; v = tmp; }
            a[2 * e] = u + 1; a[2 * e + 1] = v + 1;
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
    }
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
    /* shape checked before anything else: the arc count came from
       XLENGTH()/2, which reads a matrix of any other shape as pairs */
    if (TYPEOF(em_R) != INTSXP || !isMatrix(em_R) || nrows(em_R) != 2)
        error("C_dag_set_edges: 'em' must be a 2-row integer matrix");
    int n = ncols(em_R);
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
