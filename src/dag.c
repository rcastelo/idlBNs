#include <stdio.h>
#include <R.h>                  /* R_Calloc / R_Realloc / R_Free / R_alloc */
#include "dag.h"

/*
 * MEMORY
 *
 * All persistent state uses R_Calloc / R_Realloc / R_Free rather than
 * malloc / realloc / free, per Writing R Extensions 1.6.4 ("it is almost
 * always easier to use R's R_Calloc, which does check") and 6.1.2. The two
 * pools must not be mixed -- WRE is explicit that memory from R_Calloc must
 * not be released with free -- so nothing here uses the C library
 * allocators. Transient scratch inside a single call uses R_alloc, which
 * unwinds by itself.
 *
 * The consequence that shapes the code below: R_Calloc RAISES AN R ERROR on
 * failure instead of returning NULL, i.e. it longjmps. Two places have to be
 * ordered around that.
 *
 *   Construction. idl_dag_alloc() takes only the zeroed struct;
 *   idl_dag_init() then takes the interior. The caller must wrap the struct
 *   in an external pointer AND register its finalizer between the two, so a
 *   longjmp out of idl_dag_init() leaves an object R already owns and will
 *   reclaim. Every interior pointer is NULL until it is live, which is what
 *   makes the finalizer safe to run against a partially built DAG.
 *
 *   Mutation. link_edge() reserves capacity in all three lists it touches
 *   BEFORE writing anything, so a longjmp out of a reserve leaves the DAG
 *   exactly as it was. Past the reserves no allocation can fail, so a move
 *   is all-or-nothing -- which tests/test_c_dag.R asserts by checking the
 *   state is byte-identical after a rejected move.
 */

/* ------------------------------------------------------------------------ */
/* growable int vectors                                                     */
/* ------------------------------------------------------------------------ */

/* R_Calloc raises an R error on failure, so on return x->v is live */
static void
ivec_init(idl_ivec *x, int cap) {
    x->n = 0;
    x->cap = cap;
    x->v = R_Calloc((size_t) cap, int);
}

/* make room for at least 'want' elements. This is the ONLY place that can
   longjmp during a mutation, which is why link_edge() calls it for every
   list it is about to touch before touching any of them. */
static void
ivec_reserve(idl_ivec *x, int want) {
    if (want <= x->cap)
        return;
    int cap = x->cap < 4 ? 4 : x->cap;
    while (cap < want)
        cap *= 2;
    x->v = R_Realloc(x->v, (size_t) cap, int);
    x->cap = cap;
}

/* append, preserving insertion order (add.pasets / graphNEL addEdge).
   the caller must have reserved capacity, so this cannot fail. */
static void
ivec_push(idl_ivec *x, int e) {
    x->v[x->n++] = e;
}

/* delete one element, preserving the order of the survivors
   (remove.pasets' setdiff / graphNEL's order-preserving edge delete) */
static void
ivec_erase(idl_ivec *x, int e) {
    int k = 0;
    while (k < x->n && x->v[k] != e)
        k++;
    if (k == x->n)
        return;
    for (int j = k + 1; j < x->n; j++)
        x->v[j - 1] = x->v[j];
    x->n--;
}

/* insert keeping the vector ascending: binary search for the position,
   then one memmove. O(k), which is what lets the sorted mirror be
   maintained without a qsort per cache lookup. */
static void
ivec_insert_sorted(idl_ivec *x, int e) {
    int lo = 0, hi = x->n;
    while (lo < hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (x->v[mid] < e)
            lo = mid + 1;
        else
            hi = mid;
    }
    for (int j = x->n; j > lo; j--)
        x->v[j] = x->v[j - 1];
    x->v[lo] = e;
    x->n++;
}

/* ------------------------------------------------------------------------ */
/* lifecycle                                                                */
/* ------------------------------------------------------------------------ */

void
idl_dag_free(idl_dag *d) {
    if (d == NULL)
        return;
    /* every pointer is either a live R_Calloc allocation or NULL, because
       idl_dag_alloc() zeroes the struct before anything else is taken. so
       this is safe against a DAG whose interior was only partly built when
       an R_Calloc longjmped out of idl_dag_init(). R_Free tolerates NULL
       and sets its argument to NULL. */
    if (d->pa != NULL)
        for (int i = 0; i < d->p; i++) R_Free(d->pa[i].v);
    if (d->pas != NULL)
        for (int i = 0; i < d->p; i++) R_Free(d->pas[i].v);
    if (d->ch != NULL)
        for (int i = 0; i < d->p; i++) R_Free(d->ch[i].v);
    R_Free(d->pa); R_Free(d->pas); R_Free(d->ch);
    R_Free(d->pab); R_Free(d->adj); R_Free(d->anc); R_Free(d->desc);
    R_Free(d->sc_delta); R_Free(d->sc_set); R_Free(d->sc_ancsave);
    R_Free(d->sc_a); R_Free(d->sc_b); R_Free(d->sc_order);
    R_Free(d->sc_indeg); R_Free(d->sc_queue); R_Free(d->sc_inD);
    R_Free(d->pav_stamp);
    R_Free(d);
}

idl_dag *
idl_dag_alloc(int p) {
    if (p < 1 || p > IDLBNS_PMAX)
        error("idl_dag_alloc: 'p' must be between 1 and %d", IDLBNS_PMAX);

    /* the zeroed struct and nothing else. the caller wraps this in an
       external pointer and registers the finalizer before calling
       idl_dag_init(), so that the interior allocations there are already
       covered if one of them longjmps. */
    idl_dag *d = R_Calloc(1, idl_dag);
    d->p = p;
    d->W = IDL_NW(p);
    d->nedges = 0;

    return d;
}

void
idl_dag_init(idl_dag *d) {
    int p = d->p;
    size_t pw = (size_t) p * d->W;

    d->pa   = R_Calloc((size_t) p, idl_ivec);
    d->pas  = R_Calloc((size_t) p, idl_ivec);
    d->ch   = R_Calloc((size_t) p, idl_ivec);
    d->pab  = R_Calloc(pw, idl_word);
    d->adj  = R_Calloc(pw, idl_word);
    d->anc  = R_Calloc(pw, idl_word);
    d->desc = R_Calloc(pw, idl_word);
    d->sc_delta   = R_Calloc(d->W, idl_word);
    d->sc_set     = R_Calloc(d->W, idl_word);
    d->sc_ancsave = R_Calloc(pw, idl_word);
    d->sc_a     = R_Calloc((size_t) p, int);
    d->sc_b     = R_Calloc((size_t) p, int);
    d->sc_order = R_Calloc((size_t) p, int);
    d->sc_indeg = R_Calloc((size_t) p, int);
    d->sc_queue = R_Calloc((size_t) p, int);
    d->sc_inD   = R_Calloc((size_t) p, int);
    d->pav_stamp = R_Calloc((size_t) p, int);
    for (int i = 0; i < p; i++)
        d->pav_stamp[i] = 1;         /* 0 is the "never computed" sentinel */

    for (int i = 0; i < p; i++) {
        ivec_init(&d->pa[i], 4);
        ivec_init(&d->pas[i], 4);
        ivec_init(&d->ch[i], 4);
    }
}

/* ------------------------------------------------------------------------ */
/* legality                                                                 */
/* ------------------------------------------------------------------------ */

/* adding u -> v closes a cycle iff v is already an ancestor of u. this is
   the same O(1) test nr.nh() makes as `!anc[na, i]` */
int
idl_dag_can_add(const idl_dag *d, int u, int v) {
    if (u == v || idl_dag_has_edge(d, u, v))
        return 0;

    return !idl_dag_is_anc(d, v, u);
}

/* reversing u -> v closes a cycle iff some OTHER child of u is an ancestor
   of v, i.e. there is a path u -> ... -> v that avoids the direct edge.
   this is .reversible()'s any(anc[a[-j], a[j]]) */
int
idl_dag_can_reverse(const idl_dag *d, int u, int v) {
    if (!idl_dag_has_edge(d, u, v))
        return 0;
    const idl_word *anc_v = d->anc + (size_t) v * d->W;
    const idl_ivec *chu = &d->ch[u];
    for (int k = 0; k < chu->n; k++) {
        int w = chu->v[k];
        if (w != v && idl_bs_test(anc_v, w))
            return 0;
    }

    return 1;
}

/* ------------------------------------------------------------------------ */
/* mutations                                                                */
/* ------------------------------------------------------------------------ */

/* the adjacency half of adding u -> v, shared by add and reverse */
static void
link_edge(idl_dag *d, int u, int v) {
    /* reserve every list first: R_Realloc longjmps on failure, and doing
       this up front is what keeps a move all-or-nothing. past this point
       nothing below can fail. */
    ivec_reserve(&d->pa[v], d->pa[v].n + 1);
    ivec_reserve(&d->pas[v], d->pas[v].n + 1);
    ivec_reserve(&d->ch[u], d->ch[u].n + 1);

    ivec_push(&d->pa[v], u);              /* INSERTION order  */
    ivec_insert_sorted(&d->pas[v], u);    /* ascending mirror */
    d->pav_stamp[v]++;                    /* pa(v) changed    */
    ivec_push(&d->ch[u], v);              /* edgeL order      */
    idl_bs_set(d->pab + (size_t) v * d->W, u);
    idl_bs_set(d->adj + (size_t) u * d->W, v);
    d->nedges++;
}

static void
unlink_edge(idl_dag *d, int u, int v) {
    ivec_erase(&d->pa[v], u);
    ivec_erase(&d->pas[v], u);
    d->pav_stamp[v]++;                    /* pa(v) changed    */
    ivec_erase(&d->ch[u], v);
    idl_bs_clear(d->pab + (size_t) v * d->W, u);
    idl_bs_clear(d->adj + (size_t) u * d->W, v);
    d->nedges--;
}

/*
 * idl_dag_add_edge -- a direct port of add.ancestors() plus add.pasets().
 *
 *   delta = anc(u) | {u}            the ancestry u contributes
 *   D     = desc(v) | {v}           everything that gains it
 *   for d in D:  anc[d] |= delta
 *   for a in delta: desc[a] |= D    the transpose, kept in lockstep
 *
 * anc[u, v] <- TRUE falls out of this, since u is in delta and v is in D.
 * Cost O((|delta| + |D|) * W).
 */
void
idl_dag_add_edge(idl_dag *d, int u, int v) {
    size_t W = d->W;

    link_edge(d, u, v);

    idl_bs_copy(d->sc_delta, d->anc + (size_t) u * W, W);
    idl_bs_set(d->sc_delta, u);
    idl_bs_copy(d->sc_set, d->desc + (size_t) v * W, W);
    idl_bs_set(d->sc_set, v);

    int nD = idl_bs_to_ints(d->sc_set, W, d->sc_a);
    for (int k = 0; k < nD; k++)
        idl_bs_or(d->anc + (size_t) d->sc_a[k] * W, d->sc_delta, W);

    int nA = idl_bs_to_ints(d->sc_delta, W, d->sc_b);
    for (int k = 0; k < nA; k++)
        idl_bs_or(d->desc + (size_t) d->sc_b[k] * W, d->sc_set, W);
}

/*
 * idl_dag_remove_edge -- a port of remove.ancestors() that needs no
 * external topological sort, and so drops the RBGL dependency.
 *
 * The R version recomputes the ancestor column of every vertex in
 * D = {v} u desc(v), walking D in a topological order of the whole graph so
 * that each vertex's parents are already correct when it is reached. Two
 * facts make a local Kahn sort over the induced subgraph D sufficient:
 *
 *   - D is descendant-closed: if k is in D and k -> m, then m is in D.
 *   - pa(v) n D = empty, because a parent of v cannot also be a descendant
 *     of v in a DAG. So after the edit v is the unique in-degree-0 vertex of
 *     D, and every other k in D has a parent in D (its predecessor on some
 *     path from v).
 *
 * And D may be read off the PRE-removal desc[v], because removing an
 * INCOMING edge of v does not change v's descendants.
 *
 * The result is independent of which valid topological order is used --
 * asserted directly in tests/test_search_internals.R over 213 orders --
 * which is what makes swapping RBGL::tsort for this provably safe.
 */
void
idl_dag_remove_edge(idl_dag *d, int u, int v) {
    size_t W = d->W;
    int p = d->p;

    /* D from the pre-removal descendants of v */
    idl_bs_copy(d->sc_set, d->desc + (size_t) v * W, W);
    idl_bs_set(d->sc_set, v);
    int nD = idl_bs_to_ints(d->sc_set, W, d->sc_a);

    unlink_edge(d, u, v);

    for (int k = 0; k < nD; k++)
        d->sc_inD[d->sc_a[k]] = 1;

    /* Kahn over the subgraph induced on D */
    for (int k = 0; k < nD; k++) {
        int x = d->sc_a[k];
        int deg = 0;
        const idl_ivec *pax = &d->pa[x];
        for (int j = 0; j < pax->n; j++)
            if (d->sc_inD[pax->v[j]])
                deg++;
        d->sc_indeg[x] = deg;
    }
    int qh = 0, qt = 0, no = 0;
    for (int k = 0; k < nD; k++)
        if (d->sc_indeg[d->sc_a[k]] == 0)
            d->sc_queue[qt++] = d->sc_a[k];
    while (qh < qt) {
        int x = d->sc_queue[qh++];
        d->sc_order[no++] = x;
        const idl_ivec *chx = &d->ch[x];
        for (int j = 0; j < chx->n; j++) {
            int c = chx->v[j];
            if (d->sc_inD[c] && --d->sc_indeg[c] == 0)
                d->sc_queue[qt++] = c;
        }
    }
    /* no == nD unless D contained a cycle, which cannot happen in a DAG;
       if it somehow did, fall back to the ascending order so that the
       recomputation below is still well defined rather than partial */
    if (no != nD)
        for (no = 0; no < nD; no++)
            d->sc_order[no] = d->sc_a[no];

    /* save the old ancestor rows of D, then recompute them in topo order */
    for (int k = 0; k < nD; k++) {
        int x = d->sc_a[k];
        idl_bs_copy(d->sc_ancsave + (size_t) x * W, d->anc + (size_t) x * W, W);
    }
    for (int k = 0; k < nD; k++) {
        int x = d->sc_order[k];
        idl_word *anc_x = d->anc + (size_t) x * W;
        idl_bs_zero(anc_x, W);
        const idl_ivec *pax = &d->pa[x];
        for (int j = 0; j < pax->n; j++) {
            int a = pax->v[j];
            idl_bs_set(anc_x, a);
            idl_bs_or(anc_x, d->anc + (size_t) a * W, W);
        }
    }

    /* repair the transpose. removing an edge can only LOSE ancestry, never
       gain it, so it is enough to clear the bits that went away. */
    for (int k = 0; k < nD; k++) {
        int x = d->sc_a[k];
        const idl_word *old = d->sc_ancsave + (size_t) x * W;
        const idl_word *now = d->anc + (size_t) x * W;
        for (size_t w = 0; w < W; w++) {
            idl_word lost = old[w] & ~now[w];
            while (lost) {
                int a = (int) (w * IDL_WBITS) + idl_ctz64(lost);
                lost &= lost - 1;
                idl_bs_clear(d->desc + (size_t) a * W, x);
            }
        }
    }

    for (int k = 0; k < nD; k++)
        d->sc_inD[d->sc_a[k]] = 0;
    (void) p;
}

/* reverse.ancestors(): remove then add, in that order */
void
idl_dag_reverse_edge(idl_dag *d, int u, int v) {
    idl_dag_remove_edge(d, u, v);
    idl_dag_add_edge(d, v, u);
}

/* ------------------------------------------------------------------------ */
/* covered arcs                                                             */
/* ------------------------------------------------------------------------ */

/*
 * idl_dag_arc_is_covered
 *
 * Is the arc i -> w covered, i.e. pa(i) == pa(w) \ {i}?
 *
 * Both parent lists are held ascending and duplicate-free (the 'pas'
 * mirror), and i is necessarily in pa(w) because the arc exists, so the two
 * sizes must differ by exactly one and the comparison is a single linear
 * merge. That makes this O(|pa(i)| + |pa(w)|) with no sorting and no
 * dependence on p -- which is the whole reason the mirror is maintained
 * incrementally rather than sorted on demand.
 *
 * Used by the NCR neighbourhood (an I-covered arc is excluded from it) and
 * by rcar() (which reverses only covered arcs).
 */
int
idl_dag_arc_is_covered(const idl_dag *d, int i, int w) {
    const idl_ivec *pi = &d->pas[i];
    const idl_ivec *pw = &d->pas[w];

    if (pi->n != pw->n - 1)
        return 0;

    int a = 0;
    for (int b = 0; b < pw->n; b++) {
        int x = pw->v[b];
        if (x == i)
            continue;                       /* the setdiff */
        if (a >= pi->n || pi->v[a] != x)
            return 0;
        a++;
    }

    return a == pi->n;
}

/* ------------------------------------------------------------------------ */
/* self-consistency                                                         */
/* ------------------------------------------------------------------------ */

#define CHK(cond, ...)                                                    \
    do {                                                                  \
        if (!(cond)) {                                                    \
            snprintf(msg, msglen, __VA_ARGS__);                           \
            rc = 1;                                                       \
            goto done;                                                    \
        }                                                                 \
    } while (0)

int
idl_dag_check(const idl_dag *d, char *msg, size_t msglen) {
    int p = d->p;
    size_t W = d->W;
    int rc = 0;
    /* transient scratch: R_alloc unwinds by itself, so there is nothing to
       release on the error paths below (WRE 6.1.1) */
    void *vmax = vmaxget();
    idl_word *ref = (idl_word *) R_alloc((size_t) p * W, sizeof(idl_word));
    int *order = (int *) R_alloc((size_t) p, sizeof(int));
    int *indeg = (int *) R_alloc((size_t) p, sizeof(int));
    int *queue = (int *) R_alloc((size_t) p, sizeof(int));

    memset(ref, 0, (size_t) p * W * sizeof(idl_word));
    msg[0] = '\0';

    /* 1. lists are duplicate-free and agree with their bit views */
    int esum = 0;
    for (int v = 0; v < p; v++) {
        const idl_ivec *pav = &d->pa[v];
        esum += d->ch[v].n;
        CHK(pav->n == d->pas[v].n,
            "vertex %d: |pa| = %d but |pas| = %d", v + 1, pav->n, d->pas[v].n);
        for (int j = 0; j < pav->n; j++) {
            int a = pav->v[j];
            CHK(a >= 0 && a < p && a != v,
                "vertex %d: parent %d out of range", v + 1, a + 1);
            for (int j2 = j + 1; j2 < pav->n; j2++)
                CHK(pav->v[j2] != a, "vertex %d: duplicate parent %d",
                    v + 1, a + 1);
            CHK(idl_bs_test(d->pab + (size_t) v * W, a),
                "vertex %d: parent %d missing from pab", v + 1, a + 1);
            CHK(idl_dag_has_edge(d, a, v),
                "vertex %d: parent %d not a child of it in adj", v + 1, a + 1);
        }
        /* pas is ascending and is a permutation of pa */
        for (int j = 1; j < d->pas[v].n; j++)
            CHK(d->pas[v].v[j - 1] < d->pas[v].v[j],
                "vertex %d: pas not strictly ascending at %d", v + 1, j);
        for (int j = 0; j < d->pas[v].n; j++)
            CHK(idl_bs_test(d->pab + (size_t) v * W, d->pas[v].v[j]),
                "vertex %d: pas member %d not in pa", v + 1,
                d->pas[v].v[j] + 1);
        /* children agree with the parent lists */
        const idl_ivec *chv = &d->ch[v];
        for (int j = 0; j < chv->n; j++) {
            int c = chv->v[j];
            CHK(c >= 0 && c < p && c != v,
                "vertex %d: child %d out of range", v + 1, c + 1);
            CHK(idl_bs_test(d->pab + (size_t) c * W, v),
                "vertex %d: child %d does not list it as a parent",
                v + 1, c + 1);
        }
        CHK(idl_bs_count(d->pab + (size_t) v * W, W) == pav->n,
            "vertex %d: pab popcount %d != |pa| %d", v + 1,
            idl_bs_count(d->pab + (size_t) v * W, W), pav->n);
        CHK(idl_bs_count(d->adj + (size_t) v * W, W) == chv->n,
            "vertex %d: adj popcount %d != |ch| %d", v + 1,
            idl_bs_count(d->adj + (size_t) v * W, W), chv->n);
    }
    CHK(esum == d->nedges, "nedges = %d but the child lists hold %d",
        d->nedges, esum);

    /* 2. anc equals the transitive closure of the stored adjacency */
    {
        int qh = 0, qt = 0, no = 0;
        for (int v = 0; v < p; v++) {
            indeg[v] = d->pa[v].n;
            if (indeg[v] == 0)
                queue[qt++] = v;
        }
        while (qh < qt) {
            int x = queue[qh++];
            order[no++] = x;
            const idl_ivec *chx = &d->ch[x];
            for (int j = 0; j < chx->n; j++)
                if (--indeg[chx->v[j]] == 0)
                    queue[qt++] = chx->v[j];
        }
        CHK(no == p, "the stored adjacency is cyclic (%d of %d sorted)", no, p);
        for (int k = 0; k < p; k++) {
            int x = order[k];
            idl_word *rx = ref + (size_t) x * W;
            const idl_ivec *pax = &d->pa[x];
            for (int j = 0; j < pax->n; j++) {
                int a = pax->v[j];
                idl_bs_set(rx, a);
                idl_bs_or(rx, ref + (size_t) a * W, W);
            }
        }
        for (int v = 0; v < p; v++)
            CHK(idl_bs_eq(d->anc + (size_t) v * W, ref + (size_t) v * W, W),
                "vertex %d: anc differs from the transitive closure", v + 1);
    }

    /* 3. desc is the exact transpose of anc */
    for (int v = 0; v < p; v++)
        for (int u = 0; u < p; u++)
            CHK(idl_bs_test(d->desc + (size_t) u * W, v) ==
                idl_bs_test(d->anc + (size_t) v * W, u),
                "desc[%d] and anc[%d] disagree about %d -> %d",
                u + 1, v + 1, u + 1, v + 1);

done:
    vmaxset(vmax);

    return rc;
}
