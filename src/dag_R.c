#include <R.h>
#include <Rinternals.h>
#include "dag.h"
#include "dag_R.h"
#include "nh_scores.h"          /* IDLBNS_OP_ADD / _REMOVE / _REVERSE */

/*
 * R INTERFACE TO THE C DAG
 *
 * The DAG lives behind an external pointer with a tag and a finalizer, so
 * that R owns its lifetime. Three rules that the construction order below
 * exists to satisfy:
 *
 *   - calloc the struct, wrap it, register the finalizer, and only THEN
 *     allocate the interior. An allocation failure or an error() partway
 *     through construction then still leaves a fully reclaimable object,
 *     because every interior pointer is either live or NULL.
 *   - never store a SEXP inside the malloc'd struct; the GC cannot see it.
 *     Anything SEXP-shaped belongs in the external pointer's 'prot' slot.
 *   - validate every precondition of a mutation BEFORE the first write, so
 *     an error() can never longjmp out of a half-applied move and leave the
 *     state inconsistent for the next .Call.
 */

static SEXP idl_dag_tag = NULL;

static void
idl_dag_finalizer(SEXP ext) {
    idl_dag *d = (idl_dag *) R_ExternalPtrAddr(ext);
    if (d == NULL)
        return;                    /* already finalized, or deserialized */
    idl_dag_free(d);
    R_ClearExternalPtr(ext);       /* make a double finalize a no-op */
}

/* installed from R_init_idlBNs */
void
idl_dag_R_init(void) {
    idl_dag_tag = Rf_install("idlBNs_dag");
}

/* declared in dag_R.h; not static, because src/nbhd.c uses it too */
idl_dag *
idlBNs_dag_from_extptr(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != idl_dag_tag)
        error("not an idlBNs DAG state object");
    idl_dag *d = (idl_dag *) R_ExternalPtrAddr(x);
    /* a saved and reloaded external pointer comes back with a NULL address;
       hcmc() never returns the state so this should be unreachable, but the
       check turns a segfault into an error message */
    if (d == NULL)
        error("idlBNs DAG state object is stale (finalized or deserialized)");

    return d;
}

SEXP
C_dag_new(SEXP p_R) {
    int p = asInteger(p_R);

    if (p == NA_INTEGER || p < 1 || p > IDLBNS_PMAX)
        error("C_dag_new: 'p' must be between 1 and %d", IDLBNS_PMAX);

    /* the zeroed struct first, then hand ownership to R, and only then
       allocate the interior. idl_dag_init() uses R_Calloc, which longjmps on
       failure, so it must not run until the finalizer is registered -- see
       the contract on idl_dag_alloc() in dag.h.

       The one residual exposure is R_MakeExternalPtr() itself failing
       between the two, which would lose the bare struct (a couple of hundred
       bytes, once, on an allocation failure that is about to end the
       session anyway). */
    idl_dag *d = idl_dag_alloc(p);
    SEXP ext = PROTECT(R_MakeExternalPtr(d, idl_dag_tag, R_NilValue));
    R_RegisterCFinalizerEx(ext, idl_dag_finalizer, TRUE);
    idl_dag_init(d);
    UNPROTECT(1);

    return ext;
}

/*
 * C_dag_apply_move
 *
 * Applies one (op, u, v) move, with u and v 1-based vertex indices, in the
 * same three flavours as apply.move() in R/search.R. Every precondition is
 * checked first, including acyclicity, so the state is never left partly
 * modified.
 */
SEXP
C_dag_apply_move(SEXP st, SEXP op_R, SEXP u_R, SEXP v_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int op = asInteger(op_R);
    int u = asInteger(u_R);
    int v = asInteger(v_R);

    if (op == NA_INTEGER || u == NA_INTEGER || v == NA_INTEGER)
        error("C_dag_apply_move: 'op', 'u' and 'v' must not be NA");
    if (u < 1 || u > d->p || v < 1 || v > d->p)
        error("C_dag_apply_move: vertex out of range (u=%d, v=%d, p=%d)",
              u, v, d->p);
    if (u == v)
        error("C_dag_apply_move: self-loops are not allowed (u = v = %d)", u);

    int u0 = u - 1, v0 = v - 1;

    switch (op) {
    case IDLBNS_OP_ADD:
        if (idl_dag_has_edge(d, u0, v0))
            error("C_dag_apply_move: arc %d -> %d already present", u, v);
        if (!idl_dag_can_add(d, u0, v0))
            error("C_dag_apply_move: adding %d -> %d would close a cycle",
                  u, v);
        idl_dag_add_edge(d, u0, v0);
        break;
    case IDLBNS_OP_REMOVE:
        if (!idl_dag_has_edge(d, u0, v0))
            error("C_dag_apply_move: arc %d -> %d not present", u, v);
        idl_dag_remove_edge(d, u0, v0);
        break;
    case IDLBNS_OP_REVERSE:
        if (!idl_dag_has_edge(d, u0, v0))
            error("C_dag_apply_move: arc %d -> %d not present", u, v);
        if (!idl_dag_can_reverse(d, u0, v0))
            error("C_dag_apply_move: reversing %d -> %d would close a cycle",
                  u, v);
        idl_dag_reverse_edge(d, u0, v0);
        break;
    default:
        error("C_dag_apply_move: unknown move operation code %d", op);
    }

    return R_NilValue;
}

/* the parent sets as a list of p integer vectors of 1-based indices, in
   INSERTION order -- the same object shape .build_pasets() returns, and the
   order the score functions must see */
SEXP
C_dag_pasets(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    SEXP ans = PROTECT(allocVector(VECSXP, d->p));
    for (int v = 0; v < d->p; v++) {
        const idl_ivec *pav = &d->pa[v];
        SEXP e = PROTECT(allocVector(INTSXP, pav->n));
        int *ei = INTEGER(e);
        for (int j = 0; j < pav->n; j++)
            ei[j] = pav->v[j] + 1;
        SET_VECTOR_ELT(ans, v, e);
        UNPROTECT(1);
    }
    UNPROTECT(1);

    return ans;
}

/* the same parent sets ascending -- the canonical cache-key form */
SEXP
C_dag_pasets_sorted(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    SEXP ans = PROTECT(allocVector(VECSXP, d->p));
    for (int v = 0; v < d->p; v++) {
        const idl_ivec *pv = &d->pas[v];
        SEXP e = PROTECT(allocVector(INTSXP, pv->n));
        int *ei = INTEGER(e);
        for (int j = 0; j < pv->n; j++)
            ei[j] = pv->v[j] + 1;
        SET_VECTOR_ELT(ans, v, e);
        UNPROTECT(1);
    }
    UNPROTECT(1);

    return ans;
}

/*
 * C_dag_edgeM
 *
 * The arcs as a 2 x m integer matrix with rows "from" and "to", 1-based, in
 * graph::edgeMatrix() column order: from-vertex ascending, and within a
 * vertex its children in edgeL insertion order. rcar() selects a covered
 * arc by index into this order, so reproducing it is a correctness
 * requirement, not a convenience.
 */
SEXP
C_dag_edgeM(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    SEXP ans = PROTECT(allocMatrix(INTSXP, 2, d->nedges));
    int *a = INTEGER(ans);
    int k = 0;
    for (int u = 0; u < d->p; u++) {
        const idl_ivec *chu = &d->ch[u];
        for (int j = 0; j < chu->n; j++) {
            a[2 * k]     = u + 1;
            a[2 * k + 1] = chu->v[j] + 1;
            k++;
        }
    }
    SEXP dn = PROTECT(allocVector(VECSXP, 2));
    SEXP rn = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(rn, 0, mkChar("from"));
    SET_STRING_ELT(rn, 1, mkChar("to"));
    SET_VECTOR_ELT(dn, 0, rn);
    SET_VECTOR_ELT(dn, 1, R_NilValue);
    setAttrib(ans, R_DimNamesSymbol, dn);
    UNPROTECT(3);

    return ans;
}

/*
 * The per-vertex parent-set version stamps, so a caller can memoise things
 * that depend only on a vertex's parent set and detect staleness with one
 * integer compare. Used by the addition memo in nh_scores.c.
 */
SEXP
C_dag_pastamp(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    SEXP ans = PROTECT(allocVector(INTSXP, d->p));
    memcpy(INTEGER(ans), d->pav_stamp, (size_t) d->p * sizeof(int));
    UNPROTECT(1);

    return ans;
}

SEXP
C_dag_nedges(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);

    return ScalarInteger(d->nedges);
}

/* the ancestor matrix as R holds it: anc[u, v] is TRUE iff u is an ancestor
   of v, so column v of the result is the ancestor set of v */
SEXP
C_dag_anc(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int p = d->p;
    SEXP ans = PROTECT(allocMatrix(LGLSXP, p, p));
    int *a = LOGICAL(ans);
    for (int v = 0; v < p; v++) {
        const idl_word *anc_v = d->anc + (size_t) v * d->W;
        for (int u = 0; u < p; u++)
            a[u + (size_t) v * p] = idl_bs_test(anc_v, u);
    }
    UNPROTECT(1);

    return ans;
}

/*
 * The descendant matrix in its natural reading: element (u, v) is TRUE iff
 * v is a descendant of u.
 *
 * Under R's convention that anc[u, v] means "u is an ancestor of v", that is
 * the SAME relation, so a consistent DAG returns the same logical matrix
 * here as C_dag_anc() does. The two are stored differently -- anc's rows are
 * ancestor sets, desc's rows are descendant sets -- so comparing the two
 * exports from R is exactly the transpose invariant that idl_dag_check()
 * verifies internally, and is why this accessor is worth having rather than
 * being redundant.
 */
SEXP
C_dag_desc(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int p = d->p;
    SEXP ans = PROTECT(allocMatrix(LGLSXP, p, p));
    int *a = LOGICAL(ans);
    for (int u = 0; u < p; u++) {
        const idl_word *desc_u = d->desc + (size_t) u * d->W;
        for (int v = 0; v < p; v++)
            a[u + (size_t) v * p] = idl_bs_test(desc_u, v);
    }
    UNPROTECT(1);

    return ans;
}

/* full internal self-consistency pass; errors with a description on any
   inconsistency, returns NULL otherwise */
SEXP
C_dag_check(SEXP st) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    char msg[512];

    if (idl_dag_check(d, msg, sizeof(msg)) != 0)
        error("C_dag_check: %s", msg);

    return R_NilValue;
}
