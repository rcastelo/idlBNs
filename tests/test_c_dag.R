## 2026-09-07 differential test for the C DAG state against the R
## structures it replaces.
##
## THE SHADOW HARNESS. A random trajectory of legal moves is applied to both
## the C DAG and the R state (a graphNEL, the ancestor matrix, and the parent
## sets), and after EVERY move the two are compared exactly -- identical(),
## never a tolerance, and for the parent sets identical() rather than a
## set-wise comparison, because parent ORDER reaches the score functions and
## a permuted Cholesky rounds differently.
##
## The trajectory deliberately includes removals and reversals, which is what
## makes a vertex's child list stop being ascending. That matters: nr.nh()
## emits removals in child-insertion order and edgeMatrix() enumerates arcs
## in it, so a C implementation that stored children ascending would agree
## with R on every edge SET while changing every trajectory. The harness
## asserts below that non-ascending child lists actually occurred, so this
## cannot pass vacuously.

suppressPackageStartupMessages({
  library(graph)
  library(idlBNs)
})

## the registered routines are native symbol objects, so they go through
## .Call(); these thin closures keep the test readable
C_dag_new <- function(p)
    .Call(idlBNs:::C_dag_new, as.integer(p))
C_dag_apply_move <- function(st, op, u, v)
    .Call(idlBNs:::C_dag_apply_move, st, as.integer(op), as.integer(u),
          as.integer(v))
C_dag_pasets        <- function(st) .Call(idlBNs:::C_dag_pasets, st)
C_dag_pasets_sorted <- function(st) .Call(idlBNs:::C_dag_pasets_sorted, st)
C_dag_edgeM         <- function(st) .Call(idlBNs:::C_dag_edgeM, st)
C_dag_nedges        <- function(st) .Call(idlBNs:::C_dag_nedges, st)
C_dag_anc           <- function(st) .Call(idlBNs:::C_dag_anc, st)
C_dag_desc          <- function(st) .Call(idlBNs:::C_dag_desc, st)
C_dag_check         <- function(st) invisible(.Call(idlBNs:::C_dag_check, st))

## how many vertices ended up with a non-ascending child list, across the
## whole sweep -- the non-vacuity counter
nonasc <- 0L
nsteps_done <- 0L

## compare every exported view of the C state against the R state
compare_states <- function(st, dag, anc, pas, vn) {
  ## 1. the arc list, in edgeMatrix() column order
  stopifnot(identical(unname(C_dag_edgeM(st)), unname(edgeMatrix(dag))))
  ## 2. the parent sets, in INSERTION order
  stopifnot(identical(C_dag_pasets(st), pas))
  ## 3. the ascending mirror really is the sorted parent sets
  stopifnot(identical(C_dag_pasets_sorted(st),
                      lapply(pas, function(x) sort.int(x))))
  ## 4. the ancestor matrix, in full
  stopifnot(identical(C_dag_anc(st), unname(anc)))
  ## 5. desc is the transpose of anc (same matrix under R's convention)
  stopifnot(identical(C_dag_desc(st), C_dag_anc(st)))
  ## 6. the edge count
  stopifnot(identical(C_dag_nedges(st), as.integer(numEdges(dag))))
  ## 7. and the C structure is internally self-consistent
  C_dag_check(st)
  invisible(NULL)
}

shadow_run <- function(p, nsteps, seed) {
  set.seed(seed)
  vn <- paste0("X", seq_len(p))
  dag <- new("graphNEL", nodes=vn, edgemode="directed")
  anc <- idlBNs:::init.ancestors(vn)
  pas <- idlBNs:::init.pasets(p)
  st <- C_dag_new(p)

  compare_states(st, dag, anc, pas, vn)      ## the edgeless DAG

  for (k in seq_len(nsteps)) {
    ne <- idlBNs:::ar.nh(dag, anc)           ## every legal move
    stopifnot(length(ne$op) > 0L)
    m <- sample.int(length(ne$op), 1L)
    op <- ne$op[m]
    u <- ne$u[m]
    v <- ne$v[m]

    ## R state: anc and pasets BEFORE dag, since remove/reverse.ancestors
    ## read the DAG as it stood before the move
    anc <- switch(op,
                  idlBNs:::add.ancestors(anc, vn[u], vn[v]),
                  idlBNs:::remove.ancestors(anc, dag, vn[u], vn[v]),
                  idlBNs:::reverse.ancestors(anc, dag, vn[u], vn[v]))
    pas <- idlBNs:::move.pasets(pas, op, u, v)
    dag <- idlBNs:::apply.move(dag, op, u, v, vn)

    ## C state
    C_dag_apply_move(st, op, u, v)

    compare_states(st, dag, anc, pas, vn)
    nsteps_done <<- nsteps_done + 1L

    ## non-vacuity bookkeeping
    el <- edgeL(dag)
    for (i in seq_len(p)) {
      e <- el[[i]]$edges
      if (!is.null(e) && length(e) > 1L && !identical(as.integer(e),
                                                      sort(as.integer(e))))
        nonasc <<- nonasc + 1L
    }
  }
  invisible(NULL)
}

################################################################################
## 1. the sweep
################################################################################

for (p in c(4, 6, 8, 12, 20))
  for (nsteps in c(0, 1, 5, 40, 120))
    for (seed in 1:8)
      shadow_run(p, nsteps, seed)

stopifnot(nsteps_done > 1000L)
## if the trajectories never produced a non-ascending child list, the
## ordering half of this test would be vacuous -- see the header
stopifnot(nonasc > 0L)
cat(sprintf("shadow harness: %d moves compared, %d non-ascending child lists seen\n",
            nsteps_done, nonasc))

################################################################################
## 2. hand-built removal topologies, the cases a random sweep may under-sample
################################################################################

mkboth <- function(p, arcs) {
  vn <- paste0("X", seq_len(p))
  dag <- new("graphNEL", nodes=vn, edgemode="directed")
  st <- C_dag_new(p)
  for (a in arcs) {
    dag <- addEdge(vn[a[1]], vn[a[2]], dag)
    C_dag_apply_move(st, 1L, a[1], a[2])
  }
  anc <- idlBNs:::init.ancestors(vn)
  em <- edgeMatrix(dag)
  for (k in seq_len(ncol(em)))
    anc <- idlBNs:::add.ancestors(anc, vn[em["from", k]], vn[em["to", k]])
  list(vn=vn, dag=dag, st=st, anc=anc)
}

cases <- list(
  list(p=4, arcs=list(c(1,2), c(2,3), c(3,4)), rm=c(2,3)),           ## chain
  list(p=4, arcs=list(c(1,2), c(1,3), c(2,4), c(3,4)), rm=c(2,4)),   ## diamond side
  list(p=4, arcs=list(c(1,2), c(1,3), c(2,4), c(3,4)), rm=c(1,2)),   ## into apex
  list(p=5, arcs=list(c(1,5), c(2,5), c(3,5), c(4,5)), rm=c(2,5)),   ## many parents
  list(p=6, arcs=list(c(1,2), c(2,3), c(3,4), c(3,5), c(5,6)), rm=c(2,3)),
  list(p=2, arcs=list(c(1,2)), rm=c(1,2)),                           ## single arc
  list(p=6, arcs=list(c(1,2), c(2,3), c(4,5), c(5,6)), rm=c(1,2)),   ## two components
  list(p=4, arcs=list(c(1,2), c(1,4), c(2,3), c(3,4)), rm=c(1,4)),   ## reconvergent
  ## arcs inserted OUT of ascending order, so child lists are not sorted
  list(p=5, arcs=list(c(1,4), c(1,2), c(1,3), c(2,5)), rm=c(1,2)),
  list(p=5, arcs=list(c(1,5), c(1,3), c(1,2), c(3,4)), rm=c(1,3))
)

for (cs in cases) {
  b <- mkboth(cs$p, cs$arcs)
  u <- cs$rm[1]; v <- cs$rm[2]
  ## removal
  C_dag_apply_move(b$st, 2L, u, v)
  want <- idlBNs:::remove.ancestors(b$anc, b$dag, b$vn[u], b$vn[v])
  stopifnot(identical(C_dag_anc(b$st), unname(want)))
  stopifnot(identical(unname(C_dag_edgeM(b$st)),
                      unname(edgeMatrix(removeEdge(b$vn[u], b$vn[v], b$dag)))))
  C_dag_check(b$st)

  ## reversal, where it stays acyclic
  b2 <- mkboth(cs$p, cs$arcs)
  revg <- addEdge(b2$vn[v], b2$vn[u], removeEdge(b2$vn[u], b2$vn[v], b2$dag))
  ## decide acyclicity independently, by closure on the reversed graph
  A <- matrix(FALSE, cs$p, cs$p)
  emr <- edgeMatrix(revg)
  if (ncol(emr) > 0) A[cbind(emr["from", ], emr["to", ])] <- TRUE
  R <- A
  repeat { N <- R | ((R %*% A) > 0); if (identical(N, R)) break; R <- N }
  if (!any(diag(R))) {
    C_dag_apply_move(b2$st, 3L, u, v)
    want2 <- idlBNs:::reverse.ancestors(b2$anc, b2$dag, b2$vn[u], b2$vn[v])
    stopifnot(identical(C_dag_anc(b2$st), unname(want2)))
    stopifnot(identical(unname(C_dag_edgeM(b2$st)), unname(edgeMatrix(revg))))
    C_dag_check(b2$st)
  }
}
cat(sprintf("hand-built topologies: %d removal and reversal cases agree\n",
            length(cases)))

################################################################################
## 3. legality: an illegal move must error and leave the state untouched
################################################################################

b <- mkboth(4, list(c(1,2), c(2,3)))
before <- list(e=C_dag_edgeM(b$st), a=C_dag_anc(b$st), p=C_dag_pasets(b$st))
errs <- function(expr) inherits(tryCatch(expr, error=function(e) e), "error")

stopifnot(errs(C_dag_apply_move(b$st, 1L, 1L, 2L)))   ## arc already present
stopifnot(errs(C_dag_apply_move(b$st, 1L, 3L, 1L)))   ## would close a cycle
stopifnot(errs(C_dag_apply_move(b$st, 2L, 1L, 3L)))   ## arc not present
stopifnot(errs(C_dag_apply_move(b$st, 3L, 1L, 3L)))   ## arc not present
stopifnot(errs(C_dag_apply_move(b$st, 1L, 1L, 1L)))   ## self loop
stopifnot(errs(C_dag_apply_move(b$st, 1L, 0L, 2L)))   ## out of range
stopifnot(errs(C_dag_apply_move(b$st, 1L, 1L, 9L)))   ## out of range
stopifnot(errs(C_dag_apply_move(b$st, 9L, 1L, 2L)))   ## unknown op
stopifnot(errs(.Call(idlBNs:::C_dag_apply_move, b$st, NA_integer_, 1L, 2L)))

## nothing moved
stopifnot(identical(C_dag_edgeM(b$st), before$e),
          identical(C_dag_anc(b$st), before$a),
          identical(C_dag_pasets(b$st), before$p))
C_dag_check(b$st)

## a reversal that would close a cycle: 1->2, 1->3, 3->2; reversing 1->2
## makes 2->1 while 1 -> 3 -> 2 still stands
b <- mkboth(3, list(c(1,2), c(1,3), c(3,2)))
stopifnot(errs(C_dag_apply_move(b$st, 3L, 1L, 2L)))
C_dag_check(b$st)
cat("legality: illegal moves error and leave the state untouched\n")

################################################################################
## 4. lifecycle: a stale (serialized/deserialized) pointer must error, not
## crash, and a non-DAG SEXP must be rejected
################################################################################

b <- mkboth(4, list(c(1,2)))
stale <- unserialize(serialize(b$st, NULL))
stopifnot(errs(C_dag_pasets(stale)))
stopifnot(errs(C_dag_apply_move(stale, 1L, 2L, 3L)))
stopifnot(errs(C_dag_anc(stale)))
stopifnot(errs(C_dag_check(stale)))
stopifnot(errs(C_dag_pasets(1L)))
stopifnot(errs(C_dag_pasets(list())))
## the live one is unaffected by any of that
C_dag_check(b$st)

## many short-lived DAGs, so the finalizer runs under GC pressure
for (i in 1:200) {
  s <- C_dag_new(30L)
  C_dag_apply_move(s, 1L, 1L, 2L)
  C_dag_apply_move(s, 1L, 2L, 3L)
  rm(s)
}
invisible(gc())
cat("lifecycle: stale pointers error cleanly, finalizer survives GC pressure\n")

################################################################################
## 5. p = 1 and the multi-word boundary. W = ceil(p/64), so p = 64, 65 and
## 129 are where a single-word implementation would break.
################################################################################

s1 <- C_dag_new(1L)
stopifnot(identical(C_dag_nedges(s1), 0L),
          identical(dim(C_dag_anc(s1)), c(1L, 1L)))
C_dag_check(s1)

for (p in c(63, 64, 65, 127, 128, 129)) {
  s <- C_dag_new(p)
  ## a chain across the whole vertex range forces ancestry to cross word
  ## boundaries: vertex p ends up with p-1 ancestors
  for (i in seq_len(p - 1))
    C_dag_apply_move(s, 1L, i, i + 1L)
  a <- C_dag_anc(s)
  stopifnot(sum(a[, p]) == p - 1L)       ## every earlier vertex
  stopifnot(sum(a[1L, ]) == p - 1L)      ## vertex 1 is an ancestor of all
  stopifnot(identical(C_dag_nedges(s), as.integer(p - 1L)))
  C_dag_check(s)
  ## and unwinding it must restore the empty ancestor matrix
  for (i in seq_len(p - 1))
    C_dag_apply_move(s, 2L, i, i + 1L)
  stopifnot(!any(C_dag_anc(s)), identical(C_dag_nedges(s), 0L))
  C_dag_check(s)
}
cat("multi-word bitsets: p = 1, 63, 64, 65, 127, 128, 129 all consistent\n")

cat("all C DAG tests passed\n")
