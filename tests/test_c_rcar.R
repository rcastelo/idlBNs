## 2026-09-07 differential test for the C port of rcar() against the R one.
##
## rcar() is the only RNG consumer in the search, so it is the only place a
## port can silently desynchronise R's random stream. The load-bearing
## assertion here is identical(.Random.seed) AFTER the call, not just
## equality of the resulting DAG: under the default sample.kind =
## "Rejection" a size-1 draw consumes a data-dependent 1 to 3 unif_rand()
## calls, so two implementations can agree on every arc they reverse while
## leaving the generator in different places -- and every later step of the
## search would then diverge.
##
## Also pinned: the two early-return branches consume NO randomness, and
## r = 0 still consumes exactly one draw (sample(0:0, 1) goes through
## sample.int(1, 1), and R_unif_index(1) computes bits = 0 but still calls
## unif_rand() once). A port that short-circuits either would pass a
## values-only comparison and break the stream.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

dag_new   <- function(p) .Call(idlBNs:::C_dag_new, as.integer(p))
dag_move  <- function(st, op, u, v)
    invisible(.Call(idlBNs:::C_dag_apply_move, st, as.integer(op),
                    as.integer(u), as.integer(v)))
dag_rcar  <- function(st, rlen, ut)
    .Call(idlBNs:::C_dag_rcar, st, as.integer(rlen), ut)
dag_edgeM <- function(st) .Call(idlBNs:::C_dag_edgeM, st)
dag_anc   <- function(st) .Call(idlBNs:::C_dag_anc, st)
dag_pas   <- function(st) .Call(idlBNs:::C_dag_pasets, st)
dag_ced   <- function(st, ut) .Call(idlBNs:::C_dag_cedges, st, ut)
dag_check <- function(st) invisible(.Call(idlBNs:::C_dag_check, st))

## build the same DAG in both representations, inserting arcs in the given
## order so the child lists match
both <- function(p, arcs) {
  vn <- paste0("X", seq_len(p))
  dag <- new("graphNEL", nodes=vn, edgemode="directed")
  st <- dag_new(p)
  for (a in arcs) {
    dag <- addEdge(vn[a[1]], vn[a[2]], dag)
    dag_move(st, 1L, a[1], a[2])
  }
  anc <- idlBNs:::init.ancestors(vn)
  em <- edgeMatrix(dag)
  for (k in seq_len(ncol(em)))
    anc <- idlBNs:::add.ancestors(anc, vn[em["from", k]], vn[em["to", k]])

  ## the parent sets must be built INCREMENTALLY, in arc insertion order,
  ## because that is what the search maintains -- pasets starts empty in
  ## hcmc()/hillclimbing() and grows through add.pasets(), which appends.
  ## Seeding them from .build_pasets() instead would give edgeMatrix order
  ## (sorted by from-vertex), a different permutation of the same sets, and
  ## parent ORDER is bit-load-bearing: it reaches the score functions and a
  ## permuted Cholesky rounds differently.
  pas <- idlBNs:::init.pasets(p)
  for (a in arcs)
    pas <- idlBNs:::add.pasets(pas, as.integer(a[1]), as.integer(a[2]))

  list(vn=vn, dag=dag, st=st, anc=anc, pas=pas,
       vidx=setNames(seq_len(p), vn))
}

## run both implementations from the same seed and compare everything
compare_rcar <- function(b, r, ut, seed) {
  set.seed(seed)
  out <- idlBNs:::rcar(b$dag, r, ut, b$anc, b$pas, b$vidx)
  seedR <- .Random.seed
  set.seed(seed)
  rr <- dag_rcar(b$st, length(0:r), ut)
  seedC <- .Random.seed

  ## the stream: same number of unif_rand() calls, same position
  stopifnot(identical(seedR, seedC))
  ## the resulting DAG, ancestor matrix and parent sets (exact order)
  stopifnot(identical(unname(dag_edgeM(b$st)), unname(edgeMatrix(out$dag))))
  stopifnot(identical(dag_anc(b$st), unname(out$anc)))
  stopifnot(identical(dag_pas(b$st), unname(out$pasets)))
  ## rcar only ever reverses, so the arc count is invariant
  stopifnot(identical(as.integer(numEdges(out$dag)),
                      as.integer(numEdges(b$dag))))
  dag_check(b$st)

  rr
}

################################################################################
## 1. random DAGs, every r, several seeds, with and without targets
## NOTE: these pass the target FAMILY, not the union of it. An arc is
## I-covered unless some single target separates its endpoints, so
## list(integer(0), c(1L, 2L)) leaves 1 -> 2 I-covered even though both ends
## are targeted -- a case the union of the targets cannot express, and the
## families below therefore include multi-vertex targets.

################################################################################

set.seed(20260907)
nrev <- 0L
ncmp <- 0L
## The sweep is 4 sizes x 3 arc sets x 4 values of r x 5 target families x 3
## seeds. It used 6 arc sets and 5 seeds, which repeated each CONFIGURATION
## ten times over and cost 5.3 s of this file's 6.3. Every dimension that
## selects a code path -- r = 0 and r = 20 at the boundaries, and the five
## target families including the two whose single target holds both ends of
## an arc -- is still crossed in full; only the repetition per configuration
## is lower.
for (p in c(3, 5, 8, 12)) {
  for (rep in 1:3) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.5, , drop=FALSE]
    if (nrow(cand) == 0) next
    cand <- cand[sample.int(nrow(cand)), , drop=FALSE]     ## shuffled inserts
    arcs <- lapply(seq_len(nrow(cand)), function(k) c(cand[k, 1], cand[k, 2]))
    for (r in c(0L, 1L, 3L, 20L))
      for (ut in list(list(integer(0)),                       # none
                      list(integer(0), 1L),                   # one singleton
                      list(integer(0), 1L, min(3L, p)),       # two singletons
                      list(integer(0), c(1L, min(2L, p))),    # BOTH ends
                      list(integer(0), c(1L, min(3L, p)), 2L)))
        for (seed in 1:3) {
          rr <- compare_rcar(both(p, arcs), r, ut, seed)
          nrev <- nrev + rr
          ncmp <- ncmp + 1L
        }
  }
}
stopifnot(ncmp > 500L)
## the sweep must actually have performed reversals, or it proves nothing
stopifnot(nrev > 0L)
cat(sprintf("rcar: %d comparisons, %d reversals performed, stream and state identical\n",
            ncmp, nrev))

################################################################################
## 2. the two early-return branches consume NO randomness
################################################################################

## (a) an edgeless DAG
b <- both(4, list())
set.seed(1); before <- .Random.seed
stopifnot(identical(dag_rcar(b$st, length(0:20), list()), 0L))
stopifnot(identical(.Random.seed, before))

## (b) arcs present, but no I-covered arc survives the target filter.
## X1 -> X2 -> X3: X1 -> X2 is covered, X2 -> X3 is not; a target naming X1
## alone separates the endpoints of X1 -> X2 and so removes the only covered
## arc. A target naming BOTH, list(integer(0), c(1L, 2L)), does not separate
## them and leaves it covered -- see tests/test_c_nh.R section 3.
b <- both(3, list(c(1,2), c(2,3)))
stopifnot(identical(as.logical(dag_ced(b$st, list())), c(TRUE, FALSE)))
stopifnot(!any(as.logical(dag_ced(b$st, list(integer(0), 1L)))),   ## precondition
          identical(as.logical(dag_ced(b$st, list(integer(0), c(1L, 2L)))),
                    c(TRUE, FALSE)))
set.seed(1); before <- .Random.seed
stopifnot(identical(dag_rcar(b$st, length(0:20), list(integer(0), 1L)), 0L))
stopifnot(identical(.Random.seed, before))
stopifnot(identical(unname(dag_edgeM(b$st)), unname(edgeMatrix(b$dag))))

## and R agrees on both branches, consuming nothing either
b <- both(4, list())
set.seed(1); before <- .Random.seed
invisible(idlBNs:::rcar(b$dag, 20L, integer(0), b$anc, b$pas, b$vidx))
stopifnot(identical(.Random.seed, before))
cat("early returns: both branches consume no randomness, in C and in R\n")

################################################################################
## 3. r = 0 still consumes exactly one draw
################################################################################

b <- both(3, list(c(1,2), c(2,3)))
set.seed(1); before <- .Random.seed
stopifnot(identical(dag_rcar(b$st, length(0:0), list()), 0L))
after <- .Random.seed
stopifnot(!identical(after, before))          ## it drew
## and it drew exactly what R drew
b2 <- both(3, list(c(1,2), c(2,3)))
set.seed(1); invisible(idlBNs:::rcar(b2$dag, 0L, integer(0), b2$anc, b2$pas, b2$vidx))
stopifnot(identical(.Random.seed, after))
cat("r = 0 consumes exactly one draw, matching R\n")

################################################################################
## 4. a chain of covered reversals: reversing a covered arc leaves it
## covered, so the mask never empties and rcar can always perform all rr
## reversals. This is the invariant that makes the "no covered arc remains"
## guard in C unreachable.
################################################################################

## X1 -> X2 -> X3 -> X4 ... a chain has exactly one covered arc at each step
b <- both(5, list(c(1,2), c(2,3), c(3,4), c(4,5)))
for (seed in 1:20) {
  bb <- both(5, list(c(1,2), c(2,3), c(3,4), c(4,5)))
  rr <- compare_rcar(bb, 20L, list(), seed)
  ## after any number of covered reversals there is still a covered arc
  stopifnot(any(as.logical(dag_ced(bb$st, list()))))
}
cat("covered reversals keep the covered set non-empty\n")

################################################################################
## 5. determinism, and independence from how the DAG was built
################################################################################

arcs <- list(c(1,3), c(1,2), c(3,4), c(2,5))
r1 <- compare_rcar(both(5, arcs), 20L, list(), 99L)
r2 <- compare_rcar(both(5, arcs), 20L, list(), 99L)
stopifnot(identical(r1, r2))

## invalid arguments are rejected before anything happens
b <- both(4, list(c(1,2)))
errs <- function(e) inherits(tryCatch(e, error=function(x) x), "error")
stopifnot(errs(dag_rcar(b$st, 0L, list())))       ## rlen must be >= 1
stopifnot(errs(dag_rcar(b$st, -1L, list())))
stopifnot(errs(.Call(idlBNs:::C_dag_rcar, b$st, NA_integer_, list())))
stopifnot(errs(.Call(idlBNs:::C_dag_rcar, b$st, 5L, "x")))
stopifnot(identical(unname(dag_edgeM(b$st)), unname(edgeMatrix(b$dag))))
dag_check(b$st)

## out-of-range target vertices ignored, as in R (see tests/test_c_nh.R section 4)
for (ut in list(list(0L), list(c(0L, 2L)), list(c(2L, 99L)), list(NA_integer_)))
  invisible(compare_rcar(both(4, list(c(1,2), c(2,3))), 5L, ut, 3L))
cat("determinism, argument validation and target tolerance all hold\n")

################################################################################
## 6. a realistic DAG from the package's own simulator, at the r the search
## actually uses
################################################################################

set.seed(11)
for (p in c(10, 20, 30)) {
  D <- r.gauss.pardag(p, 0.3, top.sort=TRUE, normalize=TRUE)
  em <- edgeMatrix(as(D, "graphNEL"))
  arcs <- lapply(seq_len(ncol(em)), function(k) c(em["from", k], em["to", k]))
  for (seed in 1:5)
    for (ut in list(list(integer(0)), list(integer(0), 2L, 5L),
                    list(integer(0), c(2L, 5L))))
      invisible(compare_rcar(both(p, arcs), 20L, ut, seed))
}
cat("simulated DAGs at p = 10, 20, 30 agree on stream and state\n")

cat("all C rcar tests passed\n")
