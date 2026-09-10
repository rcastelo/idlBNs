## 2026-09-07 differential test for the C neighbourhood enumeration against
## nr.nh(), ar.nh() and ncr.nh() in R/search.R.
##
## The assertion is identical() on the whole list(op=, u=, v=), so it pins
## the move ORDER as well as the move set. That is the point: hcmc() and
## hillclimbing() pick with which.max(), which breaks ties by position, so
## an enumeration that produces the right moves in the wrong order sends the
## search into a different local optimum while looking correct.
##
## The specific trap this is built to catch: a vertex's removals and
## reversals are emitted in the order the GRAPH stores its children, which is
## edge-insertion order and is routinely not ascending. An implementation
## that emitted them ascending would agree on every move set. So the sweep
## drives long random trajectories -- which produce plenty of removals and
## reversals -- and ASSERTS that non-ascending child lists actually occurred.

suppressPackageStartupMessages({
  library(graph)
  library(idlBNs)
})

dag_new   <- function(p) .Call(idlBNs:::C_dag_new, as.integer(p))
dag_move  <- function(st, op, u, v)
    invisible(.Call(idlBNs:::C_dag_apply_move, st, as.integer(op),
                    as.integer(u), as.integer(v)))
dag_nh    <- function(st, kind, ut = list())
    .Call(idlBNs:::C_dag_nh, st, as.integer(kind), ut)
dag_ced   <- function(st, ut) .Call(idlBNs:::C_dag_cedges, st, ut)
dag_check <- function(st) invisible(.Call(idlBNs:::C_dag_check, st))

## the R neighbourhood, in the same shape C returns
r_nh <- function(kind, dag, anc, ut) {
  ne <- switch(kind,
               idlBNs:::nr.nh(dag, anc),
               idlBNs:::ar.nh(dag, anc),
               idlBNs:::ncr.nh(dag, anc, ut))

  list(op=ne$op, u=ne$u, v=ne$v)
}

## target FAMILIES, not unions: the last two hold both endpoints of a
## possible arc in one target, which leaves that arc I-covered even though
## both ends are targeted -- the case a union of targets cannot express.
UTARGETS <- list(list(integer(0)),
                 list(integer(0), 1L),
                 list(integer(0), 1L, 3L),
                 list(integer(0), c(1L, 3L)),
                 list(integer(0), c(1L, 2L), 3L))

ncmp <- 0L
nmoves <- 0L
nonasc <- 0L
nrev <- 0L
ncov <- 0L

compare_nh <- function(st, dag, anc, p) {
  for (ut in UTARGETS) {
    ## ncr is the only kind that reads the targets, so only vary it there
    kinds <- if (identical(ut, list(integer(0)))) 1:3 else 3L
    for (kd in kinds) {
      got <- dag_nh(st, kd, ut)
      want <- r_nh(kd, dag, anc, ut)
      stopifnot(identical(got, want))
      ncmp <<- ncmp + 1L
      if (kd == 2L)
        nrev <<- nrev + sum(got$op == 3L)
    }
    ## and the covered-arc mask, in edgeMatrix order
    mC <- as.logical(dag_ced(st, ut))
    mR <- as.logical(idlBNs:::cedges(dag, ut))
    stopifnot(identical(mC, mR))
    ncov <<- ncov + sum(mR)
  }
  invisible(NULL)
}

################################################################################
## 1. long random trajectories, comparing every neighbourhood at every step
################################################################################

shadow_nh <- function(p, nsteps, seed) {
  set.seed(seed)
  vn <- paste0("X", seq_len(p))
  dag <- new("graphNEL", nodes=vn, edgemode="directed")
  anc <- idlBNs:::init.ancestors(vn)
  st <- dag_new(p)

  compare_nh(st, dag, anc, p)             ## the edgeless DAG

  for (k in seq_len(nsteps)) {
    ne <- idlBNs:::ar.nh(dag, anc)
    stopifnot(length(ne$op) > 0L)
    m <- sample.int(length(ne$op), 1L)
    op <- ne$op[m]; u <- ne$u[m]; v <- ne$v[m]

    anc <- switch(op,
                  idlBNs:::add.ancestors(anc, vn[u], vn[v]),
                  idlBNs:::remove.ancestors(anc, dag, vn[u], vn[v]),
                  idlBNs:::reverse.ancestors(anc, dag, vn[u], vn[v]))
    dag <- idlBNs:::apply.move(dag, op, u, v, vn)
    dag_move(st, op, u, v)

    compare_nh(st, dag, anc, p)
    nmoves <<- nmoves + 1L

    el <- edgeL(dag)
    for (i in seq_len(p)) {
      e <- el[[i]]$edges
      if (!is.null(e) && length(e) > 1L &&
          !identical(as.integer(e), sort(as.integer(e))))
        nonasc <<- nonasc + 1L
    }
  }
  dag_check(st)
  invisible(NULL)
}

## 200 steps at p = 8 and p = 12 is where child lists churn hardest
for (p in c(4, 6, 8, 12))
  for (nsteps in c(0, 1, 7, 200))
    for (seed in 1:4)
      shadow_nh(p, nsteps, seed)
## a couple of wider graphs, fewer steps, to exercise larger neighbourhoods
for (p in c(20, 30))
  for (seed in 1:2)
    shadow_nh(p, 40, seed)

stopifnot(nmoves > 2000L, ncmp > 5000L)
## the three anti-vacuity guards: the sweep must actually have produced
## non-ascending child lists, offered reversals, and seen covered arcs
stopifnot(nonasc > 0L, nrev > 0L, ncov > 0L)
cat(sprintf("neighbourhoods: %d comparisons over %d moves, all identical\n",
            ncmp, nmoves))
cat(sprintf("  non-ascending child lists seen : %d\n", nonasc))
cat(sprintf("  reversal moves offered         : %d\n", nrev))
cat(sprintf("  covered arcs seen              : %d\n", ncov))

################################################################################
## 2. the trap, directly: arcs inserted out of ascending order
################################################################################

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
  list(vn=vn, dag=dag, st=st, anc=anc)
}

## X1's children are inserted 4, 2, 3 -- so a correct enumeration emits its
## removals in that order, and an ascending one would emit 2, 3, 4
b <- both(6, list(c(1,4), c(1,2), c(1,3), c(2,5)))
stopifnot(identical(as.integer(edgeL(b$dag)[[1]]$edges), c(4L, 2L, 3L)))
got <- dag_nh(b$st, 1L)
rem <- got$op == 2L & got$u == 1L
stopifnot(identical(got$v[rem], c(4L, 2L, 3L)))
stopifnot(!identical(got$v[rem], sort(got$v[rem])))   ## really not ascending
for (kd in 1:3)
  stopifnot(identical(dag_nh(b$st, kd), r_nh(kd, b$dag, b$anc, integer(0))))
cat("insertion order: removals follow stored child order, not ascending\n")

################################################################################
## 3. the I-covered arc rule of ncr, including the target exemption
################################################################################

## a chain X1 -> X2 -> X3. X1 -> X2 is covered: pa(X1) is empty and
## pa(X2) \ {X1} is empty. X2 -> X3 is not: pa(X2) = {X1} but
## pa(X3) \ {X2} is empty.
b <- both(3, list(c(1,2), c(2,3)))
stopifnot(identical(as.logical(dag_ced(b$st, list())), c(TRUE, FALSE)))
## ncr therefore drops the reversal of X1 -> X2 but keeps X2 -> X3's
ncr <- dag_nh(b$st, 3L, list())
rv <- ncr$op == 3L
stopifnot(identical(ncr$u[rv], 2L), identical(ncr$v[rv], 3L))
## ar keeps both reversals
ar <- dag_nh(b$st, 2L)
stopifnot(sum(ar$op == 3L) == 2L)
## A target that SEPARATES the endpoints puts the reversal back: the arc is
## then target-protected, hence not I-covered.
for (ut in list(list(integer(0), 1L), list(integer(0), 2L),
                list(integer(0), c(1L, 3L)))) {
  ncr2 <- dag_nh(b$st, 3L, ut)
  stopifnot(sum(ncr2$op == 3L) == 2L)
  stopifnot(identical(dag_nh(b$st, 3L, ut), r_nh(3L, b$dag, b$anc, ut)))
}
## But a single target holding BOTH endpoints does not separate them, so the
## arc stays I-covered and the reversal stays dropped -- even though both
## endpoints are targeted. This is the case the union of the targets cannot
## express, and testing it against the union of list(integer(0), c(1L, 2L))
## would wrongly expect 2.
for (ut in list(list(integer(0), c(1L, 2L)),
                list(integer(0), c(1L, 2L, 3L)))) {
  ncr2 <- dag_nh(b$st, 3L, ut)
  stopifnot(sum(ncr2$op == 3L) == 1L,
            identical(as.logical(dag_ced(b$st, ut)), c(TRUE, FALSE)),
            identical(ncr2, r_nh(3L, b$dag, b$anc, ut)))
}
## a target elsewhere leaves it dropped
stopifnot(sum(dag_nh(b$st, 3L, list(integer(0), 3L))$op == 3L) == 1L)
cat("ncr: I-covered arcs dropped, target-touching covered arcs kept\n")

################################################################################
## 4. out-of-range utargets are IGNORED, not rejected.
##
## Not because such a value is meaningful -- it is not. A targets family says
## which vertices each intervention acts on, and "no intervention" is
## list(integer(0), ...), not list(0L, ...). But utargets is derived as
## sort(unique(unlist(targets))) and then used only through `%in%` against
## vertex indices, so R never errors on a value outside 1..p; it silently
## matches nothing. The C port has to reproduce that, because the whole
## verification strategy is that C and R agree -- including on malformed
## input. Validating targets belongs at the public boundary (hcmc() /
## hillclimbing()), not here.
################################################################################

b <- both(4, list(c(1,2), c(2,3)))
for (ut in list(list(0L), list(c(0L, 2L)), list(c(-1L, 2L)),
                list(c(2L, 99L)), list(NA_integer_))) {
  stopifnot(identical(dag_nh(b$st, 3L, ut), r_nh(3L, b$dag, b$anc, ut)))
  stopifnot(identical(as.logical(dag_ced(b$st, ut)),
                      as.logical(idlBNs:::cedges(b$dag, ut))))
}
cat("out-of-range utargets ignored, matching R\n")

################################################################################
## 5. degenerate and boundary shapes
################################################################################

## p = 1: no move of any kind is possible
s1 <- dag_new(1L)
for (kd in 1:3)
  stopifnot(length(dag_nh(s1, kd)$op) == 0L)
stopifnot(length(dag_ced(s1, list())) == 0L)

## p = 2 edgeless: two additions, no removals, no reversals
s2 <- dag_new(2L)
n2 <- dag_nh(s2, 2L)
stopifnot(identical(n2$op, c(1L, 1L)), identical(n2$u, c(1L, 2L)),
          identical(n2$v, c(2L, 1L)))

## a complete DAG on p vertices offers no additions at all
for (p in c(4, 6)) {
  vn <- paste0("X", seq_len(p))
  arcs <- list()
  for (i in seq_len(p - 1)) for (j in (i + 1):p) arcs <- c(arcs, list(c(i, j)))
  b <- both(p, arcs)
  got <- dag_nh(b$st, 2L)
  stopifnot(!any(got$op == 1L))                       ## nothing left to add
  stopifnot(sum(got$op == 2L) == p * (p - 1) / 2)     ## every arc removable
  for (kd in 1:3)
    stopifnot(identical(dag_nh(b$st, kd), r_nh(kd, b$dag, b$anc, integer(0))))
}

## the multi-word boundary: a chain across p = 65 and p = 129 vertices
for (p in c(63, 64, 65, 129)) {
  vn <- paste0("X", seq_len(p))
  b <- both(p, lapply(seq_len(p - 1), function(i) c(i, i + 1L)))
  for (kd in 1:3)
    stopifnot(identical(dag_nh(b$st, kd), r_nh(kd, b$dag, b$anc, integer(0))))
  dag_check(b$st)
}
cat("boundary shapes: p = 1, 2, complete DAGs, and p = 63/64/65/129 agree\n")

cat("all C neighbourhood tests passed\n")
