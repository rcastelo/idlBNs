## 2026-09-07 direct regression tests for the search-loop internals that had
## no direct coverage at all: remove.ancestors(), reverse.ancestors(),
## cedges() and rcar(). They were reachable only through a full hcmc() /
## hillclimbing() run, whose only assertions are the pinned .Rout.save
## outputs and the idlBNs.debug.pasets option -- and that option checks
## parent sets, never the ancestor matrix. So a bug in the ancestor
## bookkeeping that merely changed WHICH moves get enumerated was caught by
## nothing but a trajectory pin.
##
## These tests exist to be the oracle for a C port of the same functions, so
## they assert exact structural equality, never tolerances.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

## ---------------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------------

## build a graphNEL from an edge list, inserting in the given order (so child
## lists end up in that order, which is load-bearing -- see test_pasets.R)
mkdag <- function(vn, arcs) {
  g <- new("graphNEL", nodes=vn, edgemode="directed")
  for (a in arcs)
    g <- addEdge(vn[a[1]], vn[a[2]], g)
  g
}

## the ancestor matrix computed from scratch: transitive closure of the
## adjacency matrix by repeated boolean squaring. shares no machinery with
## add.ancestors()/remove.ancestors(), which is the point.
anc_closure <- function(dag) {
  vn <- nodes(dag)
  p <- length(vn)
  A <- matrix(FALSE, p, p, dimnames=list(vn, vn))
  em <- edgeMatrix(dag)
  if (ncol(em) > 0)
    A[cbind(em["from", ], em["to", ])] <- TRUE
  R <- A
  repeat {
    N <- R | ((R %*% A) > 0)
    if (identical(N, R)) break
    R <- N
  }
  ## R[u, v] is TRUE iff a directed path u -> ... -> v exists, which is
  ## exactly the documented meaning of anc[u, v]
  R
}

## build the incrementally maintained ancestor matrix for a DAG, the way the
## search does: from an edgeless graph, one add.ancestors() per arc
anc_incremental <- function(dag) {
  vn <- nodes(dag)
  anc <- idlBNs:::init.ancestors(vn)
  em <- edgeMatrix(dag)
  for (k in seq_len(ncol(em)))
    anc <- idlBNs:::add.ancestors(anc, vn[em["from", k]], vn[em["to", k]])
  anc
}

## ---------------------------------------------------------------------------
## 1. add.ancestors() agrees with a from-scratch transitive closure
## ---------------------------------------------------------------------------

set.seed(20260907)
for (p in c(2, 3, 5, 8, 12)) {
  vn <- paste0("X", seq_len(p))
  for (rep in 1:20) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.35, , drop=FALSE]
    arcs <- list()
    if (nrow(cand) > 0) {
      cand <- cand[sample.int(nrow(cand)), , drop=FALSE]
      arcs <- lapply(seq_len(nrow(cand)), function(k) c(cand[k, 1], cand[k, 2]))
    }
    dag <- mkdag(vn, arcs)
    stopifnot(identical(anc_incremental(dag), anc_closure(dag)))
  }
}

## ---------------------------------------------------------------------------
## 2. remove.ancestors() -- hand-built topologies that exercise the
## descendant recomputation, then a randomised sweep
## ---------------------------------------------------------------------------

## each case: vertex count, arcs (in insertion order), and the arc to remove
cases <- list(
  ## a chain: removing the middle arc splits ancestry in two
  list(p=4, arcs=list(c(1,2), c(2,3), c(3,4)), rm=c(2,3)),
  ## a diamond: removing one side leaves ancestry intact via the other
  list(p=4, arcs=list(c(1,2), c(1,3), c(2,4), c(3,4)), rm=c(2,4)),
  ## a diamond, removing the arc into the apex
  list(p=4, arcs=list(c(1,2), c(1,3), c(2,4), c(3,4)), rm=c(1,2)),
  ## head with several parents: removing one leaves the others
  list(p=5, arcs=list(c(1,5), c(2,5), c(3,5), c(4,5)), rm=c(2,5)),
  ## removal that disconnects a whole subtree from its ancestors
  list(p=6, arcs=list(c(1,2), c(2,3), c(3,4), c(3,5), c(5,6)), rm=c(2,3)),
  ## a single arc, so the descendant set is just the head
  list(p=2, arcs=list(c(1,2)), rm=c(1,2)),
  ## two independent components, removing from one must not touch the other
  list(p=6, arcs=list(c(1,2), c(2,3), c(4,5), c(5,6)), rm=c(1,2)),
  ## reconvergent: 1 reaches 4 by two different lengths
  list(p=4, arcs=list(c(1,2), c(1,4), c(2,3), c(3,4)), rm=c(1,4))
)

for (cs in cases) {
  vn <- paste0("X", seq_len(cs$p))
  dag <- mkdag(vn, cs$arcs)
  u <- vn[cs$rm[1]]
  v <- vn[cs$rm[2]]
  anc <- anc_incremental(dag)
  ## remove.ancestors() takes the DAG as it stood BEFORE the removal
  got <- idlBNs:::remove.ancestors(anc, dag, u, v)
  want <- anc_closure(removeEdge(u, v, dag))
  stopifnot(identical(got, want))
}

## randomised sweep: remove every arc of many random DAGs, one at a time
for (p in c(3, 5, 8, 12)) {
  vn <- paste0("X", seq_len(p))
  for (rep in 1:15) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.4, , drop=FALSE]
    if (nrow(cand) == 0) next
    cand <- cand[sample.int(nrow(cand)), , drop=FALSE]
    dag <- mkdag(vn, lapply(seq_len(nrow(cand)),
                            function(k) c(cand[k, 1], cand[k, 2])))
    anc <- anc_incremental(dag)
    em <- edgeMatrix(dag)
    for (k in seq_len(ncol(em))) {
      u <- vn[em["from", k]]
      v <- vn[em["to", k]]
      stopifnot(identical(idlBNs:::remove.ancestors(anc, dag, u, v),
                          anc_closure(removeEdge(u, v, dag))))
    }
  }
}

## ---------------------------------------------------------------------------
## 3. remove.ancestors() is independent of WHICH valid topological order it
## walks. this is the property that licenses replacing RBGL::tsort with an
## internal topological sort, so it is asserted here rather than assumed.
## ---------------------------------------------------------------------------

## a variant of remove.ancestors() parameterised by the topological order, so
## several valid orders can be compared. mirrors R/search.R:34-56 exactly
## apart from where 'to' comes from.
remove_anc_with_order <- function(anc, dag, u, v, to) {
  vnodes <- nodes(dag)
  em <- edgeMatrix(dag)
  pasets <- split(vnodes[em["from", ]],
                  factor(vnodes[em["to", ]], levels=vnodes))
  D <- c(v, rownames(anc)[anc[v, ]])
  D <- to[to %in% D]
  for (k in D) {
    pa.k <- pasets[[k]]
    if (k == v)
      pa.k <- setdiff(pa.k, u)
    newcol <- rep(FALSE, length(vnodes))
    names(newcol) <- vnodes
    if (length(pa.k) > 0) {
      newcol[pa.k] <- TRUE
      newcol <- newcol | apply(anc[, pa.k, drop=FALSE], 1, any)
    }
    anc[, k] <- newcol
  }
  anc
}

## is 'to' a valid topological order of dag?
is_topo <- function(dag, to) {
  vn <- nodes(dag)
  pos <- setNames(seq_along(to), to)
  em <- edgeMatrix(dag)
  if (ncol(em) == 0) return(TRUE)
  all(pos[vn[em["from", ]]] < pos[vn[em["to", ]]])
}

set.seed(99)
norders <- 0L
for (cs in cases) {
  vn <- paste0("X", seq_len(cs$p))
  dag <- mkdag(vn, cs$arcs)
  u <- vn[cs$rm[1]]
  v <- vn[cs$rm[2]]
  anc <- anc_incremental(dag)
  ref <- idlBNs:::remove.ancestors(anc, dag, u, v)
  ## try many random permutations, keeping the ones that are valid
  ## topological orders, and assert each gives the identical result
  for (t in 1:200) {
    to <- sample(vn)
    if (!is_topo(dag, to)) next
    norders <- norders + 1L
    stopifnot(identical(remove_anc_with_order(anc, dag, u, v, to), ref))
  }
}
stopifnot(norders > 20L)   ## the sweep actually exercised alternative orders
cat(sprintf("remove.ancestors: order-independent over %d valid topological orders\n",
            norders))

## ---------------------------------------------------------------------------
## 4. reverse.ancestors()
## ---------------------------------------------------------------------------

for (cs in cases) {
  vn <- paste0("X", seq_len(cs$p))
  dag <- mkdag(vn, cs$arcs)
  u <- vn[cs$rm[1]]
  v <- vn[cs$rm[2]]
  ## only meaningful when the reversal keeps the graph acyclic
  rev.g <- addEdge(v, u, removeEdge(u, v, dag))
  if (any(diag(anc_closure(rev.g))))
    next
  anc <- anc_incremental(dag)
  stopifnot(identical(idlBNs:::reverse.ancestors(anc, dag, u, v),
                      anc_closure(rev.g)))
}

## and over random DAGs, for every arc whose reversal stays acyclic
for (p in c(3, 5, 8)) {
  vn <- paste0("X", seq_len(p))
  for (rep in 1:15) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.4, , drop=FALSE]
    if (nrow(cand) == 0) next
    cand <- cand[sample.int(nrow(cand)), , drop=FALSE]
    dag <- mkdag(vn, lapply(seq_len(nrow(cand)),
                            function(k) c(cand[k, 1], cand[k, 2])))
    anc <- anc_incremental(dag)
    em <- edgeMatrix(dag)
    for (k in seq_len(ncol(em))) {
      u <- vn[em["from", k]]
      v <- vn[em["to", k]]
      rev.g <- addEdge(v, u, removeEdge(u, v, dag))
      if (any(diag(anc_closure(rev.g))))
        next
      stopifnot(identical(idlBNs:::reverse.ancestors(anc, dag, u, v),
                          anc_closure(rev.g)))
    }
  }
}

## ---------------------------------------------------------------------------
## 5. cedges() -- the covered / I-covered arc mask
## ---------------------------------------------------------------------------

## independent reference: arc u -> v is covered iff pa(u) == pa(v) \ {u}, and
## is excluded when either endpoint is an intervention target. Returned in
## edgeMatrix() column order, because rcar() selects by index into it.
cedges_reference <- function(dag, utargets) {
  vn <- nodes(dag)
  em <- edgeMatrix(dag)
  if (ncol(em) == 0)
    return(logical(0))
  pa <- lapply(seq_along(vn),
               function(j) sort(as.integer(em["from", em["to", ] == j])))
  vapply(seq_len(ncol(em)), function(k) {
           u <- em["from", k]
           v <- em["to", k]
           covered <- identical(pa[[u]], sort(setdiff(pa[[v]], u)))
           touches <- (u %in% utargets) || (v %in% utargets)
           covered && !touches
         }, logical(1))
}

set.seed(7)
ncov <- 0L
for (p in c(2, 3, 5, 8, 12)) {
  vn <- paste0("X", seq_len(p))
  for (rep in 1:20) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.4, , drop=FALSE]
    arcs <- list()
    if (nrow(cand) > 0) {
      cand <- cand[sample.int(nrow(cand)), , drop=FALSE]
      arcs <- lapply(seq_len(nrow(cand)), function(k) c(cand[k, 1], cand[k, 2]))
    }
    dag <- mkdag(vn, arcs)
    for (ut in list(integer(0), c(1L), c(1L, min(3L, p)), seq(1L, p, by=2L))) {
      got <- idlBNs:::cedges(dag, ut)
      want <- cedges_reference(dag, ut)
      stopifnot(identical(as.logical(got), want))
      ncov <- ncov + sum(want)
    }
  }
}
stopifnot(ncov > 0L)   ## covered arcs were actually encountered
cat(sprintf("cedges: agrees with an independent reference (%d covered arcs seen)\n",
            ncov))

## a chain X1->X2->X3 has exactly one covered arc, X1->X2 (pa(X1)=empty,
## pa(X2)\{X1}=empty); X2->X3 is not covered (pa(X2)={X1} but
## pa(X3)\{X2}=empty)
g3 <- mkdag(paste0("X", 1:3), list(c(1,2), c(2,3)))
em3 <- edgeMatrix(g3)
m3 <- idlBNs:::cedges(g3, integer(0))
stopifnot(identical(as.logical(m3), c(TRUE, FALSE)))
## and naming either endpoint of X1->X2 as a target un-covers it
stopifnot(identical(as.logical(idlBNs:::cedges(g3, 1L)), c(FALSE, FALSE)))
stopifnot(identical(as.logical(idlBNs:::cedges(g3, 2L)), c(FALSE, FALSE)))
stopifnot(identical(as.logical(idlBNs:::cedges(g3, 3L)), c(TRUE, FALSE)))

## ---------------------------------------------------------------------------
## 6. rcar() -- behaviour, invariants, and RNG stream consumption
## ---------------------------------------------------------------------------

## rcar() reverses only covered arcs, which is score-preserving under a
## decomposable score, so it must preserve the vertex set, the edge count,
## and acyclicity, and it must keep anc and pasets in sync with the DAG.
set.seed(4242)
nrev <- 0L
for (p in c(3, 5, 8)) {
  vn <- paste0("X", seq_len(p))
  vidx <- setNames(seq_len(p), vn)
  datm <- matrix(rnorm(20 * p), nrow=20, dimnames=list(NULL, vn))
  for (rep in 1:20) {
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind=TRUE)
    cand <- cand[runif(nrow(cand)) < 0.5, , drop=FALSE]
    if (nrow(cand) == 0) next
    cand <- cand[sample.int(nrow(cand)), , drop=FALSE]
    dag <- mkdag(vn, lapply(seq_len(nrow(cand)),
                            function(k) c(cand[k, 1], cand[k, 2])))
    anc <- anc_incremental(dag)
    pas <- idlBNs:::.build_pasets(dag, datm)
    for (r in c(0L, 1L, 5L, 20L)) {
      out <- idlBNs:::rcar(dag, r, integer(0), anc, pas, vidx)
      ## structure preserved
      stopifnot(identical(nodes(out$dag), vn),
                numEdges(out$dag) == numEdges(dag),
                !any(diag(anc_closure(out$dag))))
      ## anc and pasets stayed in sync with the returned DAG
      stopifnot(identical(out$anc, anc_closure(out$dag)))
      stopifnot(identical(lapply(out$pasets, function(x) unname(sort.int(x))),
                          lapply(idlBNs:::.build_pasets(out$dag, datm),
                                 function(x) unname(sort.int(x)))))
      if (!identical(unname(edgeMatrix(out$dag)), unname(edgeMatrix(dag))))
        nrev <- nrev + 1L
    }
  }
}
stopifnot(nrev > 0L)   ## reversals actually happened
cat(sprintf("rcar: %d calls actually reversed at least one arc\n", nrev))

## the two early-return branches must consume NO randomness at all, which is
## what lets a C port put them above GetRNGstate()
vn4 <- paste0("X", 1:4)
vidx4 <- setNames(1:4, vn4)
dat4 <- matrix(rnorm(20 * 4), nrow=20, dimnames=list(NULL, vn4))

## (a) an edgeless DAG: numEdges(dag) == 0
g0 <- new("graphNEL", nodes=vn4, edgemode="directed")
a0 <- idlBNs:::init.ancestors(vn4)
p0 <- idlBNs:::init.pasets(4L)
set.seed(1); before <- .Random.seed
out0 <- idlBNs:::rcar(g0, 20L, integer(0), a0, p0, vidx4)
stopifnot(identical(.Random.seed, before))
stopifnot(numEdges(out0$dag) == 0L, identical(out0$anc, a0),
          identical(out0$pasets, p0))

## (b) a DAG with edges but no covered arc left after target masking
g1 <- mkdag(vn4, list(c(1,2), c(2,3)))
a1 <- anc_incremental(g1)
p1 <- idlBNs:::.build_pasets(g1, dat4)
stopifnot(!any(idlBNs:::cedges(g1, c(1L, 2L, 3L))))   ## precondition
set.seed(1); before <- .Random.seed
out1 <- idlBNs:::rcar(g1, 20L, c(1L, 2L, 3L), a1, p1, vidx4)
stopifnot(identical(.Random.seed, before))
stopifnot(identical(unname(edgeMatrix(out1$dag)), unname(edgeMatrix(g1))))

## r = 0 still draws once (sample(0:0, 1) goes through sample.int(1, 1)), so a
## C port must not short-circuit it
set.seed(1); before <- .Random.seed
invisible(idlBNs:::rcar(g1, 0L, integer(0), a1, p1, vidx4))
stopifnot(!identical(.Random.seed, before))

## rcar() is deterministic given the seed
set.seed(5); r1 <- idlBNs:::rcar(g1, 20L, integer(0), a1, p1, vidx4)
set.seed(5); r2 <- idlBNs:::rcar(g1, 20L, integer(0), a1, p1, vidx4)
stopifnot(identical(unname(edgeMatrix(r1$dag)), unname(edgeMatrix(r2$dag))),
          identical(r1$anc, r2$anc), identical(r1$pasets, r2$pasets))

## ---------------------------------------------------------------------------
## 7. resample() consumes exactly what sample.int() does, including the
## length-1 case a C port must not short-circuit
## ---------------------------------------------------------------------------

for (len in 1:12)
  for (s in 1:20) {
    x <- seq_len(len) * 7L
    set.seed(s); a <- idlBNs:::resample(x, size=1); sa <- .Random.seed
    set.seed(s); b <- x[sample.int(len, size=1)];   sb <- .Random.seed
    stopifnot(identical(a, b), identical(sa, sb))
  }

cat("all search internals tests passed\n")
