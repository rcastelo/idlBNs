## 2026-09-07 tests for the error-bounded candidate band that finds the best
## move without summing all p vertex terms for every candidate.
##
## WHY A BAND AT ALL. total_with() sums all p terms per candidate, which is
## what makes each total bit-identical to a full re-score, but the
## neighbourhood holds O(p^2) candidates, so a step costs O(p^3) additions:
## about 125M at p = 500, roughly 125 ms per step. The band computes a cheap
## estimate per candidate, bounds its error rigorously, and pays the exact
## O(p) summation only for the candidates whose interval could still contain
## the maximum.
##
## WHAT MUST HOLD is not "approximately the right move" but exactly the move
## and exactly the total a full exact ranking would give -- including which
## of several tied candidates wins, since which.max() takes the first and the
## band scans in ascending order with a strict '>'. That is asserted here
## directly against nh.scores.fun's full vector.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

## build a random DAG plus the structures the scorers need
setup <- function(p, n = 120L, seed = 1L, dens = 0.4, shuffle = TRUE) {
  set.seed(seed)
  vn <- paste0("X", seq_len(p))
  dat <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, vn))
  attr(dat, "sanitycheck") <- TRUE
  cand <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  cand <- cand[runif(nrow(cand)) < dens, , drop = FALSE]
  g <- new("graphNEL", nodes = vn, edgemode = "directed")
  pas <- idlBNs:::init.pasets(p)
  if (nrow(cand) > 0) {
    if (shuffle)
      cand <- cand[sample.int(nrow(cand)), , drop = FALSE]
    for (k in seq_len(nrow(cand))) {
      g <- addEdge(vn[cand[k, 1]], vn[cand[k, 2]], g)
      pas <- idlBNs:::add.pasets(pas, as.integer(cand[k, 1]),
                                 as.integer(cand[k, 2]))
    }
  }
  anc <- idlBNs:::init.ancestors(vn)
  em <- edgeMatrix(g)
  for (k in seq_len(ncol(em)))
    anc <- idlBNs:::add.ancestors(anc, vn[em["from", k]], vn[em["to", k]])

  list(vn = vn, dat = dat, g = g, pas = pas, anc = anc)
}

################################################################################
## 1. the band's winner and total equal a full exact ranking's, exactly
################################################################################

ncmp <- 0L
nband <- 0L
nbandmax <- 0L
worst <- 0
set.seed(20260908)
for (p in c(4, 6, 9, 14, 20)) {
  for (rep in 1:4) {
    S <- setup(p, seed = 100L * p + rep)
    for (kd in 1:3) {
      ne <- switch(kd, idlBNs:::nr.nh(S$g, S$anc), idlBNs:::ar.nh(S$g, S$anc),
                   idlBNs:::ncr.nh(S$g, S$anc, c(1L, min(3L, p))))
      if (length(ne$op) == 0L)
        next
      for (sf in list(iBIC, iBGe)) {
        gs <- attr(sf, "global.sufstats.fun")(S$dat, list(integer(0)),
                                              rep(1L, nrow(S$dat)))
        for (cs in list(NULL,
                        .Call(idlBNs:::C_sccache_new, as.integer(p)))) {
          sco <- attr(sf, "nh.scores.fun")(ne$op, ne$u, ne$v, S$pas, gs, cs)
          am <- attr(sf, "nh.argmax.fun")(ne$op, ne$u, ne$v, S$pas, gs, cs,
                                          TRUE)
          ## the same index -- not merely an equally good one
          stopifnot(identical(am$index, which.max(sco)))
          ## and the same total, bit for bit
          stopifnot(identical(am$total, max(sco)))
          ## the bound held at every candidate (verify = TRUE would have
          ## errored otherwise) and was not close to being tight
          stopifnot(am$worst < 1)
          worst <- max(worst, am$worst)
          nband <- nband + am$band
          nbandmax <- max(nbandmax, am$band)
          ncmp <- ncmp + 1L
        }
      }
    }
  }
}
stopifnot(ncmp > 100L)
cat(sprintf("band: %d neighbourhoods, argmax and total exact in every one\n",
            ncmp))
cat(sprintf("  worst |T - est| / bound  : %.3e   (must be < 1)\n", worst))
cat(sprintf("  mean band size           : %.2f candidates (max %d)\n",
            nband / ncmp, nbandmax))

################################################################################
## 2. score-equivalence ties, which is where the band's tie-breaking has to
## match which.max()'s exactly. On purely observational data both iBIC and
## iBGe give Markov equivalent DAGs the same score, so a neighbourhood of an
## edgeless graph is full of exactly tied candidates: adding u -> v scores
## the same as adding v -> u.
################################################################################

nties <- 0L
maxties <- 0L
for (p in c(3, 5, 8, 12)) {
  S <- setup(p, seed = 7L * p, dens = 0)          ## edgeless
  ne <- idlBNs:::ar.nh(S$g, S$anc)
  for (sf in list(iBIC, iBGe)) {
    gs <- attr(sf, "global.sufstats.fun")(S$dat, list(integer(0)),
                                          rep(1L, nrow(S$dat)))
    sco <- attr(sf, "nh.scores.fun")(ne$op, ne$u, ne$v, S$pas, gs, NULL)
    am <- attr(sf, "nh.argmax.fun")(ne$op, ne$u, ne$v, S$pas, gs, NULL, TRUE)
    ## there really are exact ties here, or this proves nothing
    ties <- sum(sco == max(sco))
    stopifnot(ties >= 1L)
    maxties <- max(maxties, ties)
    nties <- nties + ties
    ## and the band picks the FIRST of them, as which.max() does
    stopifnot(identical(am$index, which.max(sco)))
    stopifnot(identical(am$total, max(sco)))
    ## the proof says the true maximum lies inside the band; when several
    ## candidates attain it exactly, EVERY one of them must therefore be in
    ## it -- a candidate excluded from the band is provably strictly worse,
    ## which an exactly tied one is not. so the band is at least as large as
    ## the tie class.
    stopifnot(am$band >= ties)
  }
}
stopifnot(nties > 8L)
stopifnot(maxties > 1L)   ## genuine multi-way ties actually occurred
cat(sprintf("score-equivalence ties: %d tied maxima (largest tie class %d), first index chosen every time\n",
            nties, maxties))

################################################################################
## 3. searches: both engines must still agree, and with the band verified at
## every single step
################################################################################

mkdata <- function(p, k = 2L, n = 120L, seed = 1L) {
  set.seed(seed)
  Mg <- r.gauss.pardag(p, 0.3, top.sort = TRUE, normalize = TRUE)
  I <- c(list(integer(0)), sample(p, size = k, replace = FALSE))
  nb <- rep(floor(n / (k + 1)), k)
  nb <- c(n - sum(nb), nb)
  dat <- do.call(rbind,
                 lapply(seq_along(I),
                        function(v) rmvnorm.ivent(nb[v], Mg, target = I[[v]],
                                                  target.value = rep(2, length(I[[v]])))))

  list(dat = dat, targets = I, target.index = rep(seq_along(nb), nb))
}

old <- getOption("idlBNs.debug.band", FALSE)
on.exit(options(idlBNs.debug.band = old), add = TRUE)
options(idlBNs.debug.band = TRUE)

nsearch <- 0L
for (p in c(6, 10, 16))
  for (seed in 1:2) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe))
      for (alg in c("hcmc", "hillclimbing")) {
        f <- match.fun(alg)
        set.seed(99L)
        rC <- f(x$dat, targets = x$targets, target.index = x$target.index,
                scorefun = sf, verbose = FALSE, engine = "C")
        sC <- .Random.seed
        set.seed(99L)
        rR <- f(x$dat, targets = x$targets, target.index = x$target.index,
                scorefun = sf, verbose = FALSE, engine = "R")
        stopifnot(identical(rC$sco, rR$sco),
                  identical(unname(edgeMatrix(rC$dag)),
                            unname(edgeMatrix(rR$dag))),
                  identical(sC, .Random.seed))
        nsearch <- nsearch + 1L
      }
  }
options(idlBNs.debug.band = old)
cat(sprintf("searches: %d cases agree with the band verified at every step\n",
            nsearch))

################################################################################
## 4. the bound under adversarial per-vertex terms.
##
## The band's correctness rests on |T_m - est_m| <= SAFETY * B_m. The hardest
## case for that is massive cancellation: term magnitudes far larger than the
## total, so that sum|base_i| dwarfs |sum base_i| and the summation error is
## as large as it can be relative to the answer. A bound that kept only the
## gamma * A_m term and dropped T_base's own error would fail here.
##
## This is driven through the real scorer by constructing data whose vertex
## terms differ by many orders of magnitude, which is what widely different
## per-column variances produce.
################################################################################

set.seed(4242)
nadv <- 0L
for (trial in 1:40) {
  p <- sample(4:10, 1)
  n <- 80
  vn <- paste0("X", seq_len(p))
  ## column scales spanning ~1e-6 to ~1e6, so the per-vertex iBIC terms span
  ## a very wide range and their absolute sum greatly exceeds the total
  sc <- 10^runif(p, -6, 6)
  dat <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, vn))
  dat <- sweep(dat, 2, sc, `*`)
  attr(dat, "sanitycheck") <- TRUE
  g <- new("graphNEL", nodes = vn, edgemode = "directed")
  pas <- idlBNs:::init.pasets(p)
  cand <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  cand <- cand[runif(nrow(cand)) < 0.5, , drop = FALSE]
  if (nrow(cand) > 0) {
    cand <- cand[sample.int(nrow(cand)), , drop = FALSE]
    for (k in seq_len(nrow(cand))) {
      g <- addEdge(vn[cand[k, 1]], vn[cand[k, 2]], g)
      pas <- idlBNs:::add.pasets(pas, as.integer(cand[k, 1]),
                                 as.integer(cand[k, 2]))
    }
  }
  anc <- idlBNs:::init.ancestors(vn)
  em <- edgeMatrix(g)
  for (k in seq_len(ncol(em)))
    anc <- idlBNs:::add.ancestors(anc, vn[em["from", k]], vn[em["to", k]])
  ne <- idlBNs:::ar.nh(g, anc)
  if (length(ne$op) == 0L)
    next
  for (sf in list(iBIC, iBGe)) {
    gs <- attr(sf, "global.sufstats.fun")(dat, list(integer(0)), rep(1L, n))
    sco <- attr(sf, "nh.scores.fun")(ne$op, ne$u, ne$v, pas, gs, NULL)
    if (!all(is.finite(sco)))
      next
    ## verify = TRUE errors if the bound is ever violated
    am <- attr(sf, "nh.argmax.fun")(ne$op, ne$u, ne$v, pas, gs, NULL, TRUE)
    stopifnot(identical(am$index, which.max(sco)),
              identical(am$total, max(sco)), am$worst < 1)
    worst <- max(worst, am$worst)
    nadv <- nadv + 1L
  }
}
stopifnot(nadv > 20L)
cat(sprintf("adversarial scales: %d neighbourhoods, bound held, worst ratio %.3e\n",
            nadv, worst))

################################################################################
## 5. degenerate neighbourhoods
################################################################################

## a single candidate: the band is that candidate
S <- setup(2, seed = 3L, dens = 0)
ne1 <- list(op = 1L, u = 1L, v = 2L)
gs <- idlBNs:::.iBIC.global.sufstats(S$dat)
am <- idlBNs:::.iBIC.nh.argmax(ne1$op, ne1$u, ne1$v, S$pas, gs, NULL, TRUE)
sco <- idlBNs:::.iBIC.nh.scores(ne1$op, ne1$u, ne1$v, S$pas, gs, NULL)
stopifnot(identical(am$index, 1L), identical(am$total, sco[1]), am$band == 1)

## an empty neighbourhood is an error, not a silent answer
stopifnot(inherits(tryCatch(idlBNs:::.iBIC.nh.argmax(integer(0), integer(0),
                                                     integer(0), S$pas, gs,
                                                     NULL, FALSE),
                            error = function(e) e), "error"))
cat("degenerate neighbourhoods handled\n")

cat("all argmax band tests passed\n")
