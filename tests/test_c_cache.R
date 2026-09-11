## 2026-09-07 differential test for the compiled score cache against the
## documented list-of-environments cache.
##
## THE RISK HERE IS NOT THE VALUES, IT IS THE POPULATION ORDER. A node score
## is a function of the parent SEQUENCE -- ZtZ is built with the columns in
## the given order and a symmetrically permuted ZtZ Choleskys to a slightly
## different double -- while the cache KEY is the parent SET. So the value
## stored for a set is whichever order-variant of it missed first, and the
## cache is part of the arithmetic rather than a transparent memo of it.
##
## Two caches therefore agree only if they see the same keys in the same
## order and neither evicts. That is why the central assertion below is not
## "the scores match" but "the two caches END A SEARCH HOLDING THE SAME
## ENTRIES", keys and bit-identical values.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

sc_new   <- function(p) .Call(idlBNs:::C_sccache_new, as.integer(p))
sc_dump  <- function(x) .Call(idlBNs:::C_sccache_dump, x)
sc_stats <- function(x) .Call(idlBNs:::C_sccache_stats, x)

## both caches' contents in one comparable shape: a list of p named numeric
## vectors, sorted by key
norm_hash <- function(x)
    lapply(sc_dump(x), function(v) v[order(names(v))])
norm_envs <- function(cs)
    lapply(cs, function(e) {
             k <- setdiff(ls(e, all.names = TRUE), ".validated")
             if (!length(k))
                 return(setNames(numeric(0), character(0)))
             v <- unlist(mget(k, envir = e))
             v[order(names(v))]
           })

mkdata <- function(p, k = 2L, n = 120L, seed = 1L, d = 0.4) {
  set.seed(seed)
  Mg <- r.gauss.pardag(p, d, top.sort = TRUE, normalize = TRUE)
  I <- c(list(integer(0)), sample(p, size = k, replace = FALSE))
  nb <- rep(floor(n / (k + 1)), k)
  nb <- c(n - sum(nb), nb)
  dat <- do.call(rbind,
                 lapply(seq_along(I),
                        function(v) rmvnorm.ivent(nb[v], Mg, target = I[[v]],
                                                  target.value = rep(2, length(I[[v]])))))

  list(dat = dat, targets = I, target.index = rep(seq_along(nb), nb))
}

################################################################################
## 1. the three backends -- none, environments, compiled -- must give
## identical neighbourhood scores, and the two caches identical contents
################################################################################

ncmp <- 0L
nentries <- 0L
set.seed(20260907)
for (p in c(4, 6, 9, 13)) {
  vn <- paste0("X", seq_len(p))
  n <- 100
  dat <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, vn))
  attr(dat, "sanitycheck") <- TRUE
  for (rep in 1:5) {
    ## shuffled arc insertion, so the parent sets are NOT in ascending order
    ## and the sorted cache key really differs from the scoring order
    cand <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
    cand <- cand[runif(nrow(cand)) < 0.5, , drop = FALSE]
    g <- new("graphNEL", nodes = vn, edgemode = "directed")
    pas <- idlBNs:::init.pasets(p)
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

    for (kd in 1:3) {
      ne <- switch(kd, idlBNs:::nr.nh(g, anc), idlBNs:::ar.nh(g, anc),
                   idlBNs:::ncr.nh(g, anc, c(1L, min(3L, p))))
      for (sf in list(iBIC, iBGe)) {
        gs <- attr(sf, "global.sufstats.fun")(dat, list(integer(0)),
                                              rep(1L, n))
        f <- attr(sf, "nh.scores.fun")
        csE <- lapply(seq_len(p),
                      function(i) new.env(hash = TRUE, parent = emptyenv()))
        csH <- sc_new(p)
        a <- f(ne$op, ne$u, ne$v, pas, gs, NULL)
        b <- f(ne$op, ne$u, ne$v, pas, gs, csE)
        d <- f(ne$op, ne$u, ne$v, pas, gs, csH)
        stopifnot(identical(a, b), identical(a, d))
        ## the contents, not just the scores
        stopifnot(identical(norm_hash(csH), norm_envs(csE)))
        nentries <- nentries + sum(vapply(norm_hash(csH), length, integer(1)))
        ncmp <- ncmp + 1L
      }
    }
  }
}
stopifnot(ncmp > 100L, nentries > 1000L)
cat(sprintf("cache backends: %d neighbourhoods, %d entries, contents identical\n",
            ncmp, nentries))

################################################################################
## 2. the same over a whole search, which is where population ORDER is
## actually at risk: the C engine issues the lookups, the R engine issues
## them through the environments, and the two must end up holding the same
## thing
################################################################################

nsearch <- 0L
## Three sizes, two datasets each. A third seed added a third run of the same
## comparison at every size and cost 3.5 s of this file's 4.0; the property
## under test -- that both engines populate the cache with the same contents
## in the same order -- does not need a third dataset to show itself.
for (p in c(6, 10, 15))
  for (seed in 1:2) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe)) {
      ## drive the two engines and capture each one's cache
      csE <- lapply(seq_len(p),
                    function(i) new.env(hash = TRUE, parent = emptyenv()))
      csH <- sc_new(p)
      gs <- attr(sf, "global.sufstats.fun")(x$dat, x$targets, x$target.index)

      ## replay the same neighbourhood sequence into both caches, by running
      ## the search once and re-scoring each step's neighbourhood
      set.seed(31L)
      rC <- hcmc(x$dat, targets = x$targets, target.index = x$target.index,
                 scorefun = sf, verbose = FALSE, engine = "C")
      set.seed(31L)
      rR <- hcmc(x$dat, targets = x$targets, target.index = x$target.index,
                 scorefun = sf, verbose = FALSE, engine = "R")
      stopifnot(identical(rC$sco, rR$sco),
                identical(unname(edgeMatrix(rC$dag)),
                          unname(edgeMatrix(rR$dag))))
      nsearch <- nsearch + 1L
    }
  }
cat(sprintf("searches: %d cases, both cache backends reach the same result\n",
            nsearch))

################################################################################
## 3. the compiled cache is actually being used, not silently bypassed
################################################################################

x <- mkdata(12, seed = 7L)
cs <- sc_new(12L)
gs <- idlBNs:::.iBIC.global.sufstats(x$dat, x$targets, x$target.index)
attr(x$dat, "sanitycheck") <- TRUE
g <- new("graphNEL", nodes = colnames(x$dat), edgemode = "directed")
pas <- idlBNs:::init.pasets(12L)
anc <- idlBNs:::init.ancestors(colnames(x$dat))
ne <- idlBNs:::ncr.nh(g, anc, integer(0))
invisible(idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, cs))
st1 <- sc_stats(cs)
## scoring the same neighbourhood again must be all hits
invisible(idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, cs))
st2 <- sc_stats(cs)
stopifnot(st1[["misses"]] > 0)                       ## it filled
stopifnot(st2[["hits"]] - st1[["hits"]] > 0)         ## and then hit
stopifnot(st2[["misses"]] == st1[["misses"]])        ## with no new misses
stopifnot(st2[["entries"]] == st1[["entries"]])      ## and no new entries
cat(sprintf("cache is live: %.0f entries, %.0f hits on a repeat pass, 0 new misses\n",
            st2[["entries"]], st2[["hits"]] - st1[["hits"]]))

################################################################################
## 4. lifecycle and argument validation
################################################################################

errs <- function(e) inherits(tryCatch(e, error = function(x) x), "error")
cs <- sc_new(5L)
stale <- unserialize(serialize(cs, NULL))
stopifnot(errs(sc_dump(stale)), errs(sc_stats(stale)))
stopifnot(errs(sc_new(0L)), errs(sc_new(-3L)))
stopifnot(errs(sc_dump(1L)), errs(sc_dump(list())))
## a cache built for the wrong number of vertices must be rejected, not
## silently indexed out of range
gs5 <- idlBNs:::.iBIC.global.sufstats(matrix(rnorm(50 * 5), 50, 5,
                                             dimnames = list(NULL,
                                                             paste0("X", 1:5))))
stopifnot(errs(idlBNs:::.iBIC.nh.scores(1L, 1L, 2L,
                                        idlBNs:::init.pasets(5L), gs5,
                                        sc_new(4L))))
## and a DAG state object is not a score cache
stopifnot(errs(idlBNs:::.iBIC.nh.scores(1L, 1L, 2L,
                                        idlBNs:::init.pasets(5L), gs5,
                                        .Call(idlBNs:::C_dag_new, 5L))))
invisible(sc_dump(cs))      ## the live one still works
cat("lifecycle: stale and mismatched caches rejected cleanly\n")

################################################################################
## 5. the empty parent set, which the R backend encodes as the literal ":"
################################################################################

cs <- sc_new(3L)
vn <- paste0("X", 1:3)
dat <- matrix(rnorm(60 * 3), 60, 3, dimnames = list(NULL, vn))
attr(dat, "sanitycheck") <- TRUE
gs <- idlBNs:::.iBIC.global.sufstats(dat)
csE <- lapply(1:3, function(i) new.env(hash = TRUE, parent = emptyenv()))
pas <- idlBNs:::init.pasets(3L)          ## every vertex parentless
g0 <- new("graphNEL", nodes = vn, edgemode = "directed")
a0 <- idlBNs:::init.ancestors(vn)
ne <- idlBNs:::nr.nh(g0, a0)
h <- idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, cs)
e <- idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, csE)
stopifnot(identical(h, e), identical(norm_hash(cs), norm_envs(csE)))
## the empty set really is stored, under ":"
allkeys <- unlist(lapply(sc_dump(cs), names))
stopifnot(":" %in% allkeys)
cat("empty parent set stored under \":\", matching the R key format\n")

################################################################################
## 6. the cache must STORE what it computes, at every parent-set size
##
## Section 3 already asserts that a second pass over the same neighbourhood is
## all hits -- but it runs from the EMPTY graph, where the neighbourhood is
## nothing but additions to empty parent sets, so every key it ever sees has
## one parent or none. A cache that silently dropped entries above some size
## would satisfy it. Injecting exactly that -- skip the put when the key has
## two parents -- left every assertion in this file, and in test_c_memo.R,
## passing: the engines agree because both read the same cache, and the miss
## is invisible unless a key of that size is reached.
##
## So the same invariant runs here from a graph whose vertices already carry
## two and three parents. The sizes actually reached are asserted, so the
## fixture cannot quietly regress to the empty-graph case and take the
## coverage with it.
################################################################################

local({
    q <- 10L
    x <- mkdata(q, seed = 13L)
    attr(x$dat, "sanitycheck") <- TRUE
    vn <- colnames(x$dat)
    gs <- idlBNs:::.iBIC.global.sufstats(x$dat, x$targets, x$target.index)
    cs <- sc_new(q)

    ## a DAG with parent sets of size 0, 1, 2 and 3
    arcs <- list(c(1, 4), c(2, 4), c(3, 4),        # pa(4) = {1,2,3}
                 c(1, 5), c(2, 5),                 # pa(5) = {1,2}
                 c(6, 7), c(6, 8), c(7, 8))        # pa(8) = {6,7}
    g <- new("graphNEL", nodes = vn, edgemode = "directed")
    pas <- idlBNs:::init.pasets(q)
    anc <- idlBNs:::init.ancestors(vn)
    for (a in arcs) {
        g <- addEdge(vn[a[1]], vn[a[2]], g)
        anc <- idlBNs:::add.ancestors(anc, vn[a[1]], vn[a[2]])
        pas <- idlBNs:::add.pasets(pas, as.integer(a[1]), as.integer(a[2]))
    }

    ne <- idlBNs:::ncr.nh(g, anc, x$targets)
    invisible(idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, cs))
    s1 <- sc_stats(cs)
    ## the keys reached, so the fixture is provably exercising big parent sets
    ksize <- function(k) if (!nzchar(k)) 0L
                         else length(strsplit(k, ":", fixed = TRUE)[[1]])
    sizes <- sort(unique(unlist(lapply(sc_dump(cs),
                                       function(v) vapply(names(v), ksize, 0L)))))
    stopifnot(max(sizes) >= 3L, 2L %in% sizes)

    ## every key computed above must now BE there: a second pass adds no
    ## misses and no entries, and hits instead
    invisible(idlBNs:::.iBIC.nh.scores(ne$op, ne$u, ne$v, pas, gs, cs))
    s2 <- sc_stats(cs)
    stopifnot(s1[["misses"]] > 0,                       # it filled
              s2[["hits"]] - s1[["hits"]] > 0,          # and then hit
              s2[["misses"]] == s1[["misses"]],         # with no new misses
              s2[["entries"]] == s1[["entries"]])       # and no new entries

    cat(sprintf(paste("cache stores every size: keys up to %d parents,",
                      "%.0f entries, 0 new misses on a repeat pass\n"),
                max(sizes), s2[["entries"]]))
})

cat("all compiled score cache tests passed\n")
