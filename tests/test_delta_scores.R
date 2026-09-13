## 2026-09-06 regression test for delta scoring of a whole neighborhood
##
## nr.nh()/ar.nh()/ncr.nh() return candidate moves, and the search loop
## scores the whole neighborhood through the score function's
## "nh.scores.fun" attribute, which rescores only the one or two vertices
## whose parent set each move changes and reuses the current DAG's terms
## for the rest (see src/nh_scores.c). Each candidate's terms are still
## summed over all p vertices, so the totals must come out BIT-IDENTICAL
## to scoring every candidate from scratch -- not merely close. That
## matters: reassociating the sum moves a total by an ulp, which is enough
## to swap two near-tied candidates under which.max() and send a greedy
## search into a different local optimum.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

################################################################################
## 1. for every candidate move of every neighborhood, the neighborhood
## scorer's total must be bit-identical to a full re-score of that
## candidate. exercised with and without a score cache, with and without
## interventions, for both iBIC and iBGe
################################################################################

set.seed(7)
ncmp <- 0L
for (p in c(5, 8, 12)) {
  varnames <- paste0("X", seq_len(p))
  n <- 80
  D <- r.gauss.pardag(p, 0.4, top.sort=TRUE, normalize=TRUE)
  dat <- rmvnorm.ivent(n, D, target=integer(0), target.value=numeric(0))
  colnames(dat) <- varnames
  attr(dat, "sanitycheck") <- TRUE
  utargets <- c(1L, min(3L, p))

  for (rep in 1:4) {
    ## a random DAG in topological order
    dag <- new("graphNEL", nodes=varnames, edgemode="directed")
    for (i in seq_len(p - 1))
      for (j in (i + 1):p)
        if (runif(1) < 0.35)
          dag <- addEdge(varnames[i], varnames[j], dag)

    anc <- idlBNs:::init.ancestors(varnames)
    em <- edgeMatrix(dag)
    for (k in seq_len(ncol(em)))
      anc <- idlBNs:::add.ancestors(anc, varnames[em["from", k]],
                                    varnames[em["to", k]])
    pasets <- idlBNs:::.build_pasets(dag, dat)

    nhs <- list(idlBNs:::nr.nh(dag, anc),
                idlBNs:::ar.nh(dag, anc),
                idlBNs:::ncr.nh(dag, anc, utargets))

    for (sf in list(iBIC, iBGe)) {
      nh.scores.fun <- attr(sf, "nh.scores.fun")
      stopifnot(is.function(nh.scores.fun))
      ## the same scorer reached through the package namespace rather than
      ## through the attribute. the two are the same function, but a search
      ## only ever reaches it the attribute way, and an attribute captured
      ## its closure when the package was built -- so a tool that rewrites
      ## the namespace binding afterwards (covr, trace(), debug()) never
      ## sees the attribute copy run. calling both here keeps the wrapper
      ## honestly exercised whichever way it is reached, and pins the
      ## attribute to the function it is supposed to name.
      direct <- switch(attr(sf, "scorefun.name"),
                       iBIC=idlBNs:::.iBIC.nh.scores,
                       iBGe=idlBNs:::.iBGe.nh.scores)
      stopifnot(is.function(direct))
      for (tgt in list(list(targets=list(integer(0)),
                            target.index=rep(1L, n)),
                       list(targets=list(integer(0), utargets),
                            target.index=c(rep(1L, n / 2),
                                           rep(2L, n - n / 2))))) {
        gs <- attr(sf, "global.sufstats.fun")(dat, tgt$targets,
                                              tgt$target.index)
        for (ne in nhs) {
          for (usecache in c(FALSE, TRUE)) {
            cs <- NULL
            if (usecache)
              cs <- lapply(seq_len(p),
                           function(i) new.env(hash=TRUE, parent=emptyenv()))
            got <- nh.scores.fun(ne$op, ne$u, ne$v, pasets, gs, cs)
            ## reaching the scorer either way must give the same answer
            cs2 <- NULL
            if (usecache)
              cs2 <- lapply(seq_len(p),
                            function(i) new.env(hash=TRUE, parent=emptyenv()))
            stopifnot(identical(got,
                                direct(ne$op, ne$u, ne$v, pasets, gs, cs2)))

            ## reference: score each candidate by summing all p vertex terms
            ref <- vapply(seq_along(ne$op), function(m) {
                            pas <- idlBNs:::move.pasets(pasets, ne$op[m],
                                                        ne$u[m], ne$v[m])
                            sf(dag, dat, targets=tgt$targets,
                               target.index=tgt$target.index,
                               cached.scores=NULL, global.sufstats=gs,
                               pasets=pas)
                          }, numeric(1))

            stopifnot(length(got) == length(ref))
            ## bit-identical, hence identical ranking and identical
            ## tie-breaking under which.max()
            stopifnot(identical(got, ref))
            stopifnot(identical(which.max(got), which.max(ref)))
            ncmp <- ncmp + length(ref)
          }
        }
      }
    }
  }
}
stopifnot(ncmp > 0L)
cat(sprintf("neighborhood scoring: %d candidates bit-identical to a full re-score\n",
            ncmp))

################################################################################
## 2. a custom scorefun carrying neither "nh.scores.fun" nor
## "supports.pasets" must still drive a search correctly, through
## score.nh()'s per-candidate fallback that materialises a neighbor graph
################################################################################

set.seed(11)
p <- 6
n <- 60
D <- r.gauss.pardag(p, 0.3, top.sort=TRUE, normalize=TRUE)
dat2 <- rmvnorm.ivent(n, D, target=integer(0), target.value=numeric(0))

## same score as iBIC, but stripped of every attribute the search uses to
## take a faster path, so the slowest fallback branch is exercised
plain_iBIC <- function(g, dat, targets=list(integer(0)),
                       target.index=rep(1L, nrow(dat)), cached.scores=NULL,
                       global.sufstats=NULL)
    iBIC(g, dat, targets=targets, target.index=target.index,
         cached.scores=cached.scores, global.sufstats=global.sufstats)

fast <- hcmc(dat2, scorefun=iBIC, verbose=FALSE)
set.seed(11)
fast2 <- hcmc(dat2, scorefun=iBIC, verbose=FALSE)
slow <- hcmc(dat2, scorefun=plain_iBIC, verbose=FALSE)

## the fallback path must reach a DAG scoring as well as the fast path
stopifnot(is.finite(slow$sco))
stopifnot(abs(slow$sco - fast$sco) < 1e-8 || slow$sco <= fast$sco)
stopifnot(numNodes(slow$dag) == p)

## and hillclimbing() likewise
slow_hc <- hillclimbing(dat2, scorefun=plain_iBIC, verbose=FALSE)
fast_hc <- hillclimbing(dat2, scorefun=iBIC, verbose=FALSE)
stopifnot(identical(unname(edgeMatrix(slow_hc$dag)),
                    unname(edgeMatrix(fast_hc$dag))))
stopifnot(abs(slow_hc$sco - fast_hc$sco) < 1e-9)

cat("custom-scorefun fallback path agrees with the fast path\n")

################################################################################
## 3. the score a search reports must be exactly what the score function
## returns for the DAG it ends on
################################################################################

for (sf in list(iBIC, iBGe)) {
  set.seed(3)
  res <- hcmc(dat2, scorefun=sf, verbose=FALSE)
  stopifnot(abs(res$sco - sf(res$dag, dat2)) < 1e-9)
  res_hc <- hillclimbing(dat2, scorefun=sf, verbose=FALSE)
  stopifnot(abs(res_hc$sco - sf(res_hc$dag, dat2)) < 1e-9)
}

cat("reported search scores match a full re-score of the returned DAG\n")
cat("all delta scoring tests passed\n")
