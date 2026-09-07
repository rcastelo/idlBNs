## Shared harness for the Tier B golden trajectory reference.
##
## Trajectory bit-exactness rests on LAPACK rounding, which differs between
## BLAS implementations, so this reference is NOT part of the package test
## suite (regression/ is .Rbuildignore'd). It is a development gate: it
## catches drift that the in-suite C-vs-R differential tests cannot, namely a
## C and R implementation that agree with each other but have BOTH moved
## relative to the pre-port baseline. A change to the score cache's
## population order does exactly that.
##
## Run regression/make_reference.R once, on a known-good tree, to create
## regression/reference.csv. Run regression/check_reference.R at every gate.

suppressPackageStartupMessages({
    library(graph)
    library(pcalg)
    library(idlBNs)
})

## the grid: 48 runs
GOLDEN_P     <- c(10L, 20L, 30L, 40L)
GOLDEN_SCORE <- c("iBIC", "iBGe")
GOLDEN_ALG   <- c("hcmc", "hillclimbing")
GOLDEN_SEED  <- 1:3

## simulate one data set per p: a random DAG, a family of 3 intervention
## targets plus the empty target, and mixed observational/interventional
## Gaussian data. deterministic in 'seed' alone.
golden_data <- function(p, k=3L, n=200L, seed=123L, d=0.3) {
    set.seed(seed)
    Mg <- r.gauss.pardag(p, d, top.sort=TRUE, normalize=TRUE)
    I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))
    nbytgts <- rep(floor(n / (k + 1)), k)
    nbytgts <- c(n - sum(nbytgts), nbytgts)
    dat <- do.call(rbind,
                   lapply(seq_along(I),
                          function(v) rmvnorm.ivent(nbytgts[v], Mg,
                                                    target=I[[v]],
                                                    target.value=rep(2, length(I[[v]])))))

    list(dat=dat, targets=I, target.index=rep(seq_along(nbytgts), nbytgts))
}

## the learned DAG as a single canonical string: sorted "u>v" arcs, so the
## comparison is on the arc SET and is insensitive to edgeL storage order
## (which is an internal detail, unlike the arc set, which is the result)
golden_edges <- function(dag) {
    em <- edgeMatrix(dag)
    if (ncol(em) == 0L)
        return("")
    v <- nodes(dag)

    paste(sort(paste0(v[em["from", ]], ">", v[em["to", ]])), collapse=",")
}

## run the whole grid, returning a data.frame ready to write or compare.
## scores are formatted %.17g -- full double precision, so a one-ulp
## trajectory divergence shows up as a text difference
golden_run <- function(verbose=TRUE) {
    rows <- list()
    for (p in GOLDEN_P) {
        x <- golden_data(p)
        for (sf in GOLDEN_SCORE)
            for (alg in GOLDEN_ALG)
                for (sd in GOLDEN_SEED) {
                    set.seed(1000L + sd)
                    tm <- system.time(
                        res <- do.call(alg,
                                       list(x$dat, targets=x$targets,
                                            target.index=x$target.index,
                                            scorefun=get(sf), verbose=FALSE)))
                    rows[[length(rows) + 1L]] <-
                        data.frame(p=p, score=sf, alg=alg, seed=sd,
                                   sco=sprintf("%.17g", res$sco),
                                   nedges=numEdges(res$dag),
                                   edges=golden_edges(res$dag),
                                   stringsAsFactors=FALSE)
                    if (verbose)
                        cat(sprintf("  %-4s %-12s p=%-3d seed=%d  %7.2fs  sco=%s\n",
                                    sf, alg, p, sd, tm[["elapsed"]],
                                    sprintf("%.17g", res$sco)))
                }
    }

    do.call(rbind, rows)
}

GOLDEN_FILE <- file.path("regression", "reference.csv")
