## 2026-09-11 the population ("large-sample limit") input to the scores.
##
## A population model stands in for a data matrix: the same score arithmetic
## is fed the sufficient statistics a sample of that size would have IN
## EXPECTATION. Three things are pinned here.
##
##   1. The statistics are the ones analysis/01_population_score.R derives
##      independently, so the limit is the documented one and not whatever
##      this implementation happens to compute.
##   2. The finite-sample statistics CONVERGE to them, for both scores.
##   3. The whole search runs on a population model, through both engines
##      and both scores, and the object behaves like the matrix it replaces.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

P  <- 5L
set.seed(11)
G  <- pcalg::r.gauss.pardag(P, 0.5); G$set.intercept(runif(P, -1, 1))
TG <- list(integer(0), 2L, c(3L, 5L))
CNT <- c(600, 200, 200); IV <- 2; TV <- 0
N  <- sum(CNT)
PM <- population(G, ivent.value = IV, ivent.var = TV)

mkg <- function(A) { nn <- as.character(seq_len(P)); dimnames(A) <- list(nn, nn)
                     as(as(A * 1, "graphAM"), "graphNEL") }
set.seed(3)
DAGS <- replicate(10, { A <- matrix(FALSE, P, P); o <- sample(P)
    for (i in 1:(P-1)) for (j in (i+1):P) if (runif(1) < 0.45) A[o[i], o[j]] <- TRUE
    A }, simplify = FALSE)

################################################################################
## 1. the statistics are analysis/01_population_score.R's, independently derived
################################################################################

B <- G$weight.mat(); OM <- G$err.var(); B0 <- G$intercept()
int.mom <- function(I) { Bk <- B; omk <- OM; b0k <- B0
    if (length(I)) { Bk[, I] <- 0; omk[I] <- TV; b0k[I] <- IV }
    M <- solve(diag(P) - t(Bk))
    list(mu = as.vector(M %*% b0k), Sigma = M %*% diag(omk, P) %*% t(M)) }
pool <- function(envs, w) { w <- w / sum(w)
    mu <- Reduce(`+`, Map(function(e, wi) wi * e$mu, envs, w))
    list(mu = mu, Sigma = Reduce(`+`, Map(function(e, wi)
             wi * (e$Sigma + tcrossprod(e$mu - mu)), envs, w))) }
resid.var <- function(S, j, pa) if (!length(pa)) S[j, j] else
    max(S[j, j] - drop(S[j, pa, drop = FALSE] %*%
                       solve(S[pa, pa, drop = FALSE], S[pa, j])), .Machine$double.eps)
ss <- lapply(seq_len(P), function(j) {
    keep <- which(!vapply(TG, function(I) j %in% I, TRUE)); w <- CNT / sum(CNT)
    list(W = sum(w[keep]), Sigma = pool(lapply(TG[keep], int.mom), w[keep])$Sigma) })
pop.score <- function(A) { fit <- 0; dim <- 0
    for (j in seq_len(P)) { pa <- which(A[, j])
        fit <- fit - 0.5 * ss[[j]]$W * (1 + log(resid.var(ss[[j]]$Sigma, j, pa)))
        dim <- dim + 2 + length(pa) }
    N * fit - 0.5 * log(N) * dim }

## iBIC on the population model differs from it only by a constant that does
## not depend on the DAG -- 0.5*log(n)*p, the p intercept terms pop.score()
## counts in its dimension and iBIC does not -- so the two RANK DAGs alike
dummy <- matrix(0, 2, P, dimnames = list(NULL, as.character(seq_len(P))))
attr(dummy, "sanitycheck") <- TRUE
gsf <- idlBNs:::.iBIC.population.sufstats(PM, TG, CNT)
offs <- vapply(DAGS, function(A) {
    pas <- lapply(seq_len(P), function(j) as.integer(which(A[, j])))
    iBIC(mkg(A), dummy, TG, NULL, global.sufstats = gsf, pasets = pas) - pop.score(A)
}, 0)
stopifnot(max(abs(offs - 0.5 * log(N) * P)) < 1e-8,   # the predicted constant
          diff(range(offs)) < 1e-8)                   # and it IS constant
cat("test_population.R: statistics match 01_population_score.R exactly\n")

################################################################################
## 2. the finite-sample statistics converge to them, for BOTH scores
################################################################################

rel <- function(n) {
    cnt <- round(n * CNT / sum(CNT))
    set.seed(42)
    dat <- do.call(rbind, Map(function(I, m)
        rmvnorm.ivent(m, G, target = I, target.value = rep(IV, length(I))), TG, cnt))
    colnames(dat) <- as.character(seq_len(P)); ti <- rep(seq_along(cnt), cnt)
    attr(dat, "sanitycheck") <- TRUE
    fb <- idlBNs:::.iBIC.global.sufstats(dat, TG, ti)
    pb <- idlBNs:::.iBIC.population.sufstats(PM, TG, cnt)
    fg <- idlBNs:::.iBGe.global.sufstats(dat, TG, ti)
    pg <- idlBNs:::.iBGe.population.sufstats(PM, TG, cnt)
    v <- vapply(DAGS[1:4], function(A) {
        pas <- lapply(seq_len(P), function(j) as.integer(which(A[, j]))); g <- mkg(A)
        c(abs(iBIC(g, dat, TG, ti, global.sufstats = fb, pasets = pas) -
              iBIC(g, dat, TG, ti, global.sufstats = pb, pasets = pas)) /
          abs(iBIC(g, dat, TG, ti, global.sufstats = pb, pasets = pas)),
          abs(iBGe(g, dat, TG, ti, global.sufstats = fg, pasets = pas) -
              iBGe(g, dat, TG, ti, global.sufstats = pg, pasets = pas)) /
          abs(iBGe(g, dat, TG, ti, global.sufstats = pg, pasets = pas)))
    }, numeric(2))
    apply(v, 1, max)
}
r4 <- rel(1e4); r6 <- rel(1e6)
stopifnot(r6[1] < r4[1] / 5, r6[2] < r4[2] / 5,   # shrinking, both scores
          r6[1] < 5e-3, r6[2] < 5e-3)
cat(sprintf("test_population.R: finite-sample converges (iBIC %.1e -> %.1e, iBGe %.1e -> %.1e)\n",
            r4[1], r6[1], r4[2], r6[2]))

################################################################################
## 3. the object stands in for the data matrix, and the search runs on it
################################################################################

stopifnot(ncol(PM) == P, identical(colnames(PM), as.character(seq_len(P))),
          is.na(nrow(PM)))
## a bare GaussParDAG is accepted too, as a hard intervention to 0
stopifnot(inherits(idlBNs:::.as.population(G), "idlBNsPopulation"),
          idlBNs:::.as.population(G)$ivent.value == 0)
## target.index must now be one count per environment
stopifnot(inherits(tryCatch(hcmc(PM, targets = TG, target.index = rep(1L, 10),
                                 verbose = FALSE), error = function(e) e), "error"))

## Omitting target.index. The formal default used to be rep(1L, nrow(x)),
## which a population object cannot satisfy -- nrow() is NA by construction --
## so it died inside rep() with "invalid 'times' argument" before any check
## could say anything useful. There is no defensible default count either,
## because the notional sample size sets the score's penalty and so decides
## between nested models; it has to be stated, on the object or in the call.
for (f in list(hcmc, hillclimbing)) {
    e <- tryCatch(f(PM, targets = TG, verbose = FALSE), error = function(e) e)
    stopifnot(inherits(e, "error"),
              grepl("sample size", conditionMessage(e)))   # and says why
}
e <- tryCatch(iBIC(mkg(DAGS[[1]]), PM, targets = TG), error = function(e) e)
stopifnot(inherits(e, "error"), grepl("sample size", conditionMessage(e)))

## with n on the object it is split equally, and explicit counts still win
PMn <- population(G, n = 4 * sum(CNT), ivent.value = IV, ivent.var = TV)
set.seed(5); a <- hcmc(PMn, targets = TG, verbose = FALSE)
set.seed(5); b <- hcmc(PMn, targets = TG, verbose = FALSE,
                       target.index = rep(4 * sum(CNT) / length(TG), length(TG)))
set.seed(5); d <- hcmc(PMn, targets = TG, target.index = CNT, verbose = FALSE)
stopifnot(identical(a$sco, b$sco), !identical(a$sco, d$sco))
## C: how much data an observational environment carries relative to each
## interventional one, the quantity Wang, Solus, Yang and Uhler's
## counterexample turns on. n is split in proportion to (C, 1, ..., 1), which
## is the weighting analysis/33, /43 and /45 use, so the population path
## reproduces their convention exactly. (The finite-sample allocator in
## analysis/34 and /44 rounds the same proportions to whole rows; here there
## are no rows, so the split stays real-valued.)
ti <- function(k, Cv, n = 8000) {
    tg <- if (k == 0) list(integer(0))
          else c(list(integer(0)), lapply(seq_len(k), function(i) as.integer(i)))
    idlBNs:::.resolve.target.index(population(G, n = n, C = Cv), tg, NULL)
}
stopifnot(isTRUE(all.equal(ti(3, 5), 8000 * c(5, 1, 1, 1) / 8)),
          isTRUE(all.equal(ti(1, 50), 8000 * c(50, 1) / 51)),
          isTRUE(all.equal(ti(3, 1), rep(2000, 4))),   # C = 1 is the equal split
          isTRUE(all.equal(ti(0, 50), 8000)))          # inert with no interventions
## an environment is observational because its target is empty, not because of
## where it sits in the list
a <- idlBNs:::.resolve.target.index(population(G, n = 900, C = 7),
                                    list(integer(0), 2L, 3L), NULL)
b <- idlBNs:::.resolve.target.index(population(G, n = 900, C = 7),
                                    list(2L, integer(0), 3L), NULL)
stopifnot(identical(sort(a), sort(b)), which.max(a) == 1L, which.max(b) == 2L)
## C only proportions; without n there is still nothing to scale
stopifnot(inherits(tryCatch(hcmc(population(G, C = 5), targets = TG,
                                 verbose = FALSE), error = function(e) e), "error"),
          inherits(tryCatch(population(G, C = 0), error = function(e) e), "error"))
## and explicit counts still win over C
set.seed(5); u <- hcmc(population(G, n = 1000, C = 5, ivent.value = IV),
                       targets = TG, verbose = FALSE)
set.seed(5); v <- hcmc(population(G, n = 1000, C = 5, ivent.value = IV),
                       targets = TG, target.index = c(400, 300, 300), verbose = FALSE)
stopifnot(!identical(u$sco, v$sco))

## Zero-count environments. A count of 0 says an environment contributes
## nothing, which is legitimate -- but it stayed in the pooling weights, and a
## vertex left alone ONLY by zero-count environments then divided 0 by 0. The
## result was not an error: iBGe returned NaN and the search died later on
## "missing value where TRUE/FALSE needed", while iBIC failed inside the
## Cholesky with "dpotrf failed". Worse, a FRACTIONAL count never went near
## either path and came back as a finite, plausible score for a vertex with
## less than one observation's worth of mass, where the finite-sample path
## refuses anything under two rows.
local({
    G0 <- mkg(matrix(FALSE, P, P))
    tg <- list(integer(0), 1L)              # only env 1 leaves variable 1 alone
    for (cs in list(c(0, 500), c(0.4, 500), c(1.9, 500)))
        for (sf in list(iBIC, iBGe)) {
            e <- tryCatch(sf(G0, PM, targets = tg, target.index = cs),
                          error = function(e) e)
            stopifnot(inherits(e, "error"),
                      grepl("variable 1", conditionMessage(e)))   # and says which
        }
    ## the boundary is where the finite-sample path puts it
    stopifnot(is.finite(iBIC(G0, PM, targets = tg, target.index = c(2, 500))))
    ## a zero-count environment that is not the only one leaving a vertex alone
    ## is simply ignored -- dropping it must not change any score
    stopifnot(isTRUE(all.equal(
        iBIC(G0, PM, targets = list(integer(0), 1L, 2L), target.index = c(400, 100, 0)),
        iBIC(G0, PM, targets = list(integer(0), 1L),     target.index = c(400, 100)))))
    ## and a whole search still runs with one present
    stopifnot(is.finite(hcmc(PM, targets = list(integer(0), 1L, 2L),
                             target.index = c(4000, 1000, 0), verbose = FALSE)$sco))
})

## and the data path is untouched by the change of default
XX <- matrix(rnorm(120 * P), 120, P, dimnames = list(NULL, as.character(seq_len(P))))
stopifnot(is.finite(hcmc(XX, verbose = FALSE)$sco),
          is.finite(hillclimbing(XX, verbose = FALSE)$sco))

truth <- as(dag2essgraph(mkg(G$weight.mat() != 0), targets = TG), "graphNEL")
for (sf in list(iBIC, iBGe)) for (eng in c("C", "R")) {
    set.seed(5)
    h <- hcmc(PM, targets = TG, target.index = CNT, scorefun = sf,
              verbose = FALSE, engine = eng, sampler = "exact",
              escape = "exhaustive")
    stopifnot(is.finite(h$sco), h$sampler.fallbacks == 0L)
    if (identical(sf, iBIC))                      # the limit recovers the truth
        stopifnot(shd(as(dag2essgraph(h$dag, targets = TG), "graphNEL"), truth) == 0)
}
hc <- hillclimbing(PM, targets = TG, target.index = CNT, verbose = FALSE)
stopifnot(is.finite(hc$sco))
cat("test_population.R: hcmc() and hillclimbing() run in the population limit\n")
cat("all population-input tests passed\n")
