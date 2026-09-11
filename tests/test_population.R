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
