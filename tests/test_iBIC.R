## 2026-08-20 regression test for iBIC

suppressPackageStartupMessages({
  library(graph)
  library(idlBNs)
})

p <- 3
nobs <- 100
nint <- 100
n <- nobs + nint

## define a DAG structure of a Bayesian network with three vertices
## forming a Markov chain X1 -> X2 -> X3
g <- new("graphNEL", nodes=c("X1", "X2", "X3"), edgemode="directed")
g <- addEdge("X1", "X2", g)
g <- addEdge("X2", "X3", g)

## simulate observational data for the previous DAG X1 -> X2 -> X3
set.seed(123)
X1 <- rnorm(nobs, mean=0, sd=1)
X2 <- 0.5 * X1 + rnorm(nobs, mean=0, sd=1)
X3 <- 0.5 * X2 + rnorm(nobs, mean=0, sd=1)
obsdat <- data.frame(X1=X1, X2=X2, X3=X3)

## simulate interventional data for the same DAG, where X2 is intervened
X1 <- rnorm(nint, mean=0, sd=1)
X2 <- rnorm(nint, mean=0, sd=1) + 1.0
X3 <- 0.5 * X2 + rnorm(nint, mean=0, sd=1)
intdat <- data.frame(X1=X1, X2=X2, X3=X3)

## combine observational and interventional data
dat <- rbind(obsdat, intdat)

## define the targets and target indices for the interventional data
targets <- list(0L, 2L)
target.index <- c(rep(1L, nobs), rep(2L, nint))

## calculate the interventional BIC score for the DAG and data set
iBIC(g, dat, targets, target.index)

## create another Markov equivalent DAG by reversing the arc X1 -> X2
## to obtain X1 <- X2 -> X3
g2 <- g
g2 <- removeEdge("X1", "X2", g2)
g2 <- addEdge("X2", "X1", g2)

## calculate the interventional BIC score for the new DAG on the
## same data, notice that the score is different despite being a
## Markov equivalent DAG
iBIC(g2, dat, targets, target.index)

## this is not the case if we do not indicate the presence of interventions
## in the data
iBIC(g, dat)
iBIC(g2, dat)

## verify that the C engine and the R engine give numerically identical results

tol <- .Machine$double.eps^0.5

s_C  <- iBIC(g, dat, targets, target.index, engine="C")
s_R  <- iBIC(g, dat, targets, target.index, engine="R")
stopifnot(abs(s_C - s_R) < tol)

s2_C <- iBIC(g2, dat, targets, target.index, engine="C")
s2_R <- iBIC(g2, dat, targets, target.index, engine="R")
stopifnot(abs(s2_C - s2_R) < tol)

## pure observational data (no interventions)
sobs_C  <- iBIC(g,  dat, engine="C")
sobs_R  <- iBIC(g,  dat, engine="R")
stopifnot(abs(sobs_C - sobs_R) < tol)

sobs2_C <- iBIC(g2, dat, engine="C")
sobs2_R <- iBIC(g2, dat, engine="R")
stopifnot(abs(sobs2_C - sobs2_R) < tol)

## verify that .check_cached_scores()'s fast path (skip full validation on
## repeat calls with the same cached.scores object, marked internally after
## its first successful validation) does not change the calculated score,
## reusing the same cached.scores object across repeated calls, as is done
## by hcmc()/hillclimbing() during search

csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()),
                  simplify=FALSE)

s_first  <- iBIC(g, dat, targets, target.index, cached.scores=csco)
s_second <- iBIC(g, dat, targets, target.index, cached.scores=csco)
s_third  <- iBIC(g, dat, targets, target.index, cached.scores=csco)
stopifnot(abs(s_first - s_second) < tol)
stopifnot(abs(s_first - s_third) < tol)

## a different DAG, reusing the same (now fast-pathed) cached.scores object
s2_cached <- iBIC(g2, dat, targets, target.index, cached.scores=csco)
stopifnot(abs(s2_cached - s2_C) < tol)

## a malformed cached.scores must still be rejected: since it has never been
## validated before, it cannot have the internal '.validated' marker set, so
## the fast path must not accidentally bypass validation for it

bad_length <- replicate(numNodes(g) - 1L,
                        new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)
res <- tryCatch({
    iBIC(g, dat, targets, target.index, cached.scores=bad_length)
    "no error"
}, error=function(e) "error")
stopifnot(identical(res, "error"))

bad_elements <- as.list(seq_len(numNodes(g)))
res <- tryCatch({
    iBIC(g, dat, targets, target.index, cached.scores=bad_elements)
    "no error"
}, error=function(e) "error")
stopifnot(identical(res, "error"))

## verify that iBIC()'s C engine's whole-DAG-loop C function (which does
## cache lookup/write-back entirely in C, replicating .cached_scores_key()'s
## key format) is cross-compatible with the pure-R engine's cache handling:
## a cached.scores object populated by one engine must be readable, with
## identical scores, by the other engine

csco_R_populated <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()),
                              simplify=FALSE)
s_via_R <- iBIC(g, dat, targets, target.index,
                cached.scores=csco_R_populated, engine="R")
s_via_C_cached <- iBIC(g, dat, targets, target.index,
                       cached.scores=csco_R_populated, engine="C")
stopifnot(abs(s_via_R - s_via_C_cached) < tol)

csco_C_populated <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()),
                              simplify=FALSE)
s_via_C <- iBIC(g, dat, targets, target.index,
               cached.scores=csco_C_populated, engine="C")
s_via_R_cached <- iBIC(g, dat, targets, target.index,
                      cached.scores=csco_C_populated, engine="R")
stopifnot(abs(s_via_C - s_via_R_cached) < tol)

## the two cached.scores objects, populated independently by the C and R
## engines respectively, must end up with identical keys in every per-node
## environment (confirms the C and R cache-key formats are byte-identical,
## not merely that the final scores happen to match)
for (i in seq_len(numNodes(g)))
    stopifnot(identical(sort(ls(csco_R_populated[[i]])),
                        sort(ls(csco_C_populated[[i]]))))

## edge cases: a minimal 2-node DAG (one root with an empty parent set,
## exercising the literal ":" cache-key special case, plus one child), and
## repeated cache hits/misses through the C engine

g2n <- new("graphNEL", nodes=c("Y1", "Y2"), edgemode="directed")
g2n <- addEdge("Y1", "Y2", g2n)
set.seed(1)
dat2n <- data.frame(Y1=rnorm(50), Y2=rnorm(50))
s2n_C <- iBIC(g2n, dat2n, engine="C")
s2n_R <- iBIC(g2n, dat2n, engine="R")
stopifnot(abs(s2n_C - s2n_R) < tol)

csco_hitmiss <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()),
                          simplify=FALSE)
s_allmiss <- iBIC(g, dat, targets, target.index,
                  cached.scores=csco_hitmiss, engine="C")
s_allhit  <- iBIC(g, dat, targets, target.index,
                  cached.scores=csco_hitmiss, engine="C")
stopifnot(abs(s_allmiss - s_allhit) < tol)

## no caching at all (cached.scores=NULL) through the C engine, repeated
s_nocache1 <- iBIC(g, dat, targets, target.index, cached.scores=NULL, engine="C")
s_nocache2 <- iBIC(g, dat, targets, target.index, cached.scores=NULL, engine="C")
stopifnot(abs(s_nocache1 - s_nocache2) < tol)
