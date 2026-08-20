## 2026-08-20 regression test for iBGe

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
iBGe(g2, dat, targets, target.index)

## this is not the case if we do not indicate the presence of interventions
## in the data
iBGe(g, dat)
iBGe(g2, dat)
