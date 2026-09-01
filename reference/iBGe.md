# BGe score for interventional Gaussian data

Score the goodness-of-fit (GoF) of a given structure of a Bayesian
network given an interventional data set of continuous values, where
observations are assumed to be independent but not identically
distributed (not iid) multivariate Gaussian. This GoF score corresponds
to the interventional Bayesian Gaussian equivalent (iBGe) score defined
by Kuipers and Moffa (2025). By default, the arguments `targets` and
`target.index` are set so that the calculated BIC score assumes there
are no interventions in the data.

## Usage

``` r
iBGe(
  g,
  dat,
  targets = list(integer(0)),
  target.index = rep(1L, nrow(dat)),
  cached.scores = NULL,
  global.sufstats = NULL
)
```

## Arguments

- g:

  An acyclic directed graph (DAG) structure of the Bayesian network for
  which we want to calculate the score.

- dat:

  A `data.frame` object with data records in the rows.

- targets:

  (Default `list(integer(0))`) A `list` object with a family of targets
  provided as a list of integer vectors. Its default value indicates
  that there are no interventions in the data, i.e., the data is purely
  observational.

- target.index:

  (Default a unit vector) A vector of integers in one-to-one
  correspondence with the rows in `dat`, indicating which rows in the
  input data are intervened by which targets. Its default value
  indicates that there are no interventions in the data, i.e., the data
  is purely observational.

- cached.scores:

  An optional list of environment objects, containing cached scores per
  parent set for each vertex in `g`. If `NULL` (default), no cached
  scores are used. Using this argument can speed up the calculation of
  the score when the same parent sets are scored multiple times. To use
  this argument, first create an empty environment object with
  `csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)`
  and then pass it to this `cached.scores` parameter, i.e.,
  `cached.scores=csco`. This is currently not implemented for the iBGe
  score, but it is included as an API placeholder for future versions of
  the package that will enable this feature for the iBGe score.

- global.sufstats:

  (Default `NULL`) An optional list of global sufficient statistics for
  the iBGe score, as returned by the `.iBGe.global.sufstats()` function,
  which do not depend on the structure of a specific DAG, but only on
  the input data (`dat`), the target vertices (`targets`) and the target
  indices (`target.index`) of the interventions. If `NULL` (default),
  the `.iBGe.global.sufstats()` function is internally called. This is
  currently not implemented for the iBGe score, but it is included as an
  API placeholder for future versions of the package that will enable
  this feature for the iBGe score.

## Value

A single numeric value corresponding to the interventional BGe score of
the given structure of the Bayesian network for the given data set.

## References

Kuipers, J. and Moffa, G. The interventional Bayesian Gaussian
equivalent score for Bayesian causal inference with unknown soft
interventions. *Proceedings of the Fourth Conference on Causal Learning
and Reasoning (PMLR)*, 275:772-791, 2025.

## Examples

``` r

library(graph)

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
targets <- list(integer(0), 2L)
target.index <- c(rep(1L, nobs), rep(2L, nint))

## calculate the interventional BGe score for the DAG and data set
iBGe(g, dat, targets, target.index)
#> [1] -713.708

## create another Markov equivalent DAG by reversing the arc X1 -> X2
## to obtain X1 <- X2 -> X3
g2 <- g
g2 <- removeEdge("X1", "X2", g2)
g2 <- addEdge("X2", "X1", g2)

## calculate the interventional BGe score for the new DAG on the
## same data, notice that the score is different despite being a
## Markov equivalent DAG
iBGe(g2, dat, targets, target.index)
#> [1] -721.2087

## this is not the case if we do not indicate the presence of interventions
## in the data
iBGe(g, dat)
#> [1] -891.4243
iBGe(g2, dat)
#> [1] -891.4243
```
