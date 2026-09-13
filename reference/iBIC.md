# BIC score for observational and interventional Gaussian data

Score the goodness-of-fit (GoF) of a given structure of a Bayesian
network given an interventional data set of continuous values, where
observations are assumed to be independent but not identically
distributed (not iid) multivariate Gaussian. This GoF score corresponds
to the Bayesian information criterion (BIC) as implemented in the
`GaussL0penIntScore` class from the `pcalg` package (Kalisch et al.,
2012). By default, the arguments `targets` and `target.index` are set so
that the calculated BIC score assumes there are no interventions in the
data.

## Usage

``` r
iBIC(
  g,
  x,
  targets = list(integer(0)),
  target.index = NULL,
  cached.scores = NULL,
  global.sufstats = NULL,
  engine = c("C", "R"),
  pasets = NULL
)
```

## Arguments

- g:

  An acyclic directed graph (DAG) structure of the Bayesian network for
  which we want to calculate the score.

- x:

  Either the data, or the population it would have come from. A
  `data.frame` or `matrix` of Gaussian data, with observations in the
  rows and random variables in the columns; or a population model built
  with
  [`population`](https://rcastelo.github.io/idlBNs/reference/population.md),
  in which case the score is evaluated in the large-sample limit instead
  of on a sample. That is not a different score: it is the same
  arithmetic fed the sufficient statistics a sample of the given size
  would have in expectation. A bare `GaussParDAG` from the pcalg package
  is accepted as a population model with a hard intervention to zero.
  See `target.index` for how the notional sample size is supplied.

- targets:

  (Default `list(integer(0))`) A `list` object with a family of targets
  provided as a list of integer vectors. Its default value indicates
  that there are no interventions in the data, i.e., the data is purely
  observational.

- target.index:

  (Default `NULL`) How much data comes from each environment. What it
  holds, and what `NULL` resolves to, depend on whether `x` carries data
  or a population.

  With data in `x`, a vector of integers in one-to-one correspondence
  with the rows in `x`, saying which target intervened on each row.
  `NULL` resolves to a vector of ones: the data is purely observational.

  With a population model in `x` there are no rows to label, so it is
  one observation count per element of `targets` instead. Their sum is
  the notional sample size, which sets the score's penalty term and so
  decides between nested models in the limit. `NULL` resolves to
  [`population`](https://rcastelo.github.io/idlBNs/reference/population.md)'s
  own `n`, split in proportion to \\(C, 1, \ldots, 1)\\ between the
  observational and the interventional environments. If `n` was not
  given there, omitting this argument is an error rather than a default:
  any size invented on the caller's behalf would silently change which
  model is selected.

  The counts need not be whole numbers, and a count of zero is allowed –
  that environment contributes nothing and is dropped. What each
  variable does need is at least two observations' worth of mass from
  the environments that leave it alone, the same floor the row counts
  must clear when `x` carries data; below it the variable cannot be
  scored and the call fails naming it.

  An environment with no observations is removed from `targets`
  altogether, under either kind of input – a count of zero, or a target
  that no row refers to. It is not merely that it cannot inform the
  score: the target family is also what defines *I*-equivalence and
  therefore which reversals are *I*-covered, so an intervention that was
  never performed would otherwise narrow the equivalence classes the
  search moves in and change the graph returned. The refinement is
  earned by having observed the intervention.

- cached.scores:

  (Default `NULL`) An optional list of environment objects, containing
  cached scores per parent set for each vertex in `g`. If `NULL`
  (default), no cached scores are used. Using this argument can speed up
  the calculation of the score when the same parent sets are scored
  multiple times. To use this argument, first create an empty
  environment object with
  `csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)`
  and then pass it to this `cached.scores` parameter, i.e.,
  `cached.scores=csco`.

- global.sufstats:

  (Default `NULL`) An optional list of global sufficient statistics for
  the iBIC score, as returned by the `.iBIC.global.sufstats()` function,
  which do not depend on the structure of a specific DAG, but only on
  the input data (`x`), the target vertices (`targets`) and the target
  indices (`target.index`) of the interventions. If `NULL` (default),
  the `.iBIC.global.sufstats()` function is internally called.

- engine:

  (Default `"C"`) A character string selecting the computation engine:
  `"C"` (default) uses a compiled C routine for speed; `"R"` uses the
  pure-R implementation and is provided for testing and verification.

- pasets:

  (Default `NULL`) An optional list of parent sets, one per vertex in
  `g` in the order given by `colnames(x)`, as internally built by
  `iBIC()` from the structure of `g`. If `NULL` (default), it is
  internally computed from `g`. Search algorithms that maintain `pasets`
  incrementally across many calls (e.g.
  [`hcmc()`](https://rcastelo.github.io/idlBNs/reference/hcmc.md),
  [`hillclimbing()`](https://rcastelo.github.io/idlBNs/reference/hillclimbing.md))
  can pass it in directly to skip rebuilding it from `g` on every call.

## Value

A single numeric value corresponding to the interventional BIC score of
the given structure of the Bayesian network for the given data set.

## References

Hauser, A. and Buehlmann, P. Jointly interventional and observational
data: estimation of interventional Markov equivalence classes of
directed acyclic graphs. *Journal of the Royal Statistical Society
Series B: Statistical Methodology*, 77:291-318, 2015.

Kalisch, M., Maechler, M., Colombo, D., Maathuis M.H. and Buehlmann, P.
Causal inference using graphical models with the R package pcalg.
*Journal of Statistical Software*, 47:1-26, 2012.

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
x <- rbind(obsdat, intdat)

## define the targets and target indices for the interventional data
targets <- list(integer(0), 2L)
target.index <- c(rep(1L, nobs), rep(2L, nint))

## calculate the interventional BIC score for the DAG and data set
iBIC(g, x, targets, target.index)
#> [1] -242.3883

## create another Markov equivalent DAG by reversing the arc X1 -> X2
## to obtain X1 <- X2 -> X3
g2 <- g
g2 <- removeEdge("X1", "X2", g2)
g2 <- addEdge("X2", "X1", g2)

## calculate the interventional BIC score for the new DAG on the
## same data, notice that the score is different despite being a
## Markov equivalent DAG
iBIC(g2, x, targets, target.index)
#> [1] -249.149

## this is not the case if we do not indicate the presence of interventions
## in the data
iBIC(g, x)
#> [1] -326.3833
iBIC(g2, x)
#> [1] -326.3833
```
