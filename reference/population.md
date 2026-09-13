# Create a population parameter model for interventional Bayesian networks

Wraps a generative model with the intervention semantics, so that it can
be passed where a data matrix normally goes and the scores are evaluated
in the large-sample limit instead of on a sample. See
[`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md) for what
the accompanying `targets` and `target.index` then mean.

## Usage

``` r
population(x, n = NULL, C = 1, ivent.value = 0, ivent.var = 0)
```

## Arguments

- x:

  A `GaussParDAG` object from the pcalg package: the generative model.

- n:

  (Default `NULL`) The notional sample size the population stands in
  for. It is not cosmetic: the score's penalty is a function of it, so
  it is what decides between nested models in the limit, and there is no
  defensible default. Setting it lets `target.index` be omitted, in
  which case `n` is divided between the environments in `targets`
  according to `C` – in proportion to \\(C, 1, \ldots, 1)\\, which is an
  equal split only when `C` is 1. Passing `target.index` states the
  allocation outright and uses neither.

- C:

  (Default 1) How much data an observational environment carries
  relative to each interventional one, as in Wang, Solus, Yang and Uhler
  (2017). With `k` interventional environments the notional size `n` is
  split in proportion to \\(C, 1, \ldots, 1)\\, so the observational
  stratum gets \\C/(C+k)\\ of it and each interventional stratum
  \\1/(C+k)\\. An environment counts as observational when its target is
  empty, so the split does not depend on the order of `targets`, and `C`
  is inert when every environment is of one kind. The default 1 is an
  equal split.

  This is the quantity their counterexample turns on: the graph they
  exhibit is a population local maximum for 2-14\\ for 68-96\\ passing
  counts states the allocation outright.

- ivent.value:

  The value an intervened variable is set to. It affects the score of
  the variables that are *not* intervened, because pooling environments
  whose means differ widens the covariance their scores are computed
  from.

- ivent.var:

  The variance around that value. Zero, the default, is a hard
  intervention to a constant; a positive value is a stochastic one.

## Value

An object of class `idlBNsPopulation`. It reports
[`ncol()`](https://rdrr.io/pkg/BiocGenerics/man/nrow.html) and
[`colnames()`](https://rdrr.io/pkg/BiocGenerics/man/row_colnames.html)
like the data matrix it stands in for, and
[`nrow()`](https://rdrr.io/pkg/BiocGenerics/man/nrow.html) as `NA`,
there being no rows to count.

## Examples

``` r

suppressPackageStartupMessages({
    library(graph)
    library(pcalg)
    library(idlBNs)
})

p <- 5
k <- 2

## simulate a random Gaussian DAG, as in the package vignette
set.seed(123)
Mg <- r.gauss.pardag(p, 0.6, top.sort=TRUE, normalize=TRUE)
g <- as(Mg, "graphNEL")

## a family of intervention targets: none, then two single vertices
I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))

## the equivalence class the search should recover, given those targets
e <- as(dag2essgraph(g, targets=I), "graphNEL")

## Rather than simulate data, describe the population it would come from:
## a notional 10000 observations, an observational environment five times
## the size of each interventional one, and an intervention that fixes its
## variable at 2. The scores are then evaluated in the large-sample limit,
## so a run shows what the algorithm does rather than what one sample did.
pop <- population(Mg, n=10000, C=5, ivent.value=2)

## it stands in for the data matrix, and reports that shape: p columns,
## and no rows to count
c(ncol(pop), nrow(pop))
#> [1]  5 NA

## run iHCMC in that limit, sampling the I-equivalence class exactly
set.seed(42)
dhat <- hcmc(pop, targets=I, sampler="exact", escape="exhaustive")
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -23594.4544223585 Escapes 0 Trials 0
#> ✔ Score -17250.4709628608 Escapes 0 Trials 0 [11ms]
#> 

## in the limit it recovers the generative equivalence class
shd(e, as(dag2essgraph(dhat$dag, targets=I), "graphNEL"))
#> [1] 0

## omitting target.index split n by C, in proportion to (C, 1, 1); saying
## so outright is the same experiment
w <- ifelse(lengths(I) == 0, 5, 1)
10000 * w / sum(w)
#> [1] 7142.857 1428.571 1428.571
set.seed(42)
same <- hcmc(pop, targets=I, target.index=10000*w/sum(w),
             sampler="exact", escape="exhaustive")
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -23594.4544223585 Escapes 0 Trials 0
#> ✔ Score -17250.4709628608 Escapes 0 Trials 0 [11ms]
#> 
identical(same$sco, dhat$sco)
#> [1] TRUE

## whereas C=1, an equal split, is a different experiment
set.seed(42)
eq <- hcmc(population(Mg, n=10000, ivent.value=2), targets=I,
           sampler="exact", escape="exhaustive")
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -21689.6925175966 Escapes 0 Trials 0
#> ✔ Score -16127.2069646389 Escapes 0 Trials 0 [15ms]
#> 
identical(eq$sco, dhat$sco)
#> [1] FALSE
```
