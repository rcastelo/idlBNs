# The HCMC and iHCMC algorithms

Run the hill-climber Monte Carlo (HCMC) algorithm (Castelo and Kočka,
2003) on purely observational Gaussian data, or the interventional HCMC
(iHCMC) on mixed observational and interventional Gaussian data
(Castelo, 2026).

## Usage

``` r
hcmc(
  dat,
  r = 20,
  targets = list(integer(0)),
  target.index = rep(1L, nrow(dat)),
  scorefun = iBIC,
  MAXTRIALS = 5,
  verbose = TRUE
)
```

## Arguments

- dat:

  A `data.frame` or `matrix` object, containing input Gaussian data,
  with data value records in the rows and random variables in the
  columns.

- r:

  (Default 20) Maximum number of (*I*-)covered arc reversals.

- targets:

  (Default a list with an empty integer vector) Family of intervention
  targets provided as a list of integer vectors. The default value
  implies that there are no interventions and the data is purely
  observational.

- target.index:

  (Default an empty integer vector) A vector of integers in one-to-one
  correspondence with the rows in `dat`, indicating which rows in the
  input data are intervened by which targets.

- scorefun:

  (Default is
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md)) A
  function to calculate the goodness of fit (GoF) score of a DAG on a
  given data set.

- MAXTRIALS:

  (Default 5) Maximum number of trials to escape from local maxima.

- verbose:

  (Default TRUE) Show progress in the calculations.

## Value

A list containing a
[`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html) object
with the structure of the learned DAG, and its corresponding score.

## References

Castelo, R. and Kočka, T. On inclusion-driven learning of Bayesian
networks. *Journal of Machine Learning Research*, 4:527-574, 2003.

Castelo, R. Interventional idlBNs in DAG-space. In *Challenges and
Algorithms for Knowledge Discovery from Data*, M. van Leeuwen and J.
Vreeken (eds.). LNCS 16067, Festschrift, Springer, 2026.

## See also

[`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md),
[`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md)

## Examples

``` r

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

p <- 5
k <- 2
n <- 30

## simulate a random DAG
set.seed(123)
Mg <- r.gauss.pardag(p, 0.6, top.sort=TRUE, normalize=TRUE)
g <- as(Mg, "graphNEL")
e <- dag2essgraph(as(g, "graphNEL"))

## generate a random family of intervention targets
I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))

## sample size per different target (including the no-target)
nbytgts <- rep(floor(n / (k + 1)), k)
nbytgts <- c(n - sum(nbytgts), nbytgts)

## simulate mixed observational and interventional data
dat <- list()
for (v in seq_along(I)) {
  targets <- I[[v]]
  dat[[v]] <- rmvnorm.ivent(nbytgts[v], Mg, target=targets,
                            target.value=rep(2, length(targets)))
}
dat <- do.call("rbind", dat)

## store the target index for each row of the data
tindex <- rep(1:length(nbytgts), nbytgts)

## run the HCMC algorithm assuming all data were observational
dhat.hcmc <- hcmc(dat)
#> ⠙ Score -87.004214222737 Escapes 0 Trials 0
#> ⠹ Score -83.2933160642061 Escapes 0 Trials 0
#> ✔ Score -83.2933160642061 Escapes 0 Trials 0 [454ms]
#> 
dhat.hcmc
#> $dag
#> A graphNEL graph with directed edges
#> Number of Nodes = 5 
#> Number of Edges = 2 
#> 
#> $sco
#> [1] -83.29332
#> 

## calculate the structural Hamming distance (SHD) between the generative
## DAG and the estimated DAG
shd(e, dag2essgraph(dhat.hcmc$dag))
#> [1] 5

## run the iHCMC algorithm informing the presence of interventional data
## using by the default the interventional BIC score (see the iBIC()
## function).
dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex)
#> ⠙ Score -59.3908043567074 Escapes 0 Trials 0
#> ✔ Score -50.0166049861733 Escapes 0 Trials 0 [455ms]
#> 
dhat.ihcmc
#> $dag
#> A graphNEL graph with directed edges
#> Number of Nodes = 5 
#> Number of Edges = 4 
#> 
#> $sco
#> [1] -50.0166
#> 

## the estimated DAG is closer to the generative DAG (lower SHD value)
## than the one estimated by HCMC, which did not take into account the
## presence of interventional data
shd(e, dag2essgraph(dhat.ihcmc$dag))
#> [1] 3

## run it again this time using the interventional BGe score (see the
## iBGe() function), which provides an estimate closer to the generative DAG
dhat.ihcmc2 <- hcmc(dat, targets=I, target.index=tindex, scorefun=iBGe)
#> ⠙ Score -189.956899563879 Escapes 0 Trials 0
#> ✔ Score -181.640790427507 Escapes 0 Trials 0 [508ms]
#> 
shd(e, dag2essgraph(dhat.ihcmc2$dag))
#> [1] 2
```
