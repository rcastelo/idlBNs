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
  verbose = TRUE,
  engine = c("C", "R"),
  sampler = c("rcar", "exact"),
  escape = c("trials", "exhaustive"),
  escape.max = 512
)
```

## Arguments

- dat:

  A `data.frame` or `matrix` object, containing input Gaussian data,
  with data value records in the rows and random variables in the
  columns.

- r:

  (Default 20) Non-negative integer scalar indicating the maximum number
  of (*I*-)covered arc reversals by the RCAR algorithm (Castelo and
  Kočka, 2003).

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

- scorefun:

  (Default is
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md)) A
  function to calculate the goodness of fit (GoF) score of a DAG on a
  given data set.

- MAXTRIALS:

  (Default 5) Non-negative integer scalar indicating the maximum number
  of trials to escape from local maxima when `escape="trials"`. It is
  ignored when `escape="exhaustive"`, unless the exhaustive enumeration
  of the (*I*-)equivalence class exceeds `escape.max` or a chain
  component has more than 64 vertices, in which case the escape
  mechanism falls back to the RCAR algorithm for a maximum of
  `MAXTRIALS` trials. See the `escape` argument for details.

- verbose:

  (Default TRUE) Show progress in the calculations.

- engine:

  (Default `"C"`) A character string selecting the search engine: `"C"`
  (default) maintains the DAG, its ancestor relation and its parent sets
  in compiled code, and performs the (*I*-)covered arc reversals there
  too; `"R"` uses the pure-R implementation and is provided for testing
  and verification. Both consume the random number stream identically
  and so follow the same trajectory from the same seed, and both return
  the same result. `"C"` requires a `scorefun` able to score a whole
  neighbourhood at once, which
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md) and
  [`iBGe`](https://rcastelo.github.io/idlBNs/reference/iBGe.md) are;
  with any other score function the `"R"` engine is used regardless.

- sampler:

  (Default `"rcar"`) A character string selecting how the algorithm
  moves within the (*I*-)equivalence class of the current DAG. `"rcar"`
  performs the RCAR algorithm (a random walk of up to `r` (*I*-)covered
  arc reversals) of Castelo and Kočka (2003). `"exact"` instead draws a
  member of the class uniformly at random, by the Clique-Picking
  algorithm of Wienöbst *et al.* (2023), and ignores `r`. The walk
  produced by the RCAR algorithm and the exact draw are not equivalent:
  the walk is a random walk on the class, so its equilibrium is
  proportional to the number of (*I*-)covered arcs of each member and is
  not uniform for any value of `r`. The exact draw has no limit on the
  size of the class, and none on the number of vertices in the DAG. It
  declines, reverting to the walk and counting the draw in
  `sampler.fallbacks`, in two cases. First, vertex sets within an
  undirected chain component are 64-bit masks, so a component of more
  than 64 vertices is out of reach; components that large need a DAG
  with almost no immoralities, and the largest seen for random DAGs up
  to \\p = 500\\ is 16. Second, uniformity needs the per-component
  counts to be *exact*, because the draw compares a uniform integer
  against cumulative counts: doubles hold integers exactly only to
  \\2^{53}\\, and a clique of more than 18 vertices has a factorial past
  that, so beyond either bound the draw is declined rather than made
  with rounded weights. Counting is unaffected and still answers for a
  class of any size.

- escape:

  (Default `"trials"`) A character string selecting what happens at a
  local maximum. `"trials"` re-randomises the current DAG within its
  class and retries, up to `MAXTRIALS` times. `"exhaustive"` instead
  examines *every* member of the class and takes the best move available
  from any of them, which settles the question of whether the search is
  really at a local maximum of the class; `MAXTRIALS` is then unused. It
  falls back to `"trials"` in the cases where the enumeration declines,
  counted in `escape.fallbacks`: a class larger than `escape.max` or
  `INT_MAX`, and a class whose size cannot be computed at all because
  some chain component has more than 64 vertices. The second is the
  64-bit mask limit described under `sampler`, and applies here
  whichever `sampler` is in use, since the escape enumerates the class
  regardless of how the search moves within it.

- escape.max:

  (Default 512) Positive integer scalar giving the largest
  (*I*-)equivalence class the `escape="exhaustive"` enumeration will
  walk. That escape scores one whole neighbourhood per member, so its
  cost grows linearly in the size of the class, which is why it is the
  class size that is bounded. Producing the members is not itself the
  expensive part: each one is generated directly from its index in the
  class, so listing costs one step per member listed and no more. What
  `escape.max` bounds is therefore the scoring, not the enumeration.
  `NA` or `Inf` mean no limit, but the default of 512 is a safe value
  for the largest classes seen in random DAGs up to \\p = 500\\
  vertices.

## Value

A list with the following components: `dag`, a
[`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html) object
with the structure of the learned DAG; `sco`, its score;
`sampler.fallbacks`, the number of draws for which `sampler="exact"`
declined and the `"rcar"` walk was used instead; and `escape.fallbacks`,
the number of local maxima at which `escape="exhaustive"` declined and
the `MAXTRIALS` budget was used instead. Both counts are zero unless the
corresponding exact method was selected, and a non-zero count means part
of the run silently used the older machinery: the walk does not sample
the (*I*-)equivalence class uniformly, so a result obtained with
`sampler.fallbacks > 0` is not the one `sampler="exact"` promises.

## References

Castelo, R. and Kočka, T. On inclusion-driven learning of Bayesian
networks. *Journal of Machine Learning Research*, 4:527-574, 2003.

Castelo, R. Interventional idlBNs in DAG-space. In *Challenges and
Algorithms for Knowledge Discovery from Data*, M. van Leeuwen and J.
Vreeken (eds.). LNCS 16067, Festschrift, Springer, 2026.

Wienöbst, M., Bannach, M. and Liśkiewicz, M. Polynomial-time algorithms
for counting and sampling Markov equivalent DAGs with applications.
*Journal of Machine Learning Research*, 24(213):1-45, 2023.

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
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -87.004214222737 Escapes 0 Trials 0
#> ✔ Score -83.2933160642061 Escapes 0 Trials 0 [57ms]
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
#> $sampler.fallbacks
#> [1] 0
#> 
#> $escape.fallbacks
#> [1] 0
#> 

## calculate the structural Hamming distance (SHD) between the generative
## DAG and the estimated DAG
shd(e, dag2essgraph(dhat.hcmc$dag))
#> [1] 5

## run the iHCMC algorithm informing the presence of interventional data
## using by the default the interventional BIC score (see the iBIC()
## function).
dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex)
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -59.3908043567074 Escapes 0 Trials 0
#> ✔ Score -50.0166049861733 Escapes 0 Trials 0 [11ms]
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
#> $sampler.fallbacks
#> [1] 0
#> 
#> $escape.fallbacks
#> [1] 0
#> 

## the estimated DAG is closer to the generative DAG (lower SHD value)
## than the one estimated by HCMC, which did not take into account the
## presence of interventional data
shd(e, dag2essgraph(dhat.ihcmc$dag))
#> [1] 3

## run it again this time using the interventional BGe score (see the
## iBGe() function), which provides an estimate closer to the generative DAG
dhat.ihcmc2 <- hcmc(dat, targets=I, target.index=tindex, scorefun=iBGe)
#> ℹ Calculating global sufficient statistics
#> ⠙ Score -189.956899563879 Escapes 0 Trials 0
#> ✔ Score -181.640790427507 Escapes 0 Trials 0 [10ms]
#> 
shd(e, dag2essgraph(dhat.ihcmc2$dag))
#> [1] 2
```
