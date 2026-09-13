# Straightforward (classical) hill-climbing algorithm

Learn the structure of a Bayesian network from observational and
interventional data using a straightforward (classical) hill-climbing
algorithm that at each step during the search adds, removes and reverses
all possible arcs.

## Usage

``` r
hillclimbing(
  x,
  targets = list(integer(0)),
  target.index = NULL,
  scorefun = iBIC,
  verbose = TRUE,
  engine = c("C", "R")
)
```

## Arguments

- x:

  Either the data to learn from, or the population it would have come
  from. A `data.frame` or `matrix` of Gaussian data, with observations
  in the rows and random variables in the columns; or a population model
  built with
  [`population`](https://rcastelo.github.io/idlBNs/reference/population.md),
  in which case the search is scored in the large-sample limit instead
  of on a sample, which is what lets a run be read as the behaviour of
  the algorithm itself rather than of one dataset. A bare `GaussParDAG`
  from the pcalg package is accepted as a population model with a hard
  intervention to zero. See `target.index` for how the notional sample
  size is supplied.

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

- scorefun:

  (Default is
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md)) A
  function to calculate the goodness of fit (GoF) score of a DAG on a
  given data set.

- verbose:

  (Default TRUE) Show progress in the calculations.

- engine:

  (Default `"C"`) A character string selecting the search engine: `"C"`
  (default) maintains the DAG, its ancestor relation and its parent sets
  in compiled code; `"R"` uses the pure-R implementation and is provided
  for testing and verification. Both follow the same trajectory and
  return the same result. `"C"` requires a `scorefun` able to score a
  whole neighbourhood at once, which
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md) and
  [`iBGe`](https://rcastelo.github.io/idlBNs/reference/iBGe.md) are;
  with any other score function the `"R"` engine is used regardless.

## Value

A list containing a
[`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html) object
with the structure of the learned DAG, and its corresponding score.

## See also

[`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md),
[`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md)
