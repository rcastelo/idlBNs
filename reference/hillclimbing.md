# Straightforward (classical) hill-climbing algorithm

Learn the structure of a Bayesian network from observational and
interventional data using a straightforward (classical) hill-climbing
algorithm that at each step during the search adds, removes and reverses
all possible arcs.

## Usage

``` r
hillclimbing(
  dat,
  targets = list(integer(0)),
  target.index = rep(1L, nrow(dat)),
  scorefun = iBIC,
  verbose = TRUE,
  engine = c("C", "R")
)
```

## Arguments

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
