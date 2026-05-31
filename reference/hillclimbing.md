# Straightforward (classical) hill-climbing algorithm

Learn the structure of a Bayesian network from observational data using
a straightforward (classical) hill-climbing algorithm that at each step
during the search adds, removes and reverses all possible arcs.

## Usage

``` r
hillclimbing(dat, scorefun = iBIC, verbose = TRUE)
```

## Arguments

- dat:

  A `data.frame` object with data records in the rows.

- scorefun:

  (Default is
  [`iBIC`](https://rcastelo.github.io/idlBNs/reference/iBIC.md)) A
  function to calculate the goodness of fit (GoF) score of a DAG on a
  given data set.

- verbose:

  (Default TRUE) Show progress in the calculations.

## Value

A list containing a
[`graphNEL`](https://rdrr.io/pkg/graph/man/graphNEL-class.html) object
with the structure of the learned DAG, and its corresponding score.

## See also

[`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md),
[`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md)
