# idlBNs: Inclusion-driven learning of Bayesian networks

The idlBNs package implements inclusion-driven algorithms for learning
the structure of Bayesian networks. It currently provides the
hill-climber Monte Carlo (HCMC) algorithm for observational Gaussian
data, and the interventional HCMC (iHCMC) for interventional Gaussian
data.

## Details

The main functions are:

- [`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md) - the
  interventional BIC score for Gaussian data.

- [`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md) - the
  interventional BGe score for Gaussian data.

- [`hillclimbing()`](https://rcastelo.github.io/idlBNs/reference/hillclimbing.md) -
  a straightforward hill-climbing algorithm for learning the structure
  of Bayesian networks.

- [`hcmc()`](https://rcastelo.github.io/idlBNs/reference/hcmc.md) - the
  HCMC and iHCMC algorithms for learning the structure of Bayesian
  networks from observational and interventional data, respectively.

For detailed information on usage, see the package vignette, by typing
[`vignette("idlBNs")`](https://rcastelo.github.io/idlBNs/articles/idlBNs.md).

The code of the development version of the package is available at the
GitHub repository <https://github.com/rcastelo/idlBNs>.

Questions and bug reports should be posted by opening an issue in the
idlBNs GitHub repo at <https://github.com/rcastelo/idlBNs/issues>.

## See also

Useful links:

- <https://rcastelo.github.io/idlBNs>

- Report bugs at <https://github.com/rcastelo/idlBNs/issues>

## Author

**Maintainer**: Robert Castelo <robert.castelo@upf.edu>

Authors:

- Robert Castelo <robert.castelo@upf.edu>
