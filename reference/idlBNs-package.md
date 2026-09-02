# idlBNs: Inclusion-driven learning of Bayesian networks

The idlBNs package implements inclusion-driven algorithms for learning
the structure of Bayesian networks. It currently provides the
hill-climber Monte Carlo (HCMC) algorithm for observational Gaussian
data, and the interventional HCMC (iHCMC) for mixture of observational
and interventional Gaussian data.

## Details

The main functions are:

- [`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md) - the
  interventional BIC score for Gaussian data (Hauser and Bühlmann,
  2015).

- [`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md) - the
  interventional BGe score for Gaussian data (Kuipers and Moffa, 2025).

- [`hillclimbing()`](https://rcastelo.github.io/idlBNs/reference/hillclimbing.md) -
  a straightforward hill-climbing algorithm for learning the structure
  of Bayesian networks.

- [`hcmc()`](https://rcastelo.github.io/idlBNs/reference/hcmc.md) - the
  HCMC (Castelo and Kočka, 2003) and iHCMC algorithms (Castelo, 2026)
  for learning the structure of Bayesian networks from observational and
  interventional data, respectively.

For detailed information on usage, see the package vignette at
<https://rcastelo.github.io/idlBNs/articles/idlBNs.html> or by typing
[`vignette("idlBNs")`](https://rcastelo.github.io/idlBNs/articles/idlBNs.md)
at the R console.

If you use the HCMC algorithm in you research, please cite (Castelo and
Kočka, 2003). If you use the iHCMC algorithm, please cite (Castelo,
2026). If you use the iBIC score, please cite (Hauser and Bühlmann,
2015). If you use the iBGe score, please cite (Kuipers and Moffa, 2025).

The code of the development version of the package is available at the
GitHub repository <https://github.com/rcastelo/idlBNs>.

Questions and bug reports should be posted by opening an issue in the
idlBNs GitHub repo at <https://github.com/rcastelo/idlBNs/issues>.

## References

Castelo, R. and Kočka T. On inclusion-driven learning of Bayesian
networks. *Journal of Machine Learning Research*, 4:527-574, 2003.

Castelo, R. Interventional idlBNs in DAG-Space. In *Challenges and
Algorithms for Knowledge Discovery from Data*, M. van Leeuwen and J.
Vreeken, Eds., pp. 161-177, Springer Nature, ISBN 978-3-032-03028-3,
2026.

Hauser, A. and Bühlmann, P. Jointly interventional and observational
data: estimation of interventional markov equivalence classes of
directed acyclic graphs. *Journal of the Royal Statistical Society
Series B: Statistical Methodology*, 77(1):291-318, 2015.

Kuipers, J. and Moffa, G. The interventional Bayesian Gaussian
equivalent score for Bayesian causal inference with unknown soft
interventions. In *Causal Learning and Reasoning*, pages 772–791. PMLR,
2025.

## See also

Useful links:

- <https://rcastelo.github.io/idlBNs/>

- Report bugs at <https://github.com/rcastelo/idlBNs/issues>

## Author

Robert Castelo
