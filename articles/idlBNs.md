# The idlBNs package

Abstract

The `idlBNs` package implements inclusion-driven algorithms for learning
the structure of Bayesian networks. It currently provides the
hill-climber Monte Carlo (HCMC) algorithm for observational Gaussian
data, and the interventional HCMC (iHCMC) for interventional Gaussian
data.

## Introduction

Bayesian networks are multivariate models that assume that their joint
probability distribution factorizes according to an acyclic directed
graph (DAG). The DAG structure of a Bayesian network can be estimated
from data, and this is the so-called problem of *structure learning*.
The `idlBNs` package implements algorithms for learning the structure of
Bayesian networks from observational and interventional Gaussian data.

## Inclusion-driven structure learning

Structure learning algorithm for Bayesian networks can be classified
into three main categories: constraint-based, score-based, and hybrid
algorithms. The `idlBNs` package implements score-based algorithms that
attempt following an inclusion path in the search space of DAGs (Kočka
et al. 2001). Following an inclusion path confers makes the algorithm
consistent under the faithfulness assumption, that is, it converges to
the generative structure as the sample size grows large Castelo and
Kočka (2003).

Inclusion-driven learning of Bayesian networks, or *idlBNs* for
short[^1], can be implemented either in the space of the canonical
elements of the Markov equivalence classes of DAGs, also known as
*essential graphs*, as in the Greedy Equivalence Search (GES) algorithm
(Chickering 2002), or in the space of DAGs, as in the Hill-Climber Monte
Carlo (HCMC) algorithm (Kočka and Castelo 2001; Castelo and Kočka 2003).
The GES algorithm was adapted to interventional data in the Greedy
Interventional Equivalence Search (GIES) algorithm (Hauser and Bühlmann
2012), and an efficient implementation of both is available in the
[pcalg](https://cran.r-project.org/package=pcalg) package (Kalisch and
Bühlman 2007), through the functions
[`ges()`](https://rdrr.io/pkg/pcalg/man/ges.html) and
[`gies()`](https://rdrr.io/pkg/pcalg/man/gies.html), respectively.

The `idlBNs` provides an implementation in R of the HCMC algorithm for
observational Gaussian data, and of the interventional HCMC (iHCMC)
algorithm for a mixture of observational and interventional Gaussian
data (Castelo 2026) through the function
[`hcmc()`](https://rcastelo.github.io/idlBNs/reference/hcmc.md).

## Scoring functions for observational and interventional Gaussian data

The `idlBNs` package implements two different functions for scoring a
DAG structure from observational and interventional Gaussian data. The
[`iBIC()`](https://rcastelo.github.io/idlBNs/reference/iBIC.md)
function, which implements the interventional Bayesian Information
Criterion (iBIC) score (Hauser and Bühlmann 2012, 2015), and the
[`iBGe()`](https://rcastelo.github.io/idlBNs/reference/iBGe.md)
function, which implements the interventional Bayesian Gaussian
equivalent (iBGe) score (Kuipers and Moffa 2025).

## How to use the idlBNs package

Here we illustrate how to use the idlBNs package with multivariate
Gaussian data, simulated using the procedures of the
[pcalg](https://cran.r-project.org/package=pcalg) package for that
purpose. Besides this package, the idlBNs package currently depends also
on the Bioconductor package
[graph](https://bioconductor.org/packages/graph) for building and
representing DAGs through the object class `graphNEL`.

``` r

library(graph)
library(pcalg)
library(idlBNs)
```

We start by simulating a Gaussian Bayesian network with 10 vertices and
50% density using the function
[`r.gauss.pardag()`](https://rdrr.io/pkg/pcalg/man/r.gauss.pardag.html)
from the `pcalg` package, extracting its DAG strcture and converting it
into an essential graph using the function
[`dag2essgraph()`](https://rdrr.io/pkg/pcalg/man/dag2essgraph.html) from
the `pcalg` package.

``` r

p <- 5
d <- 0.6
set.seed(123)
M_G <- r.gauss.pardag(p, d, top.sort=TRUE, normalize=TRUE)
G <- as(M_G, "graphNEL")
E <- dag2essgraph(as(G, "graphNEL"))
```

Assume we want to consider two different singleton intervention targets
($`k=2`$), i.e., we select two vertices from the DAG structure uniformly
at random to intervene in their associated random variables, but we
intervene only in one random variable at a time in a given multivariate
observation. Following the terminology of Hauser and Bühlmann (2012), we
define a family of *intervention targets*
$`\mathcal{I} = \{I_1, I_2, I_3\}`$, where $`I_1=\{\emptyset\}`$ and
$`I_i\{v_i\}`$ with $`i=\{2, 3\}`$ and $`v_i\in V`$ is a randomly
selected vertex from the DAG structure $`G=(V, E)`$.

``` r
k <- 2
I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))
I
[[1]]
integer(0)

[[2]]
[1] 5

[[3]]
[1] 3
```

Consider simulating $`n=100`$ multivariate observations, with a similar
number of them for each intervention target.

``` r
n <- 30
nbytgts <- rep(floor(n / (k + 1)), k)
nbytgts <- c(n - sum(nbytgts), nbytgts)
nbytgts
[1] 10 10 10
```

Simulate mixed observational and interventional multivariate Gaussian
data from the DAG structure $`G`$ and the family of intervention targets
$`\mathcal{I}`$ using the function
[`rmvnorm.ivent()`](https://rdrr.io/pkg/pcalg/man/rmvnorm.ivent.html)
from the `pcalg` package, storing the index to the corresponding
intervention target for each row of the data.

``` r

dat <- list()
for (v in seq_along(I)) {
  targets <- I[[v]]
  dat[[v]] <- rmvnorm.ivent(nbytgts[v], M_G, target=targets,
                            target.value=rep(2, length(targets)))
}
dat <- do.call("rbind", dat)
tindex <- rep(1:length(nbytgts), nbytgts)
```

Run the iHCMC algorithm informing the presence of interventional data.
By default, the iBIC score is used in the search algorithm.

``` r
Dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex)
Dhat.ihcmc
$dag
A graphNEL graph with directed edges
Number of Nodes = 5 
Number of Edges = 4 

$sco
[1] -50.0166
```

Calculate the structural Hamming distance (SHD) to the essential graph
of the generative DAG structure $`G`$.

``` r
shd(dag2essgraph(Dhat.ihcmc$dag), E)
[1] 3
```

Repeat the search using the iBGe score instead of the iBIC score.

``` r
Dhat.ihcmc.ibge <- hcmc(dat, targets=I, target.index=tindex, scorefun=iBGe)
shd(dag2essgraph(Dhat.ihcmc.ibge$dag), E)
[1] 2
```

We can see that the iBGe score leads to an estimate of the DAG structure
that is close to the generative DAG structure, in terms of the SHD. Now
use the HCMC algorithm without informing the presence of interventional
data, i.e., assuming that all data are observational.

``` r
Dhat.hcmc <- hcmc(dat)
shd(dag2essgraph(Dhat.hcmc$dag), E)
[1] 5
Dhat.hcmc.ibge <- hcmc(dat, scorefun=iBGe)
shd(dag2essgraph(Dhat.hcmc.ibge$dag), E)
[1] 4
```

In this setting, using either the BIC or the BGe score leads to worse
estimates of the DAG structure.

## Session information

``` r
sessionInfo()
R version 4.6.0 (2026-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] idlBNs_1.0.0        pcalg_2.7-12        graph_1.90.0       
[4] BiocGenerics_0.58.1 generics_0.1.4      knitr_1.51         

loaded via a namespace (and not attached):
 [1] jsonlite_2.0.0      compiler_4.6.0      BiocManager_1.30.27
 [4] Rcpp_1.1.1-1.1      fastICA_1.2-7       cluster_2.1.8.2    
 [7] jquerylib_0.1.4     systemfonts_1.3.2   textshaping_1.0.5  
[10] yaml_2.3.12         fastmap_1.2.0       R6_2.6.1           
[13] igraph_2.3.2        RBGL_1.88.0         robustbase_0.99-7  
[16] bdsmatrix_1.3-7     desc_1.4.3          sfsmisc_1.1-24     
[19] bslib_0.11.0        rlang_1.2.0         cachem_1.1.0       
[22] xfun_0.58           fs_2.1.0            sass_0.4.10        
[25] cli_3.6.6           pkgdown_2.2.0       magrittr_2.0.5     
[28] digest_0.6.39       lifecycle_1.0.5     clue_0.3-68        
[31] DEoptimR_1.1-4      evaluate_1.0.5      ggm_2.5.2          
[34] corpcor_1.6.10      ragg_1.5.2          abind_1.4-8        
[37] stats4_4.6.0        rmarkdown_2.31      tools_4.6.0        
[40] pkgconfig_2.0.3     htmltools_0.5.9    
```

## References

Castelo, Robert. 2026. “Interventional idlBNs in DAG-Space.” In
*Challenges and Algorithms for Knowledge Discovery from Data*, edited by
J. Vreeken M. van Leeuwen. Springer.

Castelo, Robert, and Tomás Kočka. 2003. “On Inclusion-Driven Learning of
Bayesian Networks.” *Journal of Machine Learning Research* 4 (Sep):
527–74.

Chickering, David Maxwell. 2002. “Optimal Structure Identification with
Greedy Search.” *Journal of Machine Learning Research* 3 (Nov): 507–54.

Hauser, Alain, and Peter Bühlmann. 2012. “Characterization and Greedy
Learning of Interventional Markov Equivalence Classes of Directed
Acyclic Graphs.” *Journal of Machine Learning Research* 13 (1): 2409–64.

Hauser, Alain, and Peter Bühlmann. 2015. “Jointly Interventional and
Observational Data: Estimation of Interventional Markov Equivalence
Classes of Directed Acyclic Graphs.” *Journal of the Royal Statistical
Society Series B: Statistical Methodology* 77 (1): 291–318.

Kalisch, Markus, and Peter Bühlman. 2007. “Estimating High-Dimensional
Directed Acyclic Graphs with the PC-Algorithm.” *Journal of Machine
Learning Research* 8 (3).

Kočka, T., and R. Castelo. 2001. “Improved Learning of Bayesian
Networks.” In *Proc. Of the Conf. On Uncertainty in Artificial
Intelligence*, edited by J. Breese and D. Koller. Morgan Kaufmann.

Kočka, Tomás, Remco Bouckaert, and Milan Studený. 2001. “On
Characterizing Inclusion of Bayesian Networks.” In *Proc. Of the
Conf. On Uncertainty in Artificial Intelligence*, edited by J. Breese
and D. Koller. Morgan Kaufmann.

Kuipers, Jack, and Giusi Moffa. 2025. “The Interventional Bayesian
Gaussian Equivalent Score for Bayesian Causal Inference with Unknown
Soft Interventions.” *Causal Learning and Reasoning*, 772–91.

[^1]: Pronounced *ideal BNs*
