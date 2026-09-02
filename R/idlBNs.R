#' idlBNs: Inclusion-driven learning of Bayesian networks
#' 
#' The idlBNs package implements inclusion-driven algorithms for learning the
#' structure of Bayesian networks. It currently provides the hill-climber Monte
#' Carlo (HCMC) algorithm for observational Gaussian data, and the
#' interventional HCMC (iHCMC) for mixture of observational and interventional
#' Gaussian data.
#' 
#' The main functions are:
#'
#' \itemize{
#'     \item \code{\link{iBIC}()} - the interventional BIC score for Gaussian
#'           data (Hauser and Bühlmann, 2015).
#'     \item \code{\link{iBGe}()} - the interventional BGe score for Gaussian
#'           data (Kuipers and Moffa, 2025).
#'     \item \code{\link{hillclimbing}()} - a straightforward hill-climbing
#'           algorithm for learning the structure of Bayesian networks.
#'     \item \code{\link{hcmc}()} - the HCMC (Castelo and Kočka, 2003) and
#'           iHCMC algorithms (Castelo, 2026) for learning the structure of
#'           Bayesian networks from observational and interventional data,
#'           respectively.
#' }
#' 
#' For detailed information on usage, see the package vignette at
#' \url{https://rcastelo.github.io/idlBNs/articles/idlBNs.html} or by typing
#' \code{vignette("idlBNs")} at the R console.
#' 
#' If you use the HCMC algorithm in you research, please cite
#' (Castelo and Kočka, 2003). If you use the iHCMC algorithm, please cite
#' (Castelo, 2026). If you use the iBIC score, please cite
#' (Hauser and Bühlmann, 2015). If you use the iBGe score, please cite
#' (Kuipers and Moffa, 2025).
#'
#' The code of the development version of the package is available at the
#' GitHub repository \url{https://github.com/rcastelo/idlBNs}.
#' 
#' Questions and bug reports should be posted by opening an issue in the
#' idlBNs GitHub repo at \url{https://github.com/rcastelo/idlBNs/issues}.
#'
#' @author Robert Castelo
#' @name idlBNs-package
#' @aliases idlBNs-package
#' @aliases idlBNs
#'
#' @references Castelo, R. and Kočka T. On inclusion-driven learning of Bayesian networks. \emph{Journal of Machine Learning Research}, 4:527-574, 2003.
#'
#' @references Castelo, R. Interventional idlBNs in DAG-Space. In \emph{Challenges and Algorithms for Knowledge Discovery from Data}, M. van Leeuwen and J. Vreeken, Eds., pp. 161-177, Springer Nature, ISBN 978-3-032-03028-3, 2026.
#'
#' @references Hauser, A. and Bühlmann, P. Jointly interventional and observational data: estimation of interventional markov equivalence classes of directed acyclic graphs. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 77(1):291-318, 2015.
#'
#' @references Kuipers, J. and Moffa, G. The interventional {B}ayesian {G}aussian equivalent score for {B}ayesian causal inference with unknown soft interventions. In \emph{Causal Learning and Reasoning}, pages 772--791. PMLR, 2025.
#'
#' @keywords package
"_PACKAGE"

NULL
