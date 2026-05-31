#' idlBNs: Inclusion-driven learning of Bayesian networks
#' 
#' The idlBNs package implements inclusion-driven algorithms for learning the
#' structure of Bayesian networks. It currently provides the hill-climber Monte
#' Carlo (HCMC) algorithm for observational Gaussian data, and the
#' interventional HCMC (iHCMC) for interventional Gaussian data.
#' 
#' The main functions are:
#'
#' \itemize{
#'     \item \code{\link{iBIC}()} - the interventional BIC score for Gaussian data.
#'     \item \code{\link{iBGe}()} - the interventional BGe score for Gaussian data.
#'     \item \code{\link{hillclimbing}()} - a straightforward hill-climbing algorithm for learning the structure of Bayesian networks.
#'     \item \code{\link{hcmc}()} - the HCMC and iHCMC algorithms for learning the structure of Bayesian networks from observational and interventional data, respectively.
#' }
#' 
#' For detailed information on usage, see the package vignette, by typing
#' \code{vignette("idlBNs")}.
#' 
#' The code of the development version of the package is available at the
#' GitHub repository \url{https://github.com/rcastelo/idlBNs}.
#' 
#' Questions and bug reports should be posted by opening an issue in the
#' idlBNs GitHub repo at \url{https://github.com/rcastelo/idlBNs/issues}.
#'
#' @name idlBNs-package
#' @aliases idlBNs-package
#' @aliases idlBNs
#' @keywords package
"_PACKAGE"

NULL
