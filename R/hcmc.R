#' The HCMC and iHCMC algorithms
#'
#' @description Run the hill-climber Monte Carlo (HCMC) algorithm (Castelo and
#' Kočka, 2003) on purely observational Gaussian data, or the interventional
#' HCMC (iHCMC) on mixed observational and interventional Gaussian data
#' (Castelo, 2026).
#'
#' @param dat A `data.frame` or `matrix` object, containing input Gaussian data,
#' with data value records in the rows and random variables in the columns.
#'
#' @param r (Default 20) Maximum number of (\emph{I}-)covered arc reversals.
#'
#' @param targets (Default a list with an empty integer vector) Family of
#' intervention targets provided as a list of integer vectors. The default value
#' implies that there are no interventions and the data is purely observational.
#'
#' @param target.index (Default an empty integer vector) A vector of integers
#' in one-to-one correspondence with the rows in `dat`, indicating which rows
#' in the input data are intervened by which targets.
#'
#' @param MAXTRIALS (Default 5) Maximum number of trials to escape from local
#' maxima.
#'
#' @param scorefun (Default is [`iBIC`]) A function to calculate the goodness
#' of fit (GoF) score of a DAG on a given data set.
#'
#' @param verbose (Default TRUE) Show progress in the calculations.
#'
#' @return A list containing a [`graphNEL`][graph::graphNEL-class] object with
#' the structure of the learned DAG, and its corresponding score.
#' 
#' @references Castelo, R. and Kočka, T. On inclusion-driven learning of
#' Bayesian networks. *Journal of Machine Learning Research*, 4:527-574, 2003.
#'
#' @references Castelo, R. Interventional idlBNs in DAG-space. In *Challenges
#' and Algorithms for Knowledge Discovery from Data*, M. van Leeuwen and
#' J.  Vreeken (eds.). LNCS 16067, Festschrift, Springer, 2026.
#'
#' @seealso [iBIC()], [iBGe()]
#'
#' @examples
#'
#' suppressPackageStartupMessages({
#'   library(graph)
#'   library(pcalg)
#'   library(idlBNs)
#' })
#' 
#' p <- 5
#' k <- 2
#' n <- 30
#'
#' ## simulate a random DAG
#' set.seed(123)
#' Mg <- r.gauss.pardag(p, 0.6, top.sort=TRUE, normalize=TRUE)
#' g <- as(Mg, "graphNEL")
#' e <- dag2essgraph(as(g, "graphNEL"))
#' 
#' ## generate a random family of intervention targets
#' I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))
#'
#' ## sample size per different target (including the no-target)
#' nbytgts <- rep(floor(n / (k + 1)), k)
#' nbytgts <- c(n - sum(nbytgts), nbytgts)
#'
#' ## simulate mixed observational and interventional data
#' dat <- list()
#' for (v in seq_along(I)) {
#'   targets <- I[[v]]
#'   dat[[v]] <- rmvnorm.ivent(nbytgts[v], Mg, target=targets,
#'                             target.value=rep(2, length(targets)))
#' }
#' dat <- do.call("rbind", dat)
#'
#' ## store the target index for each row of the data
#' tindex <- rep(1:length(nbytgts), nbytgts)
#'
#' ## run the HCMC algorithm assuming all data were observational
#' dhat.hcmc <- hcmc(dat)
#' dhat.hcmc
#'
#' ## calculate the structural Hamming distance (SHD) between the generative
#' ## DAG and the estimated DAG
#' shd(e, dag2essgraph(dhat.hcmc$dag))
#'
#' ## run the iHCMC algorithm informing the presence of interventional data
#' ## using by the default the interventional BIC score (see the iBIC()
#' ## function).
#' dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex)
#' dhat.ihcmc
#'
#' ## the estimated DAG is closer to the generative DAG (lower SHD value)
#' ## than the one estimated by HCMC, which did not take into account the
#' ## presence of interventional data
#' shd(e, dag2essgraph(dhat.ihcmc$dag))
#'
#' ## run it again this time using the interventional BGe score (see the
#' ## iBGe() function), which provides an estimate closer to the generative DAG
#' dhat.ihcmc2 <- hcmc(dat, targets=I, target.index=tindex, scorefun=iBGe)
#' shd(e, dag2essgraph(dhat.ihcmc2$dag))
#'
#' @importFrom graph nodes edgeMatrix graphNEL
#' @importClassesFrom graph graphNEL
#' @importFrom cli cli_alert_info cli_progress_bar cli_progress_done
#' @importFrom cli cli_alert_success cli_alert_warning
#' @export
#' @rdname hcmc

hcmc <- function(dat, r=20, targets=list(integer(0)),
                 target.index=rep(1L, nrow(dat)),
                 scorefun=iBIC, MAXTRIALS=5, verbose=TRUE) {

  stopifnot(is.list(targets)) ## QC
  scorefun <- match.fun(scorefun)

  utargets <- sort(unique(unlist(targets)))
  dag <- graphNEL(colnames(dat), edgemode="directed")
  s0 <- -Inf
  s1 <- scorefun(dag, dat, targets, target.index)
  was_in_local_maximum <- local_maximum <- s1 < s0
  trials <- escapes <- avg_trials_per_escape <- 0

  if (verbose) {
    msg <- "Score {s1} Escapes {escapes} Trials {avg_trials_per_escape}"
    cli_progress_step(msg, spinner=TRUE)
  }

  while (!local_maximum) {
    s0 <- s1
    dag <- rcar(dag, r, utargets)
    ne <- ncr.nh(dag, utargets)
    s1 <- sapply(ne, function(g, d, tgts, tgts.idx)
                       scorefun(g, d, tgts, tgts.idx),
                 dat, targets, target.index)
    dag1 <- ne[[which.max(s1)]]
    s1 <- max(s1)
    local_maximum <- s1 <= s0
    if (!local_maximum) {
      dag <- dag1
      if (was_in_local_maximum) {
        escapes <- escapes + 1
        avg_trials_per_escape <- (avg_trials_per_escape *
                                  (escapes-1) + trials) / escapes
        was_in_local_maximum <- FALSE
      }
      trials <- 0
    } else if (trials < MAXTRIALS) {
      s1 <- s0
      dag <- rcar(dag, r, utargets)
      local_maximum <- FALSE
      was_in_local_maximum <- TRUE
      trials <- trials + 1
    } else
      s1 <- s0

    if (verbose)
      cli_progress_update()
  }

  list(dag=dag, sco=s1)
}
