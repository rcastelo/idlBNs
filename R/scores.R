#' @title BIC score for observational and interventional Gaussian data
#'
#' @description Score the goodness-of-fit (GoF) of a given structure of a
#' Bayesian network given an interventional data set of continuous values, where
#' observations are assumed to be independent but not identically distributed
#' (not iid) multivariate Gaussian. This GoF score corresponds to the Bayesian
#' information criterion (BIC) as implemented in the `GaussL0penIntScore` class
#' from the `pcalg` package (Kalisch et al., 2012). By default, the arguments
#' `targets` and `target.index` are set so that the calculated BIC score
#' assumes there are no interventions in the data.
#' 
#' @param g An acyclic directed graph (DAG) structure of the Bayesian network
#' for which we want to calculate the score.
#'
#' @param dat A `data.frame` object with data records in the rows.
#'
#' @param targets (Default `list(integer(0))`) A `list` object with a family of
#' targets provided as a list of integer vectors. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param target.index (Default a unit vector) A vector of integers in
#' one-to-one correspondence with the rows in `dat`, indicating which rows in
#' the input data are intervened by which targets. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param cached.scores (Default `NULL`) An optional list of environment
#' objects, containing cached scores per parent set for each vertex in `g`. If
#' `NULL` (default), no cached scores are used. Using this argument can speed
#' up the calculation of the score when the same parent sets are scored multiple
#' times. To use this argument, first create an empty environment object with
#' `csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)`
#' and then pass it to this `cached.scores` parameter, i.e.,
#' `cached.scores=csco`.
#'
#' @param global.sufstats (Default `NULL`) An optional list of global sufficient
#' statistics for the iBIC score, as returned by the `.iBIC.global.sufstats()`
#' function, which do not depend on the structure of a specific DAG, but only
#' on the input data (`dat`), the target vertices (`targets`) and the target
#' indices (`target.index`) of the interventions. If `NULL` (default), the
#' `.iBIC.global.sufstats()` function is internally called.
#'
#' @return A single numeric value corresponding to the interventional BIC score
#' of the given structure of the Bayesian network for the given data set.
#'
#' @references Hauser, A. and Buehlmann, P. Jointly interventional and
#' observational data: estimation of interventional Markov equivalence classes
#' of directed acyclic graphs. *Journal of the Royal Statistical Society Series
#' B: Statistical Methodology*, 77:291-318, 2015.
#'
#' @references Kalisch, M., Maechler, M., Colombo, D., Maathuis M.H. and
#' Buehlmann, P. Causal inference using graphical models with the R package
#' pcalg. *Journal of Statistical Software*, 47:1-26, 2012.
#'
#' @examples
#'
#' library(graph)
#'
#' p <- 3
#' nobs <- 100
#' nint <- 100
#' n <- nobs + nint
#'
#' ## define a DAG structure of a Bayesian network with three vertices
#' ## forming a Markov chain X1 -> X2 -> X3
#' g <- new("graphNEL", nodes=c("X1", "X2", "X3"), edgemode="directed")
#' g <- addEdge("X1", "X2", g)
#' g <- addEdge("X2", "X3", g)
#'
#' ## simulate observational data for the previous DAG X1 -> X2 -> X3
#' set.seed(123)
#' X1 <- rnorm(nobs, mean=0, sd=1)
#' X2 <- 0.5 * X1 + rnorm(nobs, mean=0, sd=1)
#' X3 <- 0.5 * X2 + rnorm(nobs, mean=0, sd=1)
#' obsdat <- data.frame(X1=X1, X2=X2, X3=X3)
#'
#' ## simulate interventional data for the same DAG, where X2 is intervened
#' X1 <- rnorm(nint, mean=0, sd=1)
#' X2 <- rnorm(nint, mean=0, sd=1) + 1.0
#' X3 <- 0.5 * X2 + rnorm(nint, mean=0, sd=1)
#' intdat <- data.frame(X1=X1, X2=X2, X3=X3)
#'
#' ## combine observational and interventional data
#' dat <- rbind(obsdat, intdat)
#'
#' ## define the targets and target indices for the interventional data
#' targets <- list(integer(0), 2L)
#' target.index <- c(rep(1L, nobs), rep(2L, nint))
#'
#' ## calculate the interventional BIC score for the DAG and data set
#' iBIC(g, dat, targets, target.index)
#'
#' ## create another Markov equivalent DAG by reversing the arc X1 -> X2
#' ## to obtain X1 <- X2 -> X3
#' g2 <- g
#' g2 <- removeEdge("X1", "X2", g2)
#' g2 <- addEdge("X2", "X1", g2)
#'
#' ## calculate the interventional BIC score for the new DAG on the
#' ## same data, notice that the score is different despite being a
#' ## Markov equivalent DAG
#' iBIC(g2, dat, targets, target.index)
#'
#' ## this is not the case if we do not indicate the presence of interventions
#' ## in the data
#' iBIC(g, dat)
#' iBIC(g2, dat)
#'
#' @importFrom graph numNodes edgeMatrix
#' @export
iBIC <- function(g, dat, targets=list(integer(0)),
                 target.index=rep(1L, nrow(dat)),
                 cached.scores=NULL, global.sufstats=NULL) {

  if (is.null(attr(dat, "sanitycheck"))) {
    dat <- .check_input_data(dat)
    .check_g_dat_consistency(g, dat)
  }

  v <- match(nodes(g), colnames(dat))
  p <- numNodes(g)
  n <- nrow(dat)
  em <- edgeMatrix(g)
  pasets <- split(em["from", ], factor(v[em["to", ]], levels=v))
  stopifnot(identical(names(pasets), as.character(v)))
  .check_cached_scores(g, cached.scores)

  if (is.null(global.sufstats))
    global.sufstats <- .iBIC.global.sufstats(dat, targets, target.index)
  onlyobsdata <- identical(targets, list(integer(0)))

  sco <- numeric(length(v))
  for (i in seq_along(pasets)) {
    s <- NULL
    if (!is.null(cached.scores)) {
        k <- .cached_scores_key(pasets[[i]])
        s <- cached.scores[[i]][[k]]
    }
    if (is.null(s)) {
        if (onlyobsdata) {
          Y <- dat[, v[i]]
          Z <- cbind(1, dat[, pasets[[i]], drop=FALSE])
        } else {
          Y <- dat[global.sufstats$non.int[[i]], v[i]]
          Z <- cbind(1, dat[global.sufstats$non.int[[i]], pasets[[i]],
                     drop=FALSE])
        }
        sigma2 <- sum(Y^2)

        ## scaled error covariance using QR decomposition
        Q <- qr.Q(qr(Z))
        sigma2 <- sigma2 - sum((Y %*% Q)^2)
        lambda <- 0.5 * log(n)
        sco[i] <- -0.5 * global.sufstats$data.count[i] *
                  (1 + log(sigma2 / global.sufstats$data.count[i])) -
                  lambda * (1 + length(pasets[[i]]))
        if (!is.null(cached.scores)) {
            cached.scores[[i]][[k]] <- sco[i]
        }
    } else
        sco[i] <- s
  }

  sum(sco)
}
## assign a name attribute to the iBIC() scoring function for reporting purposes
attr(iBIC, "scorefun.name") <- "iBIC"

## convert a list of targets and a vector of target indices to data
## observations into a logical matrix of observations by variables,
## where TRUE indicates that a variable has been intervened in a observation
.targets2mat <- function(p, targets, target.index) {
  res <- matrix(FALSE, nrow=length(target.index), ncol=p)
  ridx <- rep(seq_along(target.index), lengths(targets)[target.index])
  cidx <- unlist(targets[target.index])
  res[cbind(ridx, cidx)] <- TRUE
  res
}

## cached scores key for a given parent set, computed as sorted parent indices
## glued together with a colon, e.g., '1:2:3'. this type of parent set key is
## used by both, the iBIC() and iBGe() score functions, to store and retrieve
## cached scores for a given parent set in the corresponding per-node
## environment.
.cached_scores_key <- function(paset) {
  k <- paste(sort.int(paset), collapse=":")
  if (nchar(k) == 0L)
    k <- ":"
  k
}

## calculate global sufficient statistics for the iBIC score, which do not
## depend on the structure of a specific DAG, but only on the input data,
## the target vertices and the target indices of the interventions
.iBIC.global.sufstats <- function(dat, targets=list(integer(0)),
                                  target.index=rep(1L, nrow(dat))) {
    stopifnot(is.matrix(dat)) ## QC
    p <- ncol(dat)
    n <- nrow(dat)
    non.int <- NULL
    data.count <- rep(n, p)
    onlyobsdata <- identical(targets, list(integer(0)))
    if (!onlyobsdata) {
        ## index and tally the data points that have not been intervened
        A <- !.targets2mat(p, targets, target.index)
        non.int <- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
        data.count <- colSums(A)
    }
    list(non.int=non.int, data.count=data.count, n=n)
}

## assign the iBIC global sufficient statistics function as an attribute to the
## iBIC() scoring function, so that any search algorithm taking iBIC() as an
## input argument, e.g. scorefun=iBIC, can precompute the corresponding global
## sufficient statistics before iteratively calling iBIC() during search
attr(iBIC, "global.sufstats.fun") <- .iBIC.global.sufstats

#' @importFrom graph numNodes
#' @importFrom cli cli_abort
.check_cached_scores <- function(g, cached.scores) {
  if (!is.null(cached.scores)) {
    if (!is.list(cached.scores) || length(cached.scores) != numNodes(g)) {
      msg <- paste("cached.scores must be a list of length equal to the number",
                   "of nodes in g (", numNodes(g), ")")
      cli_abort(c("x"=msg))
    }
    if (any(vapply(cached.scores, function(x) !is.environment(x), logical(1))))
      cli_abort(c("x"="Each element of cached.scores must be an environment"))
  }
}

#' @importFrom cli cli_abort
.check_g_dat_consistency <- function(g, dat) {
  if (ncol(dat) != numNodes(g))
      cli_abort(c("x"="The number of columns in dat must equal the number of nodes in g"))
  if (is.null(colnames(dat)))
      cli_abort(c("x"="Input data in dat must have column names corresponding to the node names in g"))
  if (!all(nodes(g) %in% colnames(dat)))
      cli_abort(c("x"="All nodes in g must be present as column names in dat"))
  else if (!all(nodes(g) == colnames(dat)))
      cli_abort(c("x"="The order of nodes in g must match the order of column names in dat"))
}

#' @title BGe score for interventional Gaussian data
#'
#' @description Score the goodness-of-fit (GoF) of a given structure of a
#' Bayesian network given an interventional data set of continuous values,
#' where observations are assumed to be independent but not identically
#' distributed (not iid) multivariate Gaussian. This GoF score corresponds to
#' the interventional Bayesian Gaussian equivalent (iBGe) score defined by
#' Kuipers and Moffa (2025). By default, the arguments `targets` and
#' `target.index` are set so that the calculated BIC score assumes there are
#' no interventions in the data.
#' 
#' @param g An acyclic directed graph (DAG) structure of the Bayesian network
#' for which we want to calculate the score.
#'
#' @param dat A `data.frame` object with data records in the rows.
#'
#' @param targets (Default `list(integer(0))`) A `list` object with a family of
#' targets provided as a list of integer vectors. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param target.index (Default a unit vector) A vector of integers in
#' one-to-one correspondence with the rows in `dat`, indicating which rows in
#' the input data are intervened by which targets. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param cached.scores An optional list of environment objects, containing
#' cached scores per parent set for each vertex in `g`. If `NULL` (default),
#' no cached scores are used. Using this argument can speed up the calculation
#' of the score when the same parent sets are scored multiple times. To use
#' this argument, first create an empty environment object with
#' `csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)`
#' and then pass it to this `cached.scores` parameter, i.e.,
#' `cached.scores=csco`. This is currently not implemented for the iBGe score,
#' but it is included as an API placeholder for future versions of the package
#' that will enable this feature for the iBGe score.
#'
#' @param global.sufstats (Default `NULL`) An optional list of global sufficient
#' statistics for the iBGe score, as returned by the `.iBGe.global.sufstats()`
#' function, which do not depend on the structure of a specific DAG, but only
#' on the input data (`dat`), the target vertices (`targets`) and the target
#' indices (`target.index`) of the interventions. If `NULL` (default), the
#' `.iBGe.global.sufstats()` function is internally called. This is currently
#' not implemented for the iBGe score, but it is included as an API placeholder
#' for future versions of the package that will enable this feature for the
#' iBGe score.
#'
#' @return A single numeric value corresponding to the interventional BGe score
#' of the given structure of the Bayesian network for the given data set.
#'
#' @references Kuipers, J. and Moffa, G. The interventional Bayesian Gaussian
#' equivalent score for Bayesian causal inference with unknown soft
#' interventions. *Proceedings of the Fourth Conference on Causal Learning and
#' Reasoning (PMLR)*, 275:772-791, 2025.
#'
#' @examples
#'
#' library(graph)
#'
#' p <- 3
#' nobs <- 100
#' nint <- 100
#' n <- nobs + nint
#'
#' ## define a DAG structure of a Bayesian network with three vertices
#' ## forming a Markov chain X1 -> X2 -> X3
#' g <- new("graphNEL", nodes=c("X1", "X2", "X3"), edgemode="directed")
#' g <- addEdge("X1", "X2", g)
#' g <- addEdge("X2", "X3", g)
#'
#' ## simulate observational data for the previous DAG X1 -> X2 -> X3
#' set.seed(123)
#' X1 <- rnorm(nobs, mean=0, sd=1)
#' X2 <- 0.5 * X1 + rnorm(nobs, mean=0, sd=1)
#' X3 <- 0.5 * X2 + rnorm(nobs, mean=0, sd=1)
#' obsdat <- data.frame(X1=X1, X2=X2, X3=X3)
#'
#' ## simulate interventional data for the same DAG, where X2 is intervened
#' X1 <- rnorm(nint, mean=0, sd=1)
#' X2 <- rnorm(nint, mean=0, sd=1) + 1.0
#' X3 <- 0.5 * X2 + rnorm(nint, mean=0, sd=1)
#' intdat <- data.frame(X1=X1, X2=X2, X3=X3)
#'
#' ## combine observational and interventional data
#' dat <- rbind(obsdat, intdat)
#'
#' ## define the targets and target indices for the interventional data
#' targets <- list(integer(0), 2L)
#' target.index <- c(rep(1L, nobs), rep(2L, nint))
#'
#' ## calculate the interventional BGe score for the DAG and data set
#' iBGe(g, dat, targets, target.index)
#'
#' ## create another Markov equivalent DAG by reversing the arc X1 -> X2
#' ## to obtain X1 <- X2 -> X3
#' g2 <- g
#' g2 <- removeEdge("X1", "X2", g2)
#' g2 <- addEdge("X2", "X1", g2)
#'
#' ## calculate the interventional BGe score for the new DAG on the
#' ## same data, notice that the score is different despite being a
#' ## Markov equivalent DAG
#' iBGe(g2, dat, targets, target.index)
#'
#' ## this is not the case if we do not indicate the presence of interventions
#' ## in the data
#' iBGe(g, dat)
#' iBGe(g2, dat)
#'
#' @importFrom methods as
#' @importFrom graph numNodes edgeMatrix
#' @export
iBGe <- function(g, dat, targets=list(integer(0)),
                 target.index=rep(1L, nrow(dat)),
                 cached.scores=NULL, global.sufstats=NULL) {

  if (is.null(attr(dat, "sanitycheck"))) {
    dat <- .check_input_data(dat)
    .check_g_dat_consistency(g, dat)
  }

  v <- match(nodes(g), colnames(dat))
  em <- edgeMatrix(g)
  pasets <- split(em["from", ], factor(v[em["to", ]], levels=v))
  stopifnot(identical(names(pasets), as.character(v)))
  .check_cached_scores(g, cached.scores)

  if (is.null(global.sufstats))
    global.sufstats <- .iBGe.global.sufstats(dat, targets, target.index)

  sco <- numeric(length(v))
  for (i in seq_along(pasets)) {
    s <- NULL
    if (!is.null(cached.scores)) {
      k <- .cached_scores_key(pasets[[i]])
      s <- cached.scores[[i]][[k]]
    }
    if (is.null(s)) {
      TNj <- global.sufstats$TN[[i]]
      lp <- length(pasets[[i]])
      A <- TNj[i, i]
      awpNd2 <- (global.sufstats$awpN[i] - global.sufstats$p + lp + 1) / 2
      if (lp == 0L)
        s <- global.sufstats$scoreconstvec[[i]][1L] - awpNd2 * log(A)
      else {
        D <- TNj[pasets[[i]], pasets[[i]], drop=FALSE]
        R <- chol(D)
        logdetD <- 2 * sum(log(diag(R)))
        B <- TNj[i, pasets[[i]]]
        logdetpart2 <- log(A - sum(backsolve(R, B, transpose=TRUE)^2))
        s <- global.sufstats$scoreconstvec[[i]][lp + 1L] -
             awpNd2 * logdetpart2 - logdetD / 2
      }
      if (!is.null(cached.scores))
        cached.scores[[i]][[k]] <- s
    }
    sco[i] <- s
  }

  sum(sco)
}
## assign a name attribute to the iBGe() scoring function for reporting purposes
attr(iBGe, "scorefun.name") <- "iBGe"

## calculate global sufficient statistics for the iBGe score, which do not
## depend on the structure of a specific DAG, but only on the input data,
## the target vertices and the target indices of the interventions. part of
## this code is adapted from the BGe parametrisation in Kuipers & Moffa (2025)
## and the BiDAG package, but stripped down to exclude BDe, BDecat, DBN, MDAG,
## and other stuff not exposed in the iBGe() function of this package
.iBGe.global.sufstats <- function(dat, targets=list(integer(0)),
                                  target.index=rep(1L, nrow(dat))) {
  stopifnot(is.matrix(dat)) ## QC
  p <- ncol(dat)
  n <- nrow(dat)

  ## BGe equivalent sample size for the prior distribution of the mean vector
  ## set to 1, which assigns the weakest possible informative weight to this
  ## prior distribution
  ## see BiDAG::scoreparameters for further details on this parameter
  ## we might want to expose this as a user parameter in the future
  am <- 1

  ## BGe edge penalization factor, set to 1 (no penalization)
  ## see BiDAG::scoreparameters for further details on this parameter
  ## we might want to expose this as a user parameter in the future
  edgepf <- 1

  aw <- p + am + 1
  T0scale <- am * (aw - p - 1) / (am + 1) # follows from [GH2002, eqs. (19, 20)]
  T0 <- diag(T0scale, p, p)
  logedgepf <- log(edgepf)
  l <- seq_len(p) # l = number of parents + 1

  non.int <- NULL
  data.count <- rep(n, p)
  onlyobsdata <- identical(targets, list(integer(0)))
  if (!onlyobsdata) {
    ## index and tally the data points that have not been intervened
    A <- !.targets2mat(p, targets, target.index)
    non.int <- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
    data.count <- colSums(A)
  }

  TN <- vector("list", p)
  awpN <- numeric(p)
  scoreconstvec <- vector("list", p)
  for (j in seq_len(p)) {
    Xj <- dat
    if (!onlyobsdata)
      Xj <- dat[non.int[[j]], , drop=FALSE]
    Nj <- data.count[j]
    if (Nj < 2) {
      msg <- paste("Not enough observational input data in column number", j,
                   "(", Nj, "observed values)")
      cli_abort(c("x"=msg))
    }
    means <- colMeans(Xj)
    covmat <- cov(Xj) * (Nj - 1)
    TN[[j]] <- T0 + covmat + (am * Nj / (am + Nj)) * outer(means, means)
    awpN[j] <- aw + Nj
    constscorefact <- -(Nj / 2) * log(pi) + 0.5 * log(am / (am + Nj))
    awp <- aw - p + l
    scoreconstvec[[j]] <- constscorefact - lgamma(awp / 2) + lgamma((awp + Nj) / 2) +
                          ((awp + l - 1) / 2) * log(T0scale) - l * logedgepf
  }

  list(p=p, aw=aw, T0scale=T0scale, TN=TN, awpN=awpN,
       scoreconstvec=scoreconstvec, non.int=non.int, data.count=data.count, n=n)
}

## assign the iBGe global sufficient statistics function as an attribute to the
## iBGe() scoring function, so that any search algorithm taking iBGe() as an
## input argument, e.g. scorefun=iBGe, can precompute the corresponding global
## sufficient statistics before iteratively calling iBGe() during search
attr(iBGe, "global.sufstats.fun") <- .iBGe.global.sufstats

## first original version of the iBGe() function, which calls the vendored code
## of the iBGe score by Kuipers and Moffa (2025) based and adapted from the
## scripts provided at https://github.com/jackkuipers/iBGe and from the BiDAG
## package at https://cran.r-project.org/package=BiDAG this is included here
## to verify that further optimized versions of the iBGe() function produce the
## same results as the original version of the iBGe score by Kuipers and Moffa
.vendored_iBGe <- function(g, dat, targets=list(integer(0)),
                           target.index=rep(1L, nrow(dat))) { # nocov start

  .check_g_dat_consistency(g, dat)
  v <- nodes(g)
  p <- numNodes(g)
  n <- nrow(dat)
  em <- edgeMatrix(g)
  pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
  stopifnot(identical(names(pasets), as.character(v)))

  ## create intervention matrix for BiDAG
  A <- .targets2mat(p, targets, target.index)
  I <- matrix(0, nrow=n, ncol=p)
  I[A] <- 1
  A <- as(as(g, "graphAM"), "matrix")

  param <- .scoreparameters(scoretype="usr", data=dat,
                            usrpar=list(pctesttype="bge", Tmat=I))
  .DAGscore(param, A)
} # nocov end

## the code below has been copied and adapted from
## https://github.com/jackkuipers/iBGe and the BiDAG package at
## https://cran.r-project.org/package=BiDAG to enable calling it from the
## idlBNs package

### These user defined score functions are for known perfect interventions

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam,
                               usrpar = list(Tmat = NULL, pctesttype = "bge",
                                             am = 1, chi = 1, edgepf = 1,
                                             edgepmat = NULL)) { # nocov start
  n <- initparam$n
  Tmat <- usrpar$Tmat
  nodeparams <- vector("list", n)
  for (jj in 1:n){
    nint_obs <- which(Tmat[, jj] == 0)
    if (length(nint_obs) < 2) {
      stop("Not enough observational data.")
    } else {
      nodeparams[[jj]] <- .scoreparameters(scoretype = usrpar$pctesttype,
                                           data = initparam$data[nint_obs, ],
                                           weightvector = initparam$weightvector[nint_obs],
                                           bgepar = list(am = usrpar$am),
                                           bdepar = list(chi = usrpar$chi, edgepf = usrpar$edgepf),
                                           bdecatpar = list(chi = usrpar$chi, edgepf = usrpar$edgepf),
                                           edgepmat = usrpar$edgepmat)
    }
  }
  initparam$nodeparams <- nodeparams

  initparam
} # nocov end

## here we have put only the BGe part

#' @importFrom stats cov cov.wt
.scoreparameters <- function(scoretype=c("bge","bde","bdecat","usr"), data,
                          bgepar=list(am=1, aw=NULL, edgepf=1), bdepar=list(chi=0.5, edgepf=2),
                          bdecatpar=list(chi=0.5, edgepf=2), dbnpar=list(samestruct=TRUE,
                            slices=2, b=0, stationary=TRUE, rowids=NULL, datalist=NULL,
                            learninit=TRUE), usrpar=list(pctesttype=c("bge","bde","bdecat")),
                          mixedpar=list(nbin=0), MDAG=FALSE, DBN=FALSE, weightvector=NULL,
                          bgnodes=NULL, edgepmat=NULL, nodeslabels=NULL) { # nocov start

  initparam<-list()

  bgn<-length(bgnodes)
  n <- ncol(data)
  nsmall<-n-bgn #number of nodes in the network excluding root nodes

  if (ncol(data)!=nsmall+bgn)
    stop("n and the number of columns in the data do not match")

  if (!is.null(weightvector)) {
    if (length(weightvector)!=nrow(data)) {
      stop("Length of the weightvector does not match the number of rows (observations) in data")
    }
  }

  if (is.null(nodeslabels)) {
    if(all(is.character(colnames(data)))){
      nodeslabels<-colnames(data)
    } else {
      nodeslabels<-sapply(c(1:n), function(x)paste("v",x,sep=""))
    }
  }

  multwv<-NULL

  if (is.null(dbnpar$datalist)) colnames(data)<-nodeslabels

  initparam$labels<-nodeslabels
  initparam$type<-scoretype
  initparam$DBN<-DBN
  initparam$MDAG<-MDAG
  initparam$weightvector<-weightvector
  initparam$data<-data

  initparam$bgnodes<-bgnodes
  initparam$static<-bgnodes
  if(!is.null(bgnodes)) {
    initparam$mainnodes<-c(1:n)[-bgnodes]
  } else initparam$mainnodes<-c(1:n)

  initparam$bgn<-bgn
  initparam$n<-n
  initparam$nsmall<-nsmall

  initparam$labels.short<-initparam$labels

  if (is.null(edgepmat)) {
    initparam$logedgepmat <- NULL
  } else {
    if(all(edgepmat>0)) {
    initparam$logedgepmat <- log(edgepmat)
    } else
      stop("all entries of edgepmat matrix must be bigger than 0! 1 corresponds to no penalization")
  }

  if (scoretype == "bge") {

    if(is.null(bgepar$am)) {
      bgepar$am<-1
    }
    if(is.null(bgepar$aw)) {
      bgepar$aw<-n+bgepar$am+1
    }
    if(is.null(bgepar$edgepf)) {
      bgepar$edgepf<-1
    }

    if (is.null(weightvector)) {
      N<-nrow(data)
      covmat<-cov(data)*(N-1)
      means<-colMeans(data)
    } else {
      N<-sum(weightvector)
      forcov<-cov.wt(data,wt=weightvector,cor=TRUE,method="ML")
      covmat<-forcov$cov*N
      means<-forcov$center
    }

    initparam$am <- bgepar$am # store parameters
    initparam$aw <- bgepar$aw
    initparam$pf <- bgepar$edgepf

    initparam$N <- N # store effective sample size
    #initparam$covmat <- (N-1)*covmat
    initparam$means <- means # store means

    mu0<-numeric(n)
    #https://arxiv.org/pdf/1302.6808.pdf page 10
    T0scale <- bgepar$am*(bgepar$aw-n-1)/(bgepar$am+1) # This follows from equations (19) and (20) of [GH2002]
    T0<-diag(T0scale,n,n)
    initparam$TN <- T0 + covmat + ((bgepar$am*N)/(bgepar$am+N))* (mu0 - means)%*%t(mu0 - means)
    initparam$awpN<-bgepar$aw+N
    constscorefact<- -(N/2)*log(pi) + (1/2)*log(bgepar$am/(bgepar$am+N))

    initparam$muN <- (N*means + bgepar$am*mu0)/(N + bgepar$am) # posterior mean mean
    initparam$SigmaN <- initparam$TN/(initparam$awpN-n-1) # posterior mode covariance matrix

    initparam$scoreconstvec<-numeric(n)
    for (j in (1:n)) {# j represents the number of parents plus 1
      awp<-bgepar$aw-n+j
      initparam$scoreconstvec[j]<-constscorefact - lgamma(awp/2) + lgamma((awp+N)/2) + ((awp+j-1)/2)*log(T0scale) - j*log(initparam$pf)
    }

  } else if (scoretype == "usr") { ## usr
    if(is.null(usrpar$pctesttype)){usrpar$pctesttype <- "usr"}
    initparam$pctesttype <- usrpar$pctesttype
    initparam <- usrscoreparameters(initparam, usrpar)
  } else
    stop("not supported in this package.")

  attr(initparam, "class") <- "scoreparameters"
  return(initparam)
} # nocov end

### This function evaluates the log score of a node given its parents

.usrDAGcorescore <- function (j, parentnodes, n, param) { # nocov start
  .DAGcorescore(j, parentnodes, n, param$nodeparams[[j]])
} # nocov end

.DAGscore <- function(scorepar, incidence){ # nocov start
  if(scorepar$DBN) {
    stop("To calculate DBN score DBNscore should be used!")
  }
  n<-ncol(scorepar$data)
  if(scorepar$bgn==0) {
    mainnodes<-c(1:scorepar$n)
  } else {
    mainnodes<-c(1:n)[-scorepar$bgnodes]
  }
  P_local <- numeric(n)
  for (j in mainnodes)  { #j is a node at which scoring is done
    parentnodes <- which(incidence[,j]==1)
    P_local[j]<-.DAGcorescore(j,parentnodes,scorepar$n,scorepar)
  }
  return(sum(P_local))
} # nocov end


# The log of the BGe/BDe score, but simplified as much as possible
# see arXiv:1402.6863 
.DAGcorescore<-function(j,parentnodes,n,param) { # nocov start

  if (param$type=="bge") {
    TN<-param$TN
    awpN<-param$awpN
    scoreconstvec<-param$scoreconstvec
    
    lp<-length(parentnodes) #number of parents
    awpNd2<-(awpN-n+lp+1)/2
    A<-TN[j,j]
    switch(as.character(lp),
           "0"={# just a single term if no parents
             corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
           },
           
           "1"={# no need for matrices
             D<-TN[parentnodes,parentnodes]
             logdetD<-log(D)
             B<-TN[j,parentnodes]
             logdetpart2<-log(A-B^2/D)
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - param$logedgepmat[parentnodes, j]
             }
           },
           
           "2"={# can do matrix determinant and inverse explicitly
             # but this is numerically unstable for large matrices!
             # so we use the same approach as for 3 parents
             D<-TN[parentnodes,parentnodes]
             detD<-dettwobytwo(D)
             logdetD<-log(detD)
             B<-TN[j,parentnodes]
             #logdetpart2<-log(A-(D[2,2]*B[1]^2+D[1,1]*B[2]^2-2*D[1,2]*B[1]*B[2])/detD) #also using symmetry of D
             logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
             }
           },
           
           {# otherwise we use cholesky decomposition to perform both
             D<-as.matrix(TN[parentnodes,parentnodes])
             choltemp<-chol(D)
             logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
             B<-TN[j,parentnodes]
             logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
             }
           })

  } else if (param$type=="usr") {
    corescore <- .usrDAGcorescore(j,parentnodes,n,param)
  } 
  
  return(corescore)
} # nocov end

# The determinant of a 2 by 2 matrix
dettwobytwo <- function(D) { # nocov start
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
} # nocov end
