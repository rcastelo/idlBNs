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
#' @param targets A `list` object with a family of targets.
#'
#' @param target.index A vector of integer values in one-to-one correspondence
#' with the rows in `dat`, indicating what set of targets has generated each
#' row in `dat`.
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
#' targets <- list(0L, 2L)
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
iBIC <- function(g, dat, targets=list(0L),
                 target.index=rep(1L, nrow(dat))) {
  v <- nodes(g)
  p <- numNodes(g)
  n <- nrow(dat)
  em <- edgeMatrix(g)
  pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
  stopifnot(identical(names(pasets), v))

  onlyobsdata <- identical(targets, list(0L))
  data.count <- rep(n, p)
  if (!onlyobsdata) {
    ## index of the data points that have not been intervened per vertex
    A <- !.targets2mat(p, targets, target.index)
    non.int <- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
    data.count <- colSums(A)
  }
  sco <- numeric(length(v))
  for (i in seq_along(pasets)) {
    if (onlyobsdata) {
      Y <- dat[, v[i]]
      Z <- cbind(1, dat[, pasets[[i]], drop=FALSE])
    } else {
      Y <- dat[non.int[[i]], v[i]]
      Z <- cbind(1, dat[non.int[[i]], pasets[[i]], drop=FALSE])
    }
    sigma2 <- sum(Y^2)

    ## scaled error covariance using QR decomposition
    Q <- qr.Q(qr(Z))
    sigma2 <- sigma2 - sum((Y %*% Q)^2)
    lambda <- 0.5 * log(n)
    sco[i] <- -0.5 * data.count[i] * (1 + log(sigma2 / data.count[i])) -
                                      lambda * (1 + length(pasets[[i]]))
  }
  sum(sco)
}

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

#' @title BGe score for interventional Gaussian data
#'
#' @description Score the goodness-of-fit (GoF) of a given structure of a
#' Bayesian network given an interventional data set of continuous values,
#' where observations are assumed to be independent but not identically
#' distributed (not iid) multivariate Gaussian. This GoF score corresponds to
#' the interventional Bayesian Gaussian equivalent (iBGe) score defined by
#' Kuipers and Moffa (2022). By default, the arguments `targets` and
#' `target.index` are set so that the calculated BIC score assumes there are
#' no interventions in the data.
#' 
#' @param g An acyclic directed graph (DAG) structure of the Bayesian network
#' for which we want to calculate the score.
#'
#' @param dat A `data.frame` object with data records in the rows.
#'
#' @param targets A `list` object with a family of targets.
#'
#' @param target.index A vector of integer values in one-to-one correspondence
#' with the rows in `dat`, indicating what set of targets has generated each
#' row in `dat`.
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
#' targets <- list(0L, 2L)
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
iBGe <- function(g, dat, targets=list(0L), target.index=rep(1L, nrow(dat))) {
  v <- nodes(g)
  p <- numNodes(g)
  n <- nrow(dat)
  em <- edgeMatrix(g)
  pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
  stopifnot(identical(names(pasets), v))

  ## create intervention matrix for BiDAG
  A <- .targets2mat(p, targets, target.index)
  I <- matrix(0, nrow=n, ncol=p)
  I[A] <- 1
  A <- as(as(g, "graphAM"), "matrix")

  param <- .scoreparameters(scoretype="usr", data=dat,
                            usrpar=list(pctesttype="bge", Tmat=I))
  .DAGscore(param, A)
}

## the code below has been copied and adapted from
## https://github.com/jackkuipers/iBGe and the BiDAG package at
## https://cran.r-project.org/package=BiDAG to enable calling it from the
## idlBNs package

### These user defined score functions are for known perfect interventions

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam,
                               usrpar = list(Tmat = NULL, pctesttype = "bge",
                                             am = 1, chi = 1, edgepf = 1,
                                             edgepmat = NULL)) {
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
}

## here we have put only the BGe part

#' @importFrom stats cov cov.wt
.scoreparameters <- function(scoretype=c("bge","bde","bdecat","usr"), data,
                          bgepar=list(am=1, aw=NULL, edgepf=1), bdepar=list(chi=0.5, edgepf=2),
                          bdecatpar=list(chi=0.5, edgepf=2), dbnpar=list(samestruct=TRUE,
                            slices=2, b=0, stationary=TRUE, rowids=NULL, datalist=NULL,
                            learninit=TRUE), usrpar=list(pctesttype=c("bge","bde","bdecat")),
                          mixedpar=list(nbin=0), MDAG=FALSE, DBN=FALSE, weightvector=NULL,
                          bgnodes=NULL, edgepmat=NULL, nodeslabels=NULL) {

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
}

### This function evaluates the log score of a node given its parents

.usrDAGcorescore <- function (j, parentnodes, n, param) {
  .DAGcorescore(j, parentnodes, n, param$nodeparams[[j]])
}

.DAGscore <- function(scorepar, incidence){
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
}


# The log of the BGe/BDe score, but simplified as much as possible
# see arXiv:1402.6863 
.DAGcorescore<-function(j,parentnodes,n,param) {

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
}

# The determinant of a 2 by 2 matrix
dettwobytwo <- function(D) {
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}
