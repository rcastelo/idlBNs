test_scores <- function() {
  library(graph)

  p <- 3
  nobs <- 100
  nint <- 100
  n <- nobs + nint
 
  ## define a DAG structure of a Bayesian network with three vertices
  ## forming a Markov chain X1 -> X2 -> X3
  g <- new("graphNEL", nodes=c("X1", "X2", "X3"), edgemode="directed")
  g <- addEdge("X1", "X2", g)
  g <- addEdge("X2", "X3", g)
 
  ## simulate observational data for the previous DAG X1 -> X2 -> X3
  set.seed(123)
  X1 <- rnorm(nobs, mean=0, sd=1)
  X2 <- 0.5 * X1 + rnorm(nobs, mean=0, sd=1)
  X3 <- 0.5 * X2 + rnorm(nobs, mean=0, sd=1)
  obsdat <- data.frame(X1=X1, X2=X2, X3=X3)
 
  ## simulate interventional data for the same DAG, where X2 is intervened
  X1 <- rnorm(nint, mean=0, sd=1)
  X2 <- rnorm(nint, mean=0, sd=1) + 1.0
  X3 <- 0.5 * X2 + rnorm(nint, mean=0, sd=1)
  intdat <- data.frame(X1=X1, X2=X2, X3=X3)
 
  ## combine observational and interventional data
  dat <- rbind(obsdat, intdat)
 
  ## define the targets and target indices for the interventional data
  targets <- list(0L, 2L)
  target.index <- c(rep(1L, nobs), rep(2L, nint))
 
  ## calculate the interventional BIC score for the DAG and data set
  ibic <- iBIC(g, dat, targets, target.index)
 
  ## create another Markov equivalent DAG by reversing the arc X1 -> X2
  ## to obtain X1 <- X2 -> X3
  g2 <- g
  g2 <- removeEdge("X1", "X2", g2)
  g2 <- addEdge("X2", "X1", g2)
 
  ## calculate the interventional BIC score for the new DAG on the
  ## same data, the score should be different despite being a
  ## Markov equivalent DAG and with a lower, more negative, value
  ibic2 <- iBIC(g2, dat, targets, target.index)

  checkTrue(ibic > ibic2)
 
  ## in the case when we do not indicate the presence of interventions
  ## in the data, both scores should be identical
  ibic <- iBIC(g, dat)
  ibic2 <- iBIC(g2, dat)
  checkEqualsNumeric(ibic, ibic2)

  ## calculate the interventional BGe score for the DAG and data set
  ibge <- iBGe(g, dat, targets, target.index)

  ## calculate the interventional BGe score for the new DAG on the
  ## same data, the score should be different despite being a
  ## Markov equivalent DAG and with a lower, more negative, value
  ibge2 <- iBGe(g2, dat, targets, target.index)

  checkTrue(ibge > ibge2)

  ## in the case when we do not indicate the presence of interventions
  ## in the data, both scores should be identical
  ibge <- iBGe(g, dat)
  ibge2 <- iBGe(g2, dat)
  checkEqualsNumeric(ibge, ibge2)
}

test_search_algorithms <- function() {
  suppressPackageStartupMessages({
    library(graph)
    library(pcalg)
  })

  p <- 5
  k <- 2
  n <- 30
 
  ## simulate a random DAG
  set.seed(123)
  Mg <- r.gauss.pardag(p, 0.6, top.sort=TRUE, normalize=TRUE)
  g <- as(Mg, "graphNEL")
  e <- dag2essgraph(as(g, "graphNEL"))
  
  ## generate a random family of intervention targets
  I <- c(list(integer(0)), sample(p, size=k, replace=FALSE))
 
  ## sample size per different target (including the no-target)
  nbytgts <- rep(floor(n / (k + 1)), k)
  nbytgts <- c(n - sum(nbytgts), nbytgts)
 
  ## simulate mixed observational and interventional data
  dat <- list()
  for (v in seq_along(I)) {
    targets <- I[[v]]
    dat[[v]] <- rmvnorm.ivent(nbytgts[v], Mg, target=targets,
                              target.value=rep(2, length(targets)))
  }
  dat <- do.call("rbind", dat)
 
  ## store the target index for each row of the data
  tindex <- rep(1:length(nbytgts), nbytgts)
 
  ## run the HCMC algorithm assuming all data were observational
  dhat.hcmc <- hcmc(dat, verbose=TRUE)
 
  ## calculate the structural Hamming distance (SHD) between the generative
  ## DAG and the estimated DAG
  shd.dhat.hcmc <- shd(e, dag2essgraph(dhat.hcmc$dag))
 
  ## run the iHCMC algorithm informing the presence of interventional data
  ## using by the default the interventional BIC score (see the iBIC()
  ## function).
  dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex, verbose=FALSE)

  ## the estimated DAG is closer to the generative DAG (lower SHD value)
  ## than the one estimated by HCMC, which did not take into account the
  ## presence of interventional data
  shd.dhat.ihcmc <- shd(e, dag2essgraph(dhat.ihcmc$dag))
 
  checkTrue(dhat.ihcmc$sco > dhat.hcmc$sco)
 
  checkTrue(shd.dhat.ihcmc < shd.dhat.hcmc)

  ## compare with the classical hill-climbing algorithm
  dhat.hc <- hillclimbing(dat)
  shd.dhat.hc <- shd(e, dag2essgraph(dhat.hc$dag))

  checkTrue(dhat.ihcmc$sco > dhat.hc$sco)
 
  checkTrue(shd.dhat.ihcmc < shd.dhat.hc)
}
