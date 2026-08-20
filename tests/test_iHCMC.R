## 2026-08-20 regression test for iHCMC

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
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
set.seed(123)
dhat.hcmc <- hcmc(dat, verbose=FALSE)
dhat.hcmc

## calculate the structural Hamming distance (SHD) between the generative
## DAG and the estimated DAG
shd(e, dag2essgraph(dhat.hcmc$dag))

## run the iHCMC algorithm informing the presence of interventional data
## using by the default the interventional BIC score (see the iBIC()
## function).
set.seed(123)
dhat.ihcmc <- hcmc(dat, targets=I, target.index=tindex, verbose=FALSE)
dhat.ihcmc

## the estimated DAG is closer to the generative DAG (lower SHD value)
## than the one estimated by HCMC, which did not take into account the
## presence of interventional data
shd(e, dag2essgraph(dhat.ihcmc$dag))

## run it again this time using the interventional BGe score (see the
## iBGe() function), which provides an estimate closer to the generative DAG
set.seed(123)
dhat.ihcmc2 <- hcmc(dat, targets=I, target.index=tindex, scorefun=iBGe,
                    verbose=FALSE)
dhat.ihcmc2
shd(e, dag2essgraph(dhat.ihcmc2$dag))
