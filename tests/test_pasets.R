## 2026-09-05 regression test for incremental parent-set (pasets) maintenance

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

tol <- .Machine$double.eps^0.5

## helper to compare two pasets lists ignoring outer/inner names and order
## within each parent set
pasets_equal <- function(a, b) {
  a <- unname(lapply(a, function(x) unname(sort.int(x))))
  b <- unname(lapply(b, function(x) unname(sort.int(x))))
  identical(a, b)
}

################################################################################
## 1. property test: for every move returned by nr.nh()/ar.nh()/ncr.nh(),
## the pasets derived by applying the move's own (op, u, v) delta to the
## current DAG's pasets must equal a from-scratch rebuild on the graph that
## apply.move() produces for that same move. also checks that the three
## neighborhoods enumerate exactly the move sets an independent brute-force
## reference produces (see nh_reference() below), since the neighborhoods no
## longer return candidate graphs to check against
################################################################################

## independent brute-force reference for the three neighborhoods, written
## against the DAG itself (adjacency + reachability recomputed from scratch)
## rather than against the incrementally maintained ancestor matrix, so it
## does not share any machinery with the functions under test. returns a
## character vector of "op:u:v" move descriptors.
nh_reference <- function(dag, kind=c("nr", "ar", "ncr"), utargets=integer(0)) {
  kind <- match.arg(kind)
  v <- nodes(dag)
  p <- length(v)
  adj <- matrix(FALSE, p, p)
  em <- edgeMatrix(dag)
  if (ncol(em) > 0)
    adj[cbind(em["from", ], em["to", ])] <- TRUE
  ## transitive closure by repeated squaring (reach[a, b]: a path a -> b)
  reach <- adj
  repeat {
    nxt <- reach | ((reach %*% adj) > 0)
    if (identical(nxt, reach)) break
    reach <- nxt
  }
  pa <- lapply(seq_len(p), function(j) which(adj[, j]))
  ## a move is legal iff the resulting adjacency matrix has no directed cycle
  acyclic <- function(A) {
    R <- A
    repeat {
      N <- R | ((R %*% A) > 0)
      if (identical(N, R)) break
      R <- N
    }
    !any(diag(R))
  }
  mv <- character(0)
  for (i in seq_len(p)) {          ## additions, in i-major then j order
    for (j in seq_len(p)) {
      if (i == j || adj[i, j]) next
      A <- adj; A[i, j] <- TRUE
      if (acyclic(A)) mv <- c(mv, paste("1", i, j, sep=":"))
    }
  }
  ## nr.nh() emits, per vertex i, that vertex's additions then its removals
  mv2 <- character(0)
  for (i in seq_len(p)) {
    mv2 <- c(mv2, grep(paste0("^1:", i, ":"), mv, value=TRUE))
    for (j in which(adj[i, ])) mv2 <- c(mv2, paste("2", i, j, sep=":"))
  }
  mv <- mv2
  if (kind == "nr") return(mv)
  for (i in seq_len(p)) {          ## reversals
    for (j in which(adj[i, ])) {
      A <- adj; A[i, j] <- FALSE; A[j, i] <- TRUE
      if (!acyclic(A)) next
      if (kind == "ncr") {
        covered <- setequal(pa[[i]], setdiff(pa[[j]], i))
        touches <- (i %in% utargets) || (j %in% utargets)
        if (covered && !touches) next
      }
      mv <- c(mv, paste("3", i, j, sep=":"))
    }
  }
  mv
}

nh_moves <- function(ne) paste(ne$op, ne$u, ne$v, sep=":")

set.seed(42)
for (p in c(4, 6, 8)) {
  varnames <- paste0("X", seq_len(p))
  n <- 50
  dat <- matrix(rnorm(n * p), nrow=n, ncol=p, dimnames=list(NULL, varnames))

  ## build a handful of random DAGs (in topological order, no cycles) with a
  ## moderate edge density, and for each one check every candidate neighbor
  for (rep in 1:5) {
    dag <- new("graphNEL", nodes=varnames, edgemode="directed")
    for (i in seq_len(p - 1))
      for (j in (i + 1):p)
        if (runif(1) < 0.3)
          dag <- addEdge(varnames[i], varnames[j], dag)

    anc <- idlBNs:::init.ancestors(varnames)
    ## rebuild anc consistent with dag's actual edges via repeated add.ancestors
    em <- edgeMatrix(dag)
    for (k in seq_len(ncol(em)))
      anc <- idlBNs:::add.ancestors(anc, varnames[em["from", k]], varnames[em["to", k]])

    vidx <- setNames(seq_len(p), varnames)
    pasets <- idlBNs:::.build_pasets(dag, dat)

    ## exercise both the no-intervention and the with-targets branch of
    ## ncr.nh()'s I-covered arc test
    for (utargets in list(integer(0), c(1L, min(3L, p)))) {
      nhs <- list(nr=idlBNs:::nr.nh(dag, anc),
                  ar=idlBNs:::ar.nh(dag, anc),
                  ncr=idlBNs:::ncr.nh(dag, anc, utargets))

      for (kind in names(nhs)) {
        ne <- nhs[[kind]]
        ## the enumerated move set matches the independent reference exactly,
        ## including its order (which.max() breaks ties by position)
        stopifnot(identical(nh_moves(ne),
                            nh_reference(dag, kind, utargets)))

        for (m in seq_along(ne$op)) {
          uu <- vidx[[varnames[ne$u[m]]]]
          vv <- vidx[[varnames[ne$v[m]]]]
          derived <- switch(ne$op[m],
                            idlBNs:::add.pasets(pasets, uu, vv),
                            idlBNs:::remove.pasets(pasets, uu, vv),
                            idlBNs:::reverse.pasets(pasets, uu, vv))
          gm <- idlBNs:::apply.move(dag, ne$op[m], ne$u[m], ne$v[m], varnames)
          truth <- idlBNs:::.build_pasets(gm, dat)
          stopifnot(pasets_equal(derived, truth))
        }
      }
    }
  }
}

################################################################################
## 2. iBIC()/iBGe() 'pasets=' parity: a precomputed pasets list must give
## identical scores to the default from-scratch rebuild
################################################################################

p <- 3
nobs <- 100
nint <- 100
n <- nobs + nint

g <- new("graphNEL", nodes=c("X1", "X2", "X3"), edgemode="directed")
g <- addEdge("X1", "X2", g)
g <- addEdge("X2", "X3", g)

set.seed(123)
X1 <- rnorm(nobs, mean=0, sd=1)
X2 <- 0.5 * X1 + rnorm(nobs, mean=0, sd=1)
X3 <- 0.5 * X2 + rnorm(nobs, mean=0, sd=1)
obsdat <- data.frame(X1=X1, X2=X2, X3=X3)

X1 <- rnorm(nint, mean=0, sd=1)
X2 <- rnorm(nint, mean=0, sd=1) + 1.0
X3 <- 0.5 * X2 + rnorm(nint, mean=0, sd=1)
intdat <- data.frame(X1=X1, X2=X2, X3=X3)

dat <- rbind(obsdat, intdat)
targets <- list(0L, 2L)
target.index <- c(rep(1L, nobs), rep(2L, nint))

datm <- as.matrix(dat)
pre_pasets <- idlBNs:::.build_pasets(g, datm)

s_default <- iBIC(g, dat, targets, target.index)
s_pasets  <- iBIC(g, dat, targets, target.index, pasets=pre_pasets)
stopifnot(abs(s_default - s_pasets) < tol)

s_default_ge <- iBGe(g, dat, targets, target.index)
s_pasets_ge  <- iBGe(g, dat, targets, target.index, pasets=pre_pasets)
stopifnot(abs(s_default_ge - s_pasets_ge) < tol)

## also across a reversed DAG and pure observational data
g2 <- g
g2 <- removeEdge("X1", "X2", g2)
g2 <- addEdge("X2", "X1", g2)
pre_pasets2 <- idlBNs:::.build_pasets(g2, datm)

s2_default <- iBIC(g2, dat, targets, target.index)
s2_pasets  <- iBIC(g2, dat, targets, target.index, pasets=pre_pasets2)
stopifnot(abs(s2_default - s2_pasets) < tol)

sobs_default <- iBIC(g, dat)
sobs_pasets  <- iBIC(g, dat, pasets=idlBNs:::.build_pasets(g, datm))
stopifnot(abs(sobs_default - sobs_pasets) < tol)

## a mismatched-length pasets must still be rejected
res <- tryCatch({
    iBIC(g, dat, targets, target.index, pasets=pre_pasets[-1])
    "no error"
}, error=function(e) "error")
stopifnot(identical(res, "error"))

################################################################################
## 3. in-loop debug-assertion mode: confirm the incrementally-maintained
## pasets stays identical to a from-scratch rebuild at every iteration, for
## hillclimbing() (single sync point) and hcmc() at both r=0 (accept-branch
## sync point only) and the default r (also exercising rcar()'s two sync
## points), with both iBIC and iBGe, and with interventional data
################################################################################

old_opt <- getOption("idlBNs.debug.pasets", FALSE)
on.exit(options(idlBNs.debug.pasets=old_opt), add=TRUE)
options(idlBNs.debug.pasets=TRUE)

set.seed(1)
p <- 8; n <- 60
D <- r.gauss.pardag(p, 0.3, top.sort=TRUE, normalize=TRUE)
dat3 <- rmvnorm.ivent(n, D, target=integer(0), target.value=numeric(0))
targets3 <- list(integer(0), c(1L, 3L))
tindex3 <- c(rep(1L, n / 2), rep(2L, n - n / 2))

invisible(hillclimbing(dat3, scorefun=iBIC, verbose=FALSE))
invisible(hillclimbing(dat3, scorefun=iBGe, verbose=FALSE))
invisible(hcmc(dat3, r=0, scorefun=iBIC, verbose=FALSE))
invisible(hcmc(dat3, scorefun=iBIC, verbose=FALSE))
invisible(hcmc(dat3, scorefun=iBGe, verbose=FALSE))
invisible(hcmc(dat3, targets=targets3, target.index=tindex3, scorefun=iBIC, verbose=FALSE))
invisible(hcmc(dat3, targets=targets3, target.index=tindex3, scorefun=iBGe, verbose=FALSE))

options(idlBNs.debug.pasets=old_opt)

cat("all pasets tests passed\n")
