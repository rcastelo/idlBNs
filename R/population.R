## Population ("large-sample limit") sufficient statistics.
##
## The scores in scores.R read nothing from the data except per-vertex
## sufficient statistics: for iBIC the intercept-augmented cross-product
## matrix S[[j]] over the rows in which j was NOT intervened, together with
## that row count N_j and the total n; for iBGe the matrix TN[[j]], the
## posterior degrees of freedom awpN[j] and a constant vector, all of which
## are functions of the same per-vertex mean, covariance and count.
##
## So a population score is not a different score. It is the same arithmetic
## fed the statistics a sample of that size would have IN EXPECTATION: replace
## the empirical cross-product by N_j times the population second moment of
## the environments that do not intervene on j. Every consumer downstream --
## the C fast paths, the score cache, the addition memo, both search engines
## -- is untouched.
##
## What the caller has to supply beyond the model is the EXPERIMENT:
##
##   targets       which variables each environment intervenes on, as usual
##   target.index  in the population case, one COUNT per environment rather
##                 than one entry per row; the notional sample size is their
##                 sum, and N_j is the total over the environments that leave
##                 j alone. This is the same quantity target.index carries for
##                 real data -- how much data comes from each environment --
##                 stated directly because there are no rows to count.
##   ivent.value   what an intervened variable is set to, and
##   ivent.var     the variance around it: 0 is a hard intervention to a
##                 constant, positive a stochastic one.
##
## The last two matter even though a vertex is never scored in an environment
## that intervenes on it: a vertex j IS scored against its parents, and a
## parent intervened in that environment contributes the intervention's own
## distribution. Pooling environments with different means also adds a
## between-environment term to the covariance, which is why the VALUE and not
## only the variance changes the score.

## Per-environment moments of a pcalg::GaussParDAG under each target.
##
## pcalg's weight matrix is indexed [from, to], so the equation for j reads
## X_j = intercept_j + sum_i B[i, j] X_i + e_j, and intervening on I means
## cutting the incoming weights B[, I], fixing the intercepts of I to the
## intervention value and their error variances to its variance. Checked
## against GaussParDAG$cov.mat(target, ivent.var), which it reproduces
## exactly; cov.mat() is not used directly because it cannot give the MEAN.
.pop.env.moments <- function(x, targets, ivent.value = 0, ivent.var = 0) {
    p  <- x$node.count()
    B  <- x$weight.mat()
    om <- x$err.var()
    b0 <- x$intercept()
    lapply(targets, function(I) {
        Bk <- B; omk <- om; b0k <- b0
        if (length(I)) {
            Bk[, I] <- 0
            omk[I]  <- ivent.var
            b0k[I]  <- ivent.value
        }
        M <- solve(diag(p) - t(Bk))
        list(mu = as.vector(M %*% b0k), Sigma = M %*% diag(omk, p) %*% t(M))
    })
}

## Pool a set of environments into one distribution. A mixture's covariance is
## the mean of the covariances PLUS the spread of the means about the mixture
## mean, which is the term that makes the intervention value matter.
.pop.pool <- function(moments, w) {
    w  <- w / sum(w)
    mu <- Reduce(`+`, Map(function(e, wi) wi * e$mu, moments, w))
    S  <- Reduce(`+`, Map(function(e, wi) wi * (e$Sigma + tcrossprod(e$mu - mu)),
                          moments, w))
    list(mu = mu, Sigma = S)
}

## For every vertex, the pooled moments of the environments that do not
## intervene on it, and the count of observations they contribute.
.pop.pervertex <- function(moments, targets, counts, p) {
    n <- sum(counts)

    ## An environment with no observations carries no information, so it is
    ## dropped before anything is pooled. Left in, it stays in the pooling
    ## weights, and a vertex that only the zero-count environments leave alone
    ## divides 0 by 0: NaN moments, and then a score that is WRONG RATHER THAN
    ## REFUSED. iBGe returned NaN with no error at all and the search died
    ## later on "missing value where TRUE/FALSE needed"; iBIC failed inside
    ## the Cholesky with "dpotrf failed". Neither named the cause.
    live    <- which(counts > 0)
    moments <- moments[live]
    targets <- targets[live]
    counts  <- counts[live]

    out <- lapply(seq_len(p), function(j) {
        keep <- which(!vapply(targets, function(I) j %in% I, TRUE))
        if (!length(keep))
            cli_abort(c("x"=paste("Every environment with observations",
                                  "intervenes on variable {j}, so it has none",
                                  "of its own to be scored on."),
                        "i"=paste("Give an environment that leaves it alone a",
                                  "non-zero count in 'target.index'.")))
        Nj <- sum(counts[keep])
        ## the floor the finite-sample path applies to its row counts, applied
        ## here to mass: below it the score is not defined, and a fractional
        ## count used to come back as a finite, plausible-looking number
        if (Nj < 2)
            cli_abort(c("x"=paste("The environments that leave variable {j}",
                                  "alone total {format(Nj)} observations,",
                                  "and two are needed to score it."),
                        "i"="Raise 'n', or the counts in 'target.index'."))
        pooled <- .pop.pool(moments[keep], counts[keep])
        c(pooled, list(N = Nj))
    })
    list(n = n, v = out)
}

## Validate the population form of target.index and return it as counts.
.pop.counts <- function(target.index, targets) {
    if (!is.numeric(target.index) || length(target.index) != length(targets))
        cli_abort(c("x"=paste("With a population model in 'x', 'target.index'",
                              "must be a numeric vector with one observation",
                              "count per element of 'targets'.")))
    if (any(!is.finite(target.index)) || any(target.index < 0) ||
        sum(target.index) <= 0) {
        msg <- paste("'target.index' counts must be finite, non-negative and",
                     "not all zero.")
        cli_abort(c("x"=msg))
    }
    as.numeric(target.index)
}

## Is this second argument a population model rather than data?
.is.population <- function(x)
    inherits(x, "idlBNsPopulation") || inherits(x, "GaussParDAG")

.pop.nodes <- function(x) {
    nm <- x$.nodes
    if (is.null(nm)) as.character(seq_len(x$node.count())) else as.character(nm)
}

#' Create a population parameter model for interventional Bayesian networks
#'
#' Wraps a generative model with the intervention semantics, so that it can be
#' passed where a data matrix normally goes and the scores are evaluated in the
#' large-sample limit instead of on a sample. See [`iBIC`] for what the
#' accompanying `targets` and `target.index` then mean.
#'
#' @param x A `GaussParDAG` object from the \pkg{pcalg} package: the
#' generative model.
#' @param n (Default `NULL`) The notional sample size the population stands
#' in for. It is not cosmetic: the score's penalty is a function of it, so it
#' is what decides between nested models in the limit, and there is no
#' defensible default. Setting it lets `target.index` be omitted, in which
#' case the size is split equally across the environments in `targets`; pass
#' `target.index` instead to allocate it unequally.
#'
#' @param C (Default 1) How much data an observational environment carries
#' relative to each interventional one, as in Wang, Solus, Yang and Uhler
#' (2017). With `k` interventional environments the notional size `n` is split
#' in proportion to \eqn{(C, 1, \ldots, 1)}, so the observational stratum gets
#' \eqn{C/(C+k)} of it and each interventional stratum \eqn{1/(C+k)}. An
#' environment counts as observational when its target is empty, so the split
#' does not depend on the order of `targets`, and `C` is inert when every
#' environment is of one kind. The default 1 is an equal split.
#'
#' This is the quantity their counterexample turns on: the graph they exhibit
#' is a population local maximum for 2-14\% of edge-weight draws at `C = 1` and
#' for 68-96\% at `C = 50`. It is used only when `target.index` is omitted;
#' passing counts states the allocation outright.
#'
#' @param ivent.value The value an intervened variable is set to. It affects
#' the score of the variables that are \emph{not} intervened, because pooling
#' environments whose means differ widens the covariance their scores are
#' computed from.
#' @param ivent.var The variance around that value. Zero, the default, is a
#' hard intervention to a constant; a positive value is a stochastic one.
#'
#' @return An object of class `idlBNsPopulation`. It reports `ncol()` and
#' `colnames()` like the data matrix it stands in for, and `nrow()` as `NA`,
#' there being no rows to count.
#'
#' @export
population <- function(x, n = NULL, C = 1, ivent.value = 0, ivent.var = 0) {
    if (!inherits(x, "GaussParDAG"))
        cli_abort(c("x"="'x' must be a pcalg::GaussParDAG object"))
    if (!is.null(n) && (!is.numeric(n) || length(n) != 1L || !is.finite(n)
                        || n <= 0))
        cli_abort(c("x"="'n' must be NULL or a positive finite numeric scalar"))
    if (!is.numeric(C) || length(C) != 1L || !is.finite(C) || C <= 0)
        cli_abort(c("x"="'C' must be a positive finite numeric scalar"))
    if (!is.finite(ivent.value) || !is.numeric(ivent.value) ||
        length(ivent.value) != 1L || !is.finite(ivent.var) ||
        !is.numeric(ivent.var) || length(ivent.var) != 1L || ivent.var < 0)
        cli_abort(c("x"=paste("'ivent.value' must be a numeric scalar and",
                              "'ivent.var' a non-negative numeric scalar")))
    structure(list(model = x, n = n, C = C,
                   ivent.value = ivent.value, ivent.var = ivent.var,
                   nodes = .pop.nodes(x), p = x$node.count()),
              class = "idlBNsPopulation")
}

.as.population <- function(x)
    if (inherits(x, "idlBNsPopulation")) x else population(x)

## So that everything reading the second argument -- .build_pasets(),
## .check_g_dat_consistency(), and hcmc()'s node names and p -- works on a
## population object with no change at all. nrow() is NA because there are no
## rows: the sample size lives in target.index, deliberately.
#' @export
dim.idlBNsPopulation <- function(x) c(NA_integer_, x$p)
#' @export
dimnames.idlBNsPopulation <- function(x) list(NULL, x$nodes)

## iBIC's statistics from a population model. S[[j]] is the cross-product a
## sample of N_j observations would have in expectation: N_j * E[(1,X)'(1,X)]
## over the pooled environments that leave j alone.
.iBIC.population.sufstats <- function(x, targets = list(integer(0)),
                                      target.index = NULL) {
    x  <- .as.population(x); p <- x$p
    cn <- .pop.counts(target.index, targets)
    pv <- .pop.pervertex(.pop.env.moments(x$model, targets, x$ivent.value,
                                          x$ivent.var), targets, cn, p)
    S <- lapply(pv$v, function(e)
        e$N * rbind(c(1, e$mu),
                    cbind(e$mu, e$Sigma + tcrossprod(e$mu))))
    list(p = p, n = pv$n, non.int = NULL,
         data.count = vapply(pv$v, function(e) e$N, 0), S = S)
}

## iBGe's statistics from the same per-vertex moments. The empirical
## cov(Xj) * (Nj - 1) becomes (N_j - 1) * Sigma_j and colMeans(Xj) becomes
## mu_j; everything else is a function of N_j, p and the fixed hyperparameters,
## and is computed exactly as .iBGe.global.sufstats() does.
.iBGe.population.sufstats <- function(x, targets = list(integer(0)),
                                      target.index = NULL) {
    x  <- .as.population(x); p <- x$p
    cn <- .pop.counts(target.index, targets)
    pv <- .pop.pervertex(.pop.env.moments(x$model, targets, x$ivent.value,
                                          x$ivent.var), targets, cn, p)

    am <- 1                       # as in .iBGe.global.sufstats()
    edgepf <- 1
    aw <- p + am + 1
    T0scale <- am * (aw - p - 1) / (am + 1)
    T0 <- diag(T0scale, p, p)
    logedgepf <- log(edgepf)
    l <- seq_len(p)

    TN <- vector("list", p); awpN <- numeric(p)
    scoreconstvec <- vector("list", p)
    for (j in seq_len(p)) {
        Nj <- pv$v[[j]]$N; means <- pv$v[[j]]$mu
        covmat <- pv$v[[j]]$Sigma * (Nj - 1)
        TN[[j]] <- T0 + covmat + (am * Nj / (am + Nj)) * outer(means, means)
        awpN[j] <- aw + Nj
        constscorefact <- -(Nj / 2) * log(pi) + 0.5 * log(am / (am + Nj))
        awp <- aw - p + l
        scoreconstvec[[j]] <- constscorefact - lgamma(awp / 2) +
                              lgamma((awp + Nj) / 2) +
                              ((awp + l - 1) / 2) * log(T0scale) - l * logedgepf
    }
    list(p = p, n = pv$n, non.int = NULL,
         data.count = vapply(pv$v, function(e) e$N, 0), aw = aw,
         T0scale = T0scale, TN = TN, awpN = awpN,
         scoreconstvec = scoreconstvec)
}


## Resolve an omitted target.index.
##
## The formal default used to be rep(1L, nrow(x)), which a population object
## cannot satisfy: it has no rows, nrow() is NA by construction, and the
## default blew up inside rep() with "invalid 'times' argument" before any
## check could say something useful.
##
## For data the behaviour is unchanged. For a population model there is no
## defensible default, because the notional sample size sets the score's
## penalty and therefore decides between nested models -- inventing one would
## silently change the answer. So it has to be stated, either as n on the
## population object (split equally here) or as explicit counts.
.resolve.target.index <- function(x, targets, target.index) {
    if (!is.null(target.index))
        return(target.index)
    if (.is.population(x)) {
        px <- .as.population(x)
        if (is.null(px$n))
            cli_abort(c("x"=paste("With a population model in 'x' the notional",
                                  "sample size has to be given, because the",
                                  "score's penalty depends on it."),
                        "i"=paste("Either set 'n' in population(), or pass one",
                                  "observation count per element of 'targets'",
                                  "in 'target.index'.")))
        ## split n in proportion to (C, 1, ..., 1): an environment is
        ## observational when its target is empty, so this does not depend on
        ## the order of targets, and C is inert when they are all one kind
        w <- ifelse(vapply(targets, function(I) length(I) == 0L, TRUE), px$C, 1)
        return(px$n * w / sum(w))
    }
    rep(1L, nrow(x))
}
