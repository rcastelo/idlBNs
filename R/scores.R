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
#' @param x Either the data, or the population it would have come from. A
#' `data.frame` or `matrix` of Gaussian data, with observations in the rows
#' and random variables in the columns; or a population model built with
#' [`population`], in which case the score is evaluated in the large-sample
#' limit instead of on a sample. That is not a different score: it is the same
#' arithmetic fed the sufficient statistics a sample of the given size would
#' have in expectation. A bare `GaussParDAG` from the \pkg{pcalg} package is
#' accepted as a population model with a hard intervention to zero. See
#' `target.index` for how the notional sample size is supplied.
#'
#' @param targets (Default `list(integer(0))`) A `list` object with a family of
#' targets provided as a list of integer vectors. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param target.index (Default `NULL`) How much data comes from each
#' environment. What it holds, and what `NULL` resolves to, depend on whether
#' `x` carries data or a population.
#'
#' With data in `x`, a vector of integers in one-to-one correspondence with
#' the rows in `x`, saying which target intervened on each row. `NULL`
#' resolves to a vector of ones: the data is purely observational.
#'
#' With a population model in `x` there are no rows to label, so it is one
#' observation count per element of `targets` instead. Their sum is the
#' notional sample size, which sets the score's penalty term and so decides
#' between nested models in the limit. `NULL` resolves to [`population`]'s own
#' `n`, split in proportion to \eqn{(C, 1, \ldots, 1)} between the
#' observational and the interventional environments. If `n` was not given
#' there, omitting this argument is an error rather than a default: any size
#' invented on the caller's behalf would silently change which model is
#' selected.
#'
#' The counts need not be whole numbers, and a count of zero is allowed --
#' that environment contributes nothing and is dropped. What each variable
#' does need is at least two observations' worth of mass from the environments
#' that leave it alone, the same floor the row counts must clear when `x`
#' carries data; below it the variable cannot be scored and the call fails
#' naming it.
#'
#' An environment with no observations is removed from `targets` altogether,
#' under either kind of input -- a count of zero, or a target that no row
#' refers to. It is not merely that it cannot inform the score: the target
#' family is also what defines \emph{I}-equivalence and therefore which
#' reversals are \emph{I}-covered, so an intervention that was never
#' performed would otherwise narrow the equivalence classes the search moves
#' in and change the graph returned. The refinement is earned by having
#' observed the intervention.
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
#' on the input data (`x`), the target vertices (`targets`) and the target
#' indices (`target.index`) of the interventions. If `NULL` (default), the
#' `.iBIC.global.sufstats()` function is internally called.
#'
#' @param engine (Default `"C"`) A character string selecting the computation
#' engine: `"C"` (default) uses a compiled C routine for speed; `"R"` uses the
#' pure-R implementation and is provided for testing and verification.
#'
#' @param pasets (Default `NULL`) An optional list of parent sets, one per
#' vertex in `g` in the order given by `colnames(x)`, as internally built
#' by `iBIC()` from the structure of `g`. If `NULL` (default), it is
#' internally computed from `g`. Search algorithms that maintain `pasets`
#' incrementally across many calls (e.g. [hcmc()], [hillclimbing()]) can
#' pass it in directly to skip rebuilding it from `g` on every call.
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
#' x <- rbind(obsdat, intdat)
#'
#' ## define the targets and target indices for the interventional data
#' targets <- list(integer(0), 2L)
#' target.index <- c(rep(1L, nobs), rep(2L, nint))
#'
#' ## calculate the interventional BIC score for the DAG and data set
#' iBIC(g, x, targets, target.index)
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
#' iBIC(g2, x, targets, target.index)
#'
#' ## this is not the case if we do not indicate the presence of interventions
#' ## in the data
#' iBIC(g, x)
#' iBIC(g2, x)
#'
#' @importFrom graph edgeMatrix numNodes
#' @export
iBIC <- function(g, x, targets=list(integer(0)),
                 target.index=NULL,
                 cached.scores=NULL, global.sufstats=NULL,
                 engine=c("C", "R"), pasets=NULL) {

    engine <- match.arg(engine)

    if (is.null(attr(x, "sanitycheck"))) {
        x <- .check_input_data(x)
        .check_g_dat_consistency(g, x)
    }
    target.index <- .resolve.target.index(x, targets, target.index)

    if (is.null(pasets))
        pasets <- .build_pasets(g, x)
    else if (!is.list(pasets) || length(pasets) != numNodes(g) ||
            any(vapply(pasets, function(x) !is.integer(x), logical(1))))
            cli_abort(c("x"="'pasets' must be a list of integer vectors"))
    .check_cached_scores(g, cached.scores)

    if (is.null(global.sufstats))
        global.sufstats <- .iBIC.global.sufstats(x, targets, target.index)

    if (engine == "C")
        return(.Call(C_iBIC_score,
                     global.sufstats$S,
                     pasets,
                     as.double(global.sufstats$data.count),
                     as.double(global.sufstats$n),
                     cached.scores))

    sco <- numeric(length(pasets))
    for (i in seq_along(pasets)) {
        s <- NULL
        if (!is.null(cached.scores)) {
            k <- .cached_scores_key(pasets[[i]])
            s <- cached.scores[[i]][[k]]
        }
        if (is.null(s)) {
            Sj <- global.sufstats$S[[i]]
            idx <- c(1L, pasets[[i]] + 1L)
            ZtZ <- Sj[idx, idx, drop=FALSE]
            ZtY <- Sj[idx, i + 1L]
            YtY <- Sj[i + 1L, i + 1L]
            R <- chol(ZtZ)
            cc <- backsolve(R, ZtY, transpose=TRUE)
            RSS <- YtY - sum(cc^2)

            Nj <- global.sufstats$data.count[i]
            lambda <- 0.5 * log(global.sufstats$n)
            s <- -0.5 * Nj * (1 + log(RSS / Nj)) -
                 lambda * (1 + length(pasets[[i]]))

            if (!is.null(cached.scores))
                cached.scores[[i]][[k]] <- s
        }

        sco[i] <- s
    }

    sum(sco)
}
## assign a name attribute to the iBIC() scoring function for reporting purposes
attr(iBIC, "scorefun.name") <- "iBIC"

## build the list of parent sets of every vertex in g, one element per
## vertex in the order given by colnames(dat), where pasets[[i]] is an
## integer vector of the dat-column-index parents of vertex i. shared by
## iBIC() and iBGe(), and used as the "ground truth" to validate an
## incrementally-maintained pasets structure (see init.pasets() et al. in
## search.R) against.
#' @importFrom graph nodes edgeMatrix
.build_pasets <- function(g, dat) {
    v <- match(nodes(g), colnames(dat))
    em <- edgeMatrix(g)
    stopifnot(is.integer(em), nrow(em) == 2L, ncol(em) >= 0L) ## QC
    pasets <- split(em["from", ], factor(v[em["to", ]], levels=v))
    stopifnot(identical(names(pasets), as.character(v)))
    pasets
}

## mark iBIC() as able to accept a precomputed 'pasets' argument, so that
## hcmc()/hillclimbing() can pass in an incrementally-maintained pasets
## structure instead of having iBIC() rebuild it from g on every call. this
## is a deliberate opt-in attribute rather than formals() introspection: a
## custom user scorefun that happens to have an unrelated 'pasets' parameter
## will never accidentally be fed one, since the attribute is never set on
## it.
attr(iBIC, "supports.pasets") <- TRUE

## score a whole neighborhood of candidate moves against the current DAG in
## a single .Call(), returning one total iBIC score per candidate move.
## 'op'/'u'/'v' are the three parallel integer vectors of a neighborhood as
## returned by nr.nh()/ar.nh()/ncr.nh(), with 'u'/'v' already translated
## from nodes(dag) positions to the dat-column indices that 'pasets' is
## keyed by. Only the one or two vertices whose parent set a move changes
## get rescored, instead of all p of them (see src/nh_scores.c), so a
## neighborhood costs O(|NH|) node scores rather than O(p * |NH|).
.iBIC.nh.scores <- function(op, u, v, pasets, global.sufstats, cached.scores)
    .Call(C_iBIC_nh_scores,
          global.sufstats$S,
          pasets,
          as.double(global.sufstats$data.count),
          as.double(global.sufstats$n),
          cached.scores,
          op, u, v)

## expose the neighborhood scorer to the search algorithms the same way the
## global sufficient statistics function is exposed: as a deliberate opt-in
## attribute, so a custom user scorefun without one keeps being called once
## per candidate move (see score.nh() in search.R)
attr(iBIC, "nh.scores.fun") <- .iBIC.nh.scores

## find the best candidate move in a neighbourhood, returning only
## list(index=, total=, band=, worst=) rather than all |NH| scores. The
## winner is located through an error-bounded candidate band, so the O(p)
## exact summation is paid only for the provable handful of candidates that
## could still be the maximum, instead of for every one of the O(p^2) of
## them (see the band commentary in src/nh_scores.c). The move it picks, and
## the total it reports, are exactly what ranking every candidate exactly
## would give -- including which of several tied candidates wins.
##
## 'verify' scores every candidate exactly as well and checks the band
## against it; it is O(p * |NH|) and exists for tests/test_c_band.R and the
## idlBNs.debug.band option.
##
## 'stamp' is the DAG's per-vertex parent-set version vector. When it is
## supplied together with a compiled cache, additions -- about 99% of the
## candidates -- are served from a memo keyed on (u, v) and validated by one
## integer compare against stamp[v], which is what removes the hash lookup
## from the hot path. Without it the memo stays off and every candidate goes
## through the cache.
.iBIC.nh.argmax <- function(op, u, v, pasets, global.sufstats, cached.scores,
                            verify=FALSE, stamp=NULL)
    .Call(C_iBIC_nh_argmax,
          global.sufstats$S,
          pasets,
          as.double(global.sufstats$data.count),
          as.double(global.sufstats$n),
          cached.scores,
          op, u, v, stamp, verify)

attr(iBIC, "nh.argmax.fun") <- .iBIC.nh.argmax

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
## the target vertices and the target indices of the interventions.
##
## for each vertex j, precompute the (p+1) x (p+1) raw cross-product matrix
## S[[j]] = crossprod(cbind(1, X)), where X is the input data matrix restricted
## to the non-intervened rows/observations for vertex j, and the first column
## of S[[j]] corresponds to the intercept term. Any submatrix (j, parents) can
## then get its residual sum of squares from a small submatrix of S[[j]] via a
## Cholesky decomposition, without accessing the original input N-row data
## matrix again, analogously to the iBGe score's TN matrix.

#' @importFrom cli cli_abort
.iBIC.global.sufstats <- function(dat, targets=list(integer(0)),
                                  target.index=NULL) {
    ## resolved here too: this is reachable directly, and as the
    ## global.sufstats.fun attribute, so it cannot lean on its callers having
    ## done it. rep(1L, nrow(dat)) as a default is what a population object
    ## cannot satisfy -- it has no rows.
    target.index <- .resolve.target.index(dat, targets, target.index)
    ## a population model in place of data: the same statistics in the
    ## large-sample limit, built in R/population.R
    if (.is.population(dat))
        return(.iBIC.population.sufstats(dat, targets, target.index))
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

    S <- vector("list", p)
    S.full <- NULL ## cross-product matrix for all variables using all
                   ## rows/observations, computed once and shared by reference
                   ## among vertices that have no interventions in the data
    for (j in seq_len(p)) {
        Nj <- data.count[j]
        if (Nj < 2) {
            msg <- paste("Not enough observational input data in column number",
                         j, "(", Nj, "observed values)")
            cli_abort(c("x"=msg))
        }
        ## if there are no interventions in this vertex, compute once the full
        ## cross-product matrix for all variables using all rows/observations,
        ## and share it by reference among vertices that have no interventions
        if (!onlyobsdata && Nj < n) {
            Xj <- dat[non.int[[j]], , drop=FALSE]
            S[[j]] <- crossprod(cbind(1, Xj))
        } else {
            if (is.null(S.full))
                S.full <- crossprod(cbind(1, dat))
            S[[j]] <- S.full
        }
    }

    list(p=p, n=n, non.int=non.int, data.count=data.count, S=S)
}

## assign the iBIC global sufficient statistics function as an attribute to the
## iBIC() scoring function, so that any search algorithm taking iBIC() as an
## input argument, e.g. scorefun=iBIC, can precompute the corresponding global
## sufficient statistics before iteratively calling iBIC() during search
attr(iBIC, "global.sufstats.fun") <- .iBIC.global.sufstats

#' @importFrom graph numNodes
#' @importFrom cli cli_abort
.check_cached_scores <- function(g, cached.scores) {
    if (is.null(cached.scores))
        return(invisible(NULL))

    ## a compiled score cache, as the C search engine creates with
    ## C_sccache_new(). The C side checks its own external-pointer tag and
    ## that it was built for the right number of vertices, so there is
    ## nothing useful to verify here -- and nothing to validate it AS, since
    ## an external pointer is opaque from R. The documented
    ## list-of-environments form is still handled below, unchanged.
    if (typeof(cached.scores) == "externalptr")
        return(invisible(NULL))

    ## fast path: skip the full O(p) validation below on repeat calls with
    ## the same cached.scores object, e.g. once per neighbor from inside
    ## hcmc()'s or hillclimbing()'s search loop, where the same
    ## cached.scores is threaded through every call for the whole search.
    ## the marker is stored inside one of the per-node environments
    ## (a reference object in R), so it persists across calls even though
    ## the 'cached.scores' list argument itself is a fresh local binding
    ## on every call.
    if (is.list(cached.scores) &&
        length(cached.scores) == numNodes(g) &&
        length(cached.scores) > 0L &&
        is.environment(cached.scores[[1L]]) &&
        isTRUE(cached.scores[[1L]]$.validated))
        return(invisible(NULL))

    if (!is.list(cached.scores) || length(cached.scores) != numNodes(g)) {
        msg <- paste("'cached.scores' must be a list of length equal to",
                     "the number of vertices in g (", numNodes(g), ")")
        cli_abort(c("x"=msg))
    }
    if (any(vapply(cached.scores, function(x) !is.environment(x),
                   logical(1)))) {
        msg <- paste("Each element of 'cached.scores' must be",
                     "an environment")
        cli_abort(c("x"=msg))
    }

    assign(".validated", TRUE, envir=cached.scores[[1L]])
    invisible(NULL)
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
#' @param x Either the data, or the population it would have come from. A
#' `data.frame` or `matrix` of Gaussian data, with observations in the rows
#' and random variables in the columns; or a population model built with
#' [`population`], in which case the score is evaluated in the large-sample
#' limit instead of on a sample. That is not a different score: it is the same
#' arithmetic fed the sufficient statistics a sample of the given size would
#' have in expectation. A bare `GaussParDAG` from the \pkg{pcalg} package is
#' accepted as a population model with a hard intervention to zero. See
#' `target.index` for how the notional sample size is supplied.
#'
#' @param targets (Default `list(integer(0))`) A `list` object with a family of
#' targets provided as a list of integer vectors. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param target.index (Default `NULL`) How much data comes from each
#' environment. What it holds, and what `NULL` resolves to, depend on whether
#' `x` carries data or a population.
#'
#' With data in `x`, a vector of integers in one-to-one correspondence with
#' the rows in `x`, saying which target intervened on each row. `NULL`
#' resolves to a vector of ones: the data is purely observational.
#'
#' With a population model in `x` there are no rows to label, so it is one
#' observation count per element of `targets` instead. Their sum is the
#' notional sample size, which sets the score's penalty term and so decides
#' between nested models in the limit. `NULL` resolves to [`population`]'s own
#' `n`, split in proportion to \eqn{(C, 1, \ldots, 1)} between the
#' observational and the interventional environments. If `n` was not given
#' there, omitting this argument is an error rather than a default: any size
#' invented on the caller's behalf would silently change which model is
#' selected.
#'
#' The counts need not be whole numbers, and a count of zero is allowed --
#' that environment contributes nothing and is dropped. What each variable
#' does need is at least two observations' worth of mass from the environments
#' that leave it alone, the same floor the row counts must clear when `x`
#' carries data; below it the variable cannot be scored and the call fails
#' naming it.
#'
#' An environment with no observations is removed from `targets` altogether,
#' under either kind of input -- a count of zero, or a target that no row
#' refers to. It is not merely that it cannot inform the score: the target
#' family is also what defines \emph{I}-equivalence and therefore which
#' reversals are \emph{I}-covered, so an intervention that was never
#' performed would otherwise narrow the equivalence classes the search moves
#' in and change the graph returned. The refinement is earned by having
#' observed the intervention.
#'
#' @param cached.scores An optional list of environment objects, containing
#' cached scores per parent set for each vertex in `g`. If `NULL` (default),
#' no cached scores are used. Using this argument can speed up the calculation
#' of the score when the same parent sets are scored multiple times. To use
#' this argument, first create an empty environment object with
#' `csco <- replicate(numNodes(g), new.env(hash=TRUE, parent=emptyenv()), simplify=FALSE)`
#' and then pass it to this `cached.scores` parameter, i.e.,
#' `cached.scores=csco`.
#'
#' @param global.sufstats (Default `NULL`) An optional list of global sufficient
#' statistics for the iBGe score, as returned by the `.iBGe.global.sufstats()`
#' function, which do not depend on the structure of a specific DAG, but only
#' on the input data (`x`), the target vertices (`targets`) and the target
#' indices (`target.index`) of the interventions. If `NULL` (default), the
#' `.iBGe.global.sufstats()` function is internally called.
#'
#' @param pasets (Default `NULL`) An optional list of parent sets, one per
#' vertex in `g` in the order given by `colnames(x)`, as internally built
#' by `iBGe()` from the structure of `g`. If `NULL` (default), it is
#' internally computed from `g`. Search algorithms that maintain `pasets`
#' incrementally across many calls (e.g. [hcmc()], [hillclimbing()]) can
#' pass it in directly to skip rebuilding it from `g` on every call.
#'
#' @param engine (Default `"C"`) A character string selecting the computation
#' engine: `"C"` (default) uses a compiled C routine for speed; `"R"` uses the
#' pure-R implementation and is provided for testing and verification.
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
#' x <- rbind(obsdat, intdat)
#'
#' ## define the targets and target indices for the interventional data
#' targets <- list(integer(0), 2L)
#' target.index <- c(rep(1L, nobs), rep(2L, nint))
#'
#' ## calculate the interventional BGe score for the DAG and data set
#' iBGe(g, x, targets, target.index)
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
#' iBGe(g2, x, targets, target.index)
#'
#' ## this is not the case if we do not indicate the presence of interventions
#' ## in the data
#' iBGe(g, x)
#' iBGe(g2, x)
#'
#' @importFrom methods as
#' @importFrom graph numNodes edgeMatrix
#' @export
iBGe <- function(g, x, targets=list(integer(0)),
                 target.index=NULL,
                 cached.scores=NULL, global.sufstats=NULL,
                 pasets=NULL, engine=c("C", "R")) {

    engine <- match.arg(engine)

    if (is.null(attr(x, "sanitycheck"))) {
        x <- .check_input_data(x)
        .check_g_dat_consistency(g, x)
    }
    target.index <- .resolve.target.index(x, targets, target.index)

    if (is.null(pasets))
        pasets <- .build_pasets(g, x)
    else if (!is.list(pasets) || length(pasets) != numNodes(g) ||
            any(vapply(pasets, function(x) !is.integer(x), logical(1))))
            cli_abort(c("x"="'pasets' must be a list of integer vectors"))
    .check_cached_scores(g, cached.scores)

    if (is.null(global.sufstats))
        global.sufstats <- .iBGe.global.sufstats(x, targets, target.index)

    if (engine == "C")
        return(.Call(C_iBGe_score,
                     global.sufstats$TN,
                     pasets,
                     as.double(global.sufstats$awpN),
                     as.double(global.sufstats$p),
                     global.sufstats$scoreconstvec,
                     cached.scores))

    sco <- numeric(length(pasets))
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

## mark iBGe() as able to accept a precomputed 'pasets' argument -- see the
## identical attribute on iBIC() for the rationale.
attr(iBGe, "supports.pasets") <- TRUE

## score a whole neighborhood of candidate moves against the current DAG in
## a single .Call() -- see .iBIC.nh.scores() for the rationale
.iBGe.nh.scores <- function(op, u, v, pasets, global.sufstats, cached.scores)
    .Call(C_iBGe_nh_scores,
          global.sufstats$TN,
          pasets,
          as.double(global.sufstats$awpN),
          as.double(global.sufstats$p),
          global.sufstats$scoreconstvec,
          cached.scores,
          op, u, v)

attr(iBGe, "nh.scores.fun") <- .iBGe.nh.scores

## find the best candidate move in a neighbourhood -- see
## .iBIC.nh.argmax() for the rationale
.iBGe.nh.argmax <- function(op, u, v, pasets, global.sufstats, cached.scores,
                            verify=FALSE, stamp=NULL)
    .Call(C_iBGe_nh_argmax,
          global.sufstats$TN,
          pasets,
          as.double(global.sufstats$awpN),
          as.double(global.sufstats$p),
          global.sufstats$scoreconstvec,
          cached.scores,
          op, u, v, stamp, verify)

attr(iBGe, "nh.argmax.fun") <- .iBGe.nh.argmax

## calculate global sufficient statistics for the iBGe score, which do not
## depend on the structure of a specific DAG, but only on the input data,
## the target vertices and the target indices of the interventions. part of
## this code is adapted from the BGe parametrisation in Kuipers & Moffa (2025)
## and the BiDAG package, but stripped down to exclude BDe, BDecat, DBN, MDAG,
## and other stuff not exposed in the iBGe() function of this package
.iBGe.global.sufstats <- function(dat, targets=list(integer(0)),
                                  target.index=NULL) {
    ## resolved here too: this is reachable directly, and as the
    ## global.sufstats.fun attribute, so it cannot lean on its callers having
    ## done it. rep(1L, nrow(dat)) as a default is what a population object
    ## cannot satisfy -- it has no rows.
    target.index <- .resolve.target.index(dat, targets, target.index)
    ## a population model in place of data: the same statistics in the
    ## large-sample limit, built in R/population.R
    if (.is.population(dat))
        return(.iBGe.population.sufstats(dat, targets, target.index))
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
    T0scale <- am * (aw - p - 1) / (am + 1) # follows [GH2002, eqs. (19, 20)]
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
        if (!onlyobsdata && length(non.int[[j]]) < n)
            Xj <- dat[non.int[[j]], , drop=FALSE]
        Nj <- data.count[j]
        if (Nj < 2) {
            msg <- paste("Not enough observational input data in column number",
                         j, "(", Nj, "observed values)")
            cli_abort(c("x"=msg))
        }
        means <- colMeans(Xj)
        covmat <- cov(Xj) * (Nj - 1)
        TN[[j]] <- T0 + covmat + (am * Nj / (am + Nj)) * outer(means, means)
        awpN[j] <- aw + Nj
        constscorefact <- -(Nj / 2) * log(pi) + 0.5 * log(am / (am + Nj))
        awp <- aw - p + l
        scoreconstvec[[j]] <- constscorefact - lgamma(awp / 2) +
                              lgamma((awp + Nj) / 2) +
                              ((awp + l - 1) / 2) * log(T0scale) - l * logedgepf
    }

    list(p=p, n=n, non.int=non.int, data.count=data.count, aw=aw,
         T0scale=T0scale, TN=TN, awpN=awpN, scoreconstvec=scoreconstvec)
}

## assign the iBGe global sufficient statistics function as an attribute to the
## iBGe() scoring function, so that any search algorithm taking iBGe() as an
## input argument, e.g. scorefun=iBGe, can precompute the corresponding global
## sufficient statistics before iteratively calling iBGe() during search
attr(iBGe, "global.sufstats.fun") <- .iBGe.global.sufstats


##
## VENDORED CODE FROM THE iBIC SCORE BY HAUSER AND BÜHLMANN (2015)
##

## first original version of the iBIC() function, which calls the vendored R
## code of the iBIC score by Hauser and Bühlmann (2015) based and adapted from
## the pcalg pacakge at https://cran.r-project.org/package=pcalg this is
## included here to verify that further optimized versions of the iBIC()
## function produce the same results as the original version of the iBIC score
## by Hauser and Bühlmann

.vendored_iBIC <- function(g, dat, targets=list(integer(0)),
                           target.index=rep(1L, nrow(dat)),
                           cached.scores=NULL, global.sufstats=NULL) { # nocov start

    if (is.null(attr(dat, "sanitycheck"))) {
        dat <- .check_input_data(dat)
        .check_g_dat_consistency(g, dat)
    }

    v <- match(nodes(g), colnames(dat))
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
            s <- -0.5 * global.sufstats$data.count[i] *
                 (1 + log(sigma2 / global.sufstats$data.count[i])) -
                 lambda * (1 + length(pasets[[i]]))

            if (!is.null(cached.scores))
                cached.scores[[i]][[k]] <- s
        }

        sco[i] <- s
    }

    sum(sco)
} # nocov end

##
## VENDORED CODE FROM THE iBGe SCORE BY KUIPERS AND MOFFA (2025)
##

## first original version of the iBGe() function, which calls the vendored code
## of the iBGe score by Kuipers and Moffa (2025) based and adapted from the
## scripts provided at https://github.com/jackkuipers/iBGe and from the BiDAG
## package at https://cran.r-project.org/package=BiDAG this is included here
## to verify that further optimized versions of the iBGe() function produce the
## same results as the original version of the iBGe score by Kuipers and Moffa

#' @importFrom graph nodes numNodes edgeMatrix
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
    for (jj in 1:n) {
        nint_obs <- which(Tmat[, jj] == 0)
        if (length(nint_obs) < 2) {
            stop("Not enough observational data.")
        } else {
            nodeparams[[jj]] <- .scoreparameters(scoretype=usrpar$pctesttype,
                                                 data=initparam$data[nint_obs, ],
                                                 weightvector=initparam$weightvector[nint_obs],
                                                 bgepar=list(am=usrpar$am),
                                                 bdepar=list(chi=usrpar$chi,
                                                             edgepf=usrpar$edgepf),
                                                 bdecatpar=list(chi=usrpar$chi,
                                                                edgepf=usrpar$edgepf),
                                                 edgepmat=usrpar$edgepmat)
        }
    }
    initparam$nodeparams <- nodeparams

    initparam
} # nocov end

## here we have put only the BGe part

#' @importFrom stats cov cov.wt
.scoreparameters <- function(scoretype=c("bge","bde","bdecat","usr"), data,
                             bgepar=list(am=1, aw=NULL, edgepf=1),
                             bdepar=list(chi=0.5, edgepf=2),
                             bdecatpar=list(chi=0.5, edgepf=2),
                             dbnpar=list(samestruct=TRUE, slices=2, b=0,
                                         stationary=TRUE, rowids=NULL,
                                         datalist=NULL, learninit=TRUE),
                             usrpar=list(pctesttype=c("bge","bde","bdecat")),
                             mixedpar=list(nbin=0), MDAG=FALSE, DBN=FALSE,
                             weightvector=NULL, bgnodes=NULL, edgepmat=NULL,
                             nodeslabels=NULL) { # nocov start

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
    if(!is.null(bgnodes))
        initparam$mainnodes<-c(1:n)[-bgnodes]
    else
        initparam$mainnodes<-c(1:n)

    initparam$bgn<-bgn
    initparam$n<-n
    initparam$nsmall<-nsmall

    initparam$labels.short<-initparam$labels

    if (is.null(edgepmat)) {
        initparam$logedgepmat <- NULL
    } else {
        if (all(edgepmat>0))
            initparam$logedgepmat <- log(edgepmat)
        else
            stop("all entries of edgepmat matrix must be bigger than 0! 1 corresponds to no penalization")
    }

    if (scoretype == "bge") {

        if (is.null(bgepar$am)) {
            bgepar$am<-1
        }
        if (is.null(bgepar$aw)) {
            bgepar$aw<-n+bgepar$am+1
        }
        if (is.null(bgepar$edgepf)) {
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
            initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) +
                                          lgamma((awp+N)/2) +
                                          ((awp+j-1)/2)*log(T0scale) - j*log(initparam$pf)
        }

    } else if (scoretype == "usr") { ## usr
        if (is.null(usrpar$pctesttype)){usrpar$pctesttype <- "usr"}
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
    if (scorepar$DBN) {
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
