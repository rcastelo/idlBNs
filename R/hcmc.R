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
#' @param MAXTRIALS (Default 5) Maximum number of trials to escape from local
#' maxima.
#'
#' @param scorefun (Default is [`iBIC`]) A function to calculate the goodness
#' of fit (GoF) score of a DAG on a given data set.
#'
#' @param verbose (Default TRUE) Show progress in the calculations.
#'
#' @param engine (Default `"C"`) A character string selecting the search
#' engine: `"C"` (default) maintains the DAG, its ancestor relation and its
#' parent sets in compiled code, and performs the (\emph{I}-)covered arc
#' reversals there too; `"R"` uses the pure-R implementation and is provided
#' for testing and verification. Both consume the random number stream
#' identically and so follow the same trajectory from the same seed, and both
#' return the same result. `"C"` requires a `scorefun` able to score a whole
#' neighbourhood at once, which [`iBIC`] and [`iBGe`] are; with any other
#' score function the `"R"` engine is used regardless.
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
#'     library(graph)
#'     library(pcalg)
#'     library(idlBNs)
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
#'     targets <- I[[v]]
#'     dat[[v]] <- rmvnorm.ivent(nbytgts[v], Mg, target=targets,
#'                               target.value=rep(2, length(targets)))
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
#' @importFrom stats setNames
#' @export
#' @rdname hcmc
hcmc <- function(dat, r=20, targets=list(integer(0)),
                 target.index=rep(1L, nrow(dat)),
                 scorefun=iBIC, MAXTRIALS=5, verbose=TRUE,
                 engine=c("C", "R")) {

    engine <- match.arg(engine)
    dat <- .check_input_data(dat)
    dag <- graphNEL(colnames(dat), edgemode="directed")
    attr(dat, "sanitycheck") <- TRUE

    stopifnot(is.list(targets)) ## QC
    scorefun <- match.fun(scorefun)

    ## the attributes that decide which engine can run, extracted before the
    ## score cache is built because they decide which KIND of cache it is
    scorefun.name <- attr(scorefun, "scorefun.name")
    supports.pasets <- isTRUE(attr(scorefun, "supports.pasets"))
    nh.scores.fun <- attr(scorefun, "nh.scores.fun")
    nh.argmax.fun <- attr(scorefun, "nh.argmax.fun")

    ## the C engine keeps the DAG, its ancestor relation and its parent sets
    ## in compiled state, and enumerates the neighbourhood there. It needs a
    ## scorefun that can score a whole neighbourhood in one call, because it
    ## holds no graphNEL during the search and so cannot serve score.nh()'s
    ## per-candidate fallback. iBIC() and iBGe() qualify; anything else falls
    ## back to R rather than failing.
    use.c <- engine == "C" && !is.null(nh.argmax.fun)
    ## idlBNs.debug.band makes the C engine additionally score every
    ## candidate exactly and check the error-bounded band against it, so the
    ## band's argmax, its reported total and its bound are all verified at
    ## every step. O(p * |NH|) per step, i.e. it gives back exactly what the
    ## band saves, so it is for testing only.
    verify.band <- isTRUE(getOption("idlBNs.debug.band", FALSE))

    ## the score cache. The C engine uses a compiled open-addressing table
    ## keyed on the parent set's integers; the R engine uses the documented
    ## list of per-vertex environments keyed on the same set as a string.
    ## Both key identically and neither evicts, so they hold the same
    ## entries with the same values after the same search -- which matters,
    ## because a node score depends on the parent SEQUENCE while the key is
    ## the parent SET, so the cache is part of the arithmetic. See
    ## src/sccache.h and tests/test_c_cache.R.
    if (use.c)
        cached.scores <- .Call(C_sccache_new, ncol(dat))
    else {
        if (!.load_suggested_package("RBGL")) {
            msg <- paste("The R engine requires the Bioconductor package",
                         "RBGL and it cannot be loaded.")
            cli_abort(c=("x"=msg))
        }
        cached.scores <- list()
        for (i in seq_len(ncol(dat)))
            cached.scores[[i]] <- new.env(hash=TRUE, parent=emptyenv())
    }

    global.sufstats <- NULL
    global.sufstats.fun <- attr(scorefun, "global.sufstats.fun")
    if (!is.null(global.sufstats.fun)) {
        if (verbose)
          cli_alert_info("Calculating global sufficient statistics")
        global.sufstats <- global.sufstats.fun(dat, targets, target.index)
    }

    anc <- init.ancestors(colnames(dat))
    vidx <- setNames(seq_len(ncol(dat)), colnames(dat))
    ## nodes(dag) never changes during the search, only its edges do, so the
    ## vertex names and the nodes(dag)-position -> dat-column map that
    ## translate a move's integer u/v are built once here
    vnames <- nodes(dag)
    vidx.nodes <- unname(vidx[vnames])
    pasets <- init.pasets(ncol(dat))
    utargets <- sort(unique(unlist(targets)))
    ## the C entry points take an INTSXP, and utargets comes out numeric if
    ## the caller wrote targets=list(integer(0), 2) rather than 2L
    utargets.i <- as.integer(utargets)
    ## length(0:r) is computed here, not as r + 1 in C: 0:r for a
    ## non-integer r is 0:floor(r) and for a negative r counts down, so the
    ## coercion stays where it already behaves correctly (see src/rcar.c)
    rlen <- length(0:r)
    s0 <- -Inf
    s1 <- scorefun(g=dag, dat=dat, targets=targets, target.index=target.index,
                   cached.scores=cached.scores, global.sufstats=global.sufstats)
    was_in_local_maximum <- local_maximum <- s1 < s0
    trials <- escapes <- avg_trials_per_escape <- 0

    if (verbose) {
      algname <- if (identical(targets, list(integer(0)))) "HCMC" else "iHCMC"
      msg <- "Running the {algname} algorithm"
      if (!is.null(scorefun.name))
          msg <- paste(msg, "with the {scorefun.name} score function")
      cli_progress_bar(msg)
      msg <- "Score {s1} Escapes {escapes} Trials {avg_trials_per_escape}"
      cli_progress_step(msg, spinner=TRUE)
    }

    ## the C engine keeps the DAG, its ancestor relation and its parent sets
    ## in compiled state behind an external pointer, enumerates the NCR
    ## neighbourhood there, and performs the covered arc reversals there --
    ## reproducing R's random stream draw for draw, so both engines follow
    ## the same trajectory from the same seed. The R engine is the reference
    ## implementation the differential tests score it against.
    ##
    ## The escape bookkeeping (trials, escapes, avg_trials_per_escape) is
    ## identical in both and stays in R.
    if (use.c) {
        st <- .Call(C_dag_new, ncol(dat))
        while (!local_maximum) {
            s0 <- s1
            .Call(C_dag_rcar, st, rlen, utargets.i)
            ne <- .Call(C_dag_nh, st, 3L, utargets.i)    ## 3 = ncr
            pasets <- .Call(C_dag_pasets, st)
            ## only the winner is needed, so the O(p) exact summation is
            ## paid for the provable handful of candidates that could still
            ## be the maximum rather than for all O(p^2) of them
            am <- nh.argmax.fun(ne$op, vidx.nodes[ne$u], vidx.nodes[ne$v],
                                pasets, global.sufstats, cached.scores,
                                verify.band, .Call(C_dag_pastamp, st))
            b <- am$index
            s1 <- am$total
            local_maximum <- s1 <= s0
            if (!local_maximum) {
                .Call(C_dag_apply_move, st, ne$op[b], ne$u[b], ne$v[b])
                if (was_in_local_maximum) {
                    escapes <- escapes + 1
                    avg_trials_per_escape <- (avg_trials_per_escape *
                                              (escapes-1) + trials) / escapes
                    was_in_local_maximum <- FALSE
                }
                trials <- 0
            } else if (trials < MAXTRIALS) {
                s1 <- s0
                .Call(C_dag_rcar, st, rlen, utargets.i)
                local_maximum <- FALSE
                was_in_local_maximum <- TRUE
                trials <- trials + 1
            } else
                s1 <- s0

            .debug_assertions(st, NULL, NULL, NULL, dat, vnames)

            if (verbose)
                cli_progress_update()
        }
        dag <- .graphNEL_from_edgeM(vnames, .Call(C_dag_edgeM, st))
    } else {
        while (!local_maximum) {
            s0 <- s1
            rcar.out <- rcar(dag, r, utargets, anc, pasets, vidx)
            dag <- rcar.out$dag
            anc <- rcar.out$anc
            pasets <- rcar.out$pasets
            ne <- ncr.nh(dag, anc, utargets)
            sco <- score.nh(ne, dag, dat, targets, target.index, cached.scores,
                            global.sufstats, pasets, vidx.nodes,
                            supports.pasets, scorefun, nh.scores.fun)
            b <- which.max(sco)
            b.op <- ne$op[b]
            b.u <- ne$u[b]
            b.v <- ne$v[b]
            s1 <- sco[b]
            local_maximum <- s1 <= s0
            if (!local_maximum) {
                ## 'anc' and 'pasets' are updated before 'dag', since
                ## remove.ancestors()/reverse.ancestors() read the DAG as it
                ## stood BEFORE the move
                anc <- switch(b.op,
                              add.ancestors(anc, vnames[b.u], vnames[b.v]),
                              remove.ancestors(anc, dag, vnames[b.u],
                                               vnames[b.v]),
                              reverse.ancestors(anc, dag, vnames[b.u],
                                                vnames[b.v]))
                pasets <- move.pasets(pasets, b.op, vidx.nodes[b.u],
                                      vidx.nodes[b.v])
                dag <- apply.move(dag, b.op, b.u, b.v, vnames)
                if (was_in_local_maximum) {
                    escapes <- escapes + 1
                    avg_trials_per_escape <- (avg_trials_per_escape *
                                              (escapes-1) + trials) / escapes
                    was_in_local_maximum <- FALSE
                }
                trials <- 0
            } else if (trials < MAXTRIALS) {
                s1 <- s0
                rcar.out <- rcar(dag, r, utargets, anc, pasets, vidx)
                dag <- rcar.out$dag
                anc <- rcar.out$anc
                pasets <- rcar.out$pasets
                local_maximum <- FALSE
                was_in_local_maximum <- TRUE
                trials <- trials + 1
            } else
                s1 <- s0

            .debug_assertions(NULL, dag, anc, pasets, dat, vnames)

            if (verbose)
                cli_progress_update()
        }
    }

    if (verbose) {
        algname <- if (identical(targets, list(integer(0)))) "HCMC" else "iHCMC"
        cli_progress_done("{algname} algorithm completed")
    }

    list(dag=dag, sco=s1)
}

#' @importFrom cli cli_abort
.check_input_data <- function(dat) {
    if (!is.data.frame(dat) && !is.matrix(dat)) {
        msg <- paste("Input data in 'dat' must be a data.frame or a",
                     "matrix object.")
        cli_abort(c("x"=msg))
    }

    if (is.null(colnames(dat))) {
        msg <- paste("Input data in 'dat' must have column names",
                     "corresponding to the random variables of the sought DAG.")
        cli_abort(c("x"=msg))
    }

    dat <- as.matrix(dat)
    if (!is.numeric(dat))
        cli_abort(c("x"="Input data in 'dat' must be numeric."))

    if (nrow(dat) < 4)
        cli_abort(c("x"="Input data in 'dat' must have 3 or more rows."))

    if (ncol(dat) < 2)
        cli_abort(c("x"="Input data in 'dat' must have 2 or more columns."))

    dat
}
