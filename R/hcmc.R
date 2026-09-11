#' The HCMC and iHCMC algorithms
#'
#' @description Run the hill-climber Monte Carlo (HCMC) algorithm (Castelo and
#' Kočka, 2003) on purely observational Gaussian data, or the interventional
#' HCMC (iHCMC) on mixed observational and interventional Gaussian data
#' (Castelo, 2026).
#'
#' @param x Either the data to learn from, or the population it would have
#' come from. A `data.frame` or `matrix` of Gaussian data, with observations
#' in the rows and random variables in the columns; or a population model
#' built with [`population`], in which case the search is scored in the
#' large-sample limit instead of on a sample, which is what lets a run be read
#' as the behaviour of the algorithm itself rather than of one dataset. A bare
#' `GaussParDAG` from the \pkg{pcalg} package is accepted as a population
#' model with a hard intervention to zero. See `target.index` for how the
#' notional sample size is supplied.
#'
#' @param r (Default 20) Non-negative integer scalar indicating the maximum
#' number of (\emph{I}-)covered arc reversals by the RCAR algorithm (Castelo
#' and Kočka, 2003).
#'
#' @param targets (Default `list(integer(0))`) A `list` object with a family of
#' targets provided as a list of integer vectors. Its default value indicates
#' that there are no interventions in the data, i.e., the data is purely
#' observational.
#'
#' @param target.index (Default a unit vector) How much data comes from each
#' environment. With data in `x`, a vector of integers in one-to-one
#' correspondence with the rows in `x`, indicating which rows in the input
#' data are intervened by which targets; its default value indicates that
#' there are no interventions in the data, i.e., the data is purely
#' observational. With a population model in `x` there are no rows to label,
#' so it is instead one observation count per element of `targets`, and the
#' notional sample size -- which is what sets the score's penalty term -- is
#' their sum.
#'
#' @param MAXTRIALS (Default 5) Non-negative integer scalar indicating the
#' maximum number of trials to escape from local maxima when `escape="trials"`.
#' It is ignored when `escape="exhaustive"`, unless the exhaustive enumeration
#' of the (\emph{I}-)equivalence class exceeds `escape.max` or a chain component
#' has more than 64 vertices, in which case the escape mechanism falls back to
#' the RCAR algorithm for a maximum of `MAXTRIALS` trials. See the `escape`
#' argument for details.
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
#' @param sampler (Default `"rcar"`) A character string selecting how the
#' algorithm moves within the (\emph{I}-)equivalence class of the current DAG.
#' `"rcar"` performs the RCAR algorithm (a random walk of up to `r`
#' (\emph{I}-)covered arc reversals) of Castelo and Kočka (2003). `"exact"`
#' instead draws a member of the class uniformly at random, by the
#' Clique-Picking algorithm of Wienöbst \emph{et al.} (2023), and ignores `r`.
#' The walk produced by the RCAR algorithm and the exact draw are not
#' equivalent: the walk is a random walk on the class, so its equilibrium is
#' proportional to the number of (\emph{I}-)covered arcs of each member and is
#' not uniform for any value of `r`.  The exact draw has no limit on the size
#' of the class, and none on the number of vertices in the DAG. It declines,
#' reverting to the walk and counting the draw in `sampler.fallbacks`, in two
#' cases. First, vertex sets within an undirected chain component are 64-bit
#' masks, so a component of more than 64 vertices is out of reach; components
#' that large need a DAG with almost no immoralities, and the largest seen for
#' random DAGs up to \eqn{p = 500} is 16. Second, uniformity needs the
#' per-component counts to be \emph{exact}, because the draw compares a
#' uniform integer against cumulative counts: doubles hold integers exactly
#' only to \eqn{2^{53}}, and a clique of more than 18 vertices has a
#' factorial past that, so beyond either bound the draw is declined rather
#' than made with rounded weights. Counting is unaffected and still answers
#' for a class of any size.
#'
#' @param escape (Default `"trials"`) A character string selecting what happens
#' at a local maximum. `"trials"` re-randomises the current DAG within its
#' class and retries, up to `MAXTRIALS` times. `"exhaustive"` instead examines
#' \emph{every} member of the class and takes the best move available from any
#' of them, which settles the question of whether the search is really at a
#' local maximum of the class; `MAXTRIALS` is then unused. It falls back to
#' `"trials"` in the cases where the enumeration declines, counted in
#' `escape.fallbacks`: a class larger than `escape.max` or `INT_MAX`, and a
#' class whose size cannot be computed at all because some chain component has
#' more than 64 vertices. The second is the 64-bit mask limit described under
#' `sampler`, and applies here whichever `sampler` is in use, since the escape
#' enumerates the class regardless of how the search moves within it.
#'
#' @param escape.max (Default 512) Positive integer scalar giving the largest
#' (\emph{I}-)equivalence class the `escape="exhaustive"` enumeration will
#' walk. That escape scores one whole neighbourhood per member, so its cost
#' grows linearly in the size of the class, which is why it is the class size
#' that is bounded. Producing the members is not itself the expensive part:
#' each one is generated directly from its index in the class, so listing
#' costs one step per member listed and no more. What `escape.max` bounds is
#' therefore the scoring, not the enumeration. `NA` or `Inf` mean no limit,
#' but the default of 512 is a safe value for the largest classes seen in
#' random DAGs up to \eqn{p = 500} vertices.
#'
#' @return A list with the following components: `dag`, a
#' [`graphNEL`][graph::graphNEL-class] object with the structure of the learned
#' DAG; `sco`, its score; `sampler.fallbacks`, the number of draws for which
#' `sampler="exact"` declined and the `"rcar"` walk was used instead; and
#' `escape.fallbacks`, the number of local maxima at which
#' `escape="exhaustive"` declined and the `MAXTRIALS` budget was used instead.
#' Both counts are zero unless the corresponding exact method was selected, and
#' a non-zero count means part of the run silently used the older machinery:
#' the walk does not sample the (\emph{I}-)equivalence class uniformly, so a
#' result obtained with `sampler.fallbacks > 0` is not the one
#' `sampler="exact"` promises.
#' 
#' @references Castelo, R. and Kočka, T. On inclusion-driven learning of
#' Bayesian networks. *Journal of Machine Learning Research*, 4:527-574, 2003.
#'
#' @references Castelo, R. Interventional idlBNs in DAG-space. In *Challenges
#' and Algorithms for Knowledge Discovery from Data*, M. van Leeuwen and
#' J.  Vreeken (eds.). LNCS 16067, Festschrift, Springer, 2026.
#'
#' @references Wienöbst, M., Bannach, M. and Liśkiewicz, M. Polynomial-time
#' algorithms for counting and sampling Markov equivalent DAGs with
#' applications. *Journal of Machine Learning Research*, 24(213):1-45, 2023.
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
#' x <- list()
#' for (v in seq_along(I)) {
#'     targets <- I[[v]]
#'     x[[v]] <- rmvnorm.ivent(nbytgts[v], Mg, target=targets,
#'                               target.value=rep(2, length(targets)))
#' }
#' x <- do.call("rbind", x)
#'
#' ## store the target index for each row of the data
#' tindex <- rep(1:length(nbytgts), nbytgts)
#'
#' ## run the HCMC algorithm assuming all data were observational
#' dhat.hcmc <- hcmc(x)
#' dhat.hcmc
#'
#' ## calculate the structural Hamming distance (SHD) between the generative
#' ## DAG and the estimated DAG
#' shd(e, dag2essgraph(dhat.hcmc$dag))
#'
#' ## run the iHCMC algorithm informing the presence of interventional data
#' ## using by the default the interventional BIC score (see the iBIC()
#' ## function).
#' dhat.ihcmc <- hcmc(x, targets=I, target.index=tindex)
#' dhat.ihcmc
#'
#' ## the estimated DAG is closer to the generative DAG (lower SHD value)
#' ## than the one estimated by HCMC, which did not take into account the
#' ## presence of interventional data
#' shd(e, dag2essgraph(dhat.ihcmc$dag))
#'
#' ## run it again this time using the interventional BGe score (see the
#' ## iBGe() function), which provides an estimate closer to the generative DAG
#' dhat.ihcmc2 <- hcmc(x, targets=I, target.index=tindex, scorefun=iBGe)
#' shd(e, dag2essgraph(dhat.ihcmc2$dag))
#'
#' @importFrom graph nodes edgeMatrix graphNEL
#' @importClassesFrom graph graphNEL
#' @importFrom cli cli_alert_info cli_progress_bar cli_progress_done
#' @importFrom cli cli_alert_success cli_alert_warning
#' @importFrom stats setNames
#' @export
#' @rdname hcmc
hcmc <- function(x, r=20, targets=list(integer(0)),
                 target.index=rep(1L, nrow(x)),
                 scorefun=iBIC, MAXTRIALS=5, verbose=TRUE,
                 engine=c("C", "R"),
                 sampler=c("rcar", "exact"), escape=c("trials", "exhaustive"),
                 escape.max=512) {

    engine <- match.arg(engine)
    sampler <- match.arg(sampler)
    escape <- match.arg(escape)
    ## escape.max caps the CLASS SIZE the exhaustive escape will walk. That
    ## escape scores one whole neighbourhood per member, which at p = 200 is
    ## about 10 ms, so its cost is linear in |[D]_I| with a large constant and
    ## has to be bounded by the class size. The size is known exactly and
    ## cheaply before anything is enumerated, so a class over budget costs
    ## nothing to refuse; above it the search falls back to the MAXTRIALS
    ## budget. Sampling needs no limit of this kind at all.
    ##
    ## Listing the members is now output-proportional -- each is decoded
    ## directly from its index in the class, one step per member -- rather than
    ## Theta(m!), so the bound is purely about the scoring it feeds and can sit
    ## far higher than the enumeration once allowed.
    x <- .check_input_data(x)
    dag <- graphNEL(colnames(x), edgemode="directed")
    attr(x, "sanitycheck") <- TRUE

    targets <- .check_targets(targets, ncol(x))
    escape.maxD <- .check_escape.max(escape.max)
    scorefun <- match.fun(scorefun)

    ## the r argument must be a finite non-negative integer scalar
    if (!is.numeric(r) || length(r) != 1 || is.na(r) || r < 0 ||
        r != floor(r)) {
        msg <- paste("The 'r' argument must be a finite non-negative",
                     "integer scalar.")
        cli_abort(c("x"=msg))
    }

    ## the MAXTRIALS argument must be a finite non-negative integer scalar
    if (!is.numeric(MAXTRIALS) || length(MAXTRIALS) != 1 ||
        is.na(MAXTRIALS) || MAXTRIALS < 0 || MAXTRIALS != floor(MAXTRIALS)) {
        msg <- paste("The 'MAXTRIALS' argument must be a finite non-negative",
                     "integer scalar.")
        cli_abort(c("x"=msg))
    }

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
        cached.scores <- .Call(C_sccache_new, ncol(x))
    else {
        if (!.load_suggested_package("RBGL")) {
            msg <- paste("The R engine requires the Bioconductor package",
                         "RBGL and it cannot be loaded.")
            cli_abort(c=("x"=msg))
        }
        cached.scores <- list()
        for (i in seq_len(ncol(x)))
            cached.scores[[i]] <- new.env(hash=TRUE, parent=emptyenv())
    }

    global.sufstats <- NULL
    global.sufstats.fun <- attr(scorefun, "global.sufstats.fun")
    if (!is.null(global.sufstats.fun)) {
        if (verbose)
          cli_alert_info("Calculating global sufficient statistics")
        global.sufstats <- global.sufstats.fun(x, targets, target.index)
    }

    anc <- init.ancestors(colnames(x))
    vidx <- setNames(seq_len(ncol(x)), colnames(x))
    ## nodes(dag) never changes during the search, only its edges do, so the
    ## vertex names and the nodes(dag)-position -> x-column map that
    ## translate a move's integer u/v are built once here
    vnames <- nodes(dag)
    vidx.nodes <- unname(vidx[vnames])
    pasets <- init.pasets(ncol(x))
    ## assuming .check_targets() has been called
    ## length(0:r) is computed here, not as r + 1 in C: 0:r for a
    ## non-integer r is 0:floor(r) and for a negative r counts down, so the
    ## coercion stays where it already behaves correctly (see src/rcar.c)
    rlen <- length(0:r)
    s0 <- -Inf
    ## g and x go POSITIONALLY: a user-supplied scorefun may still call its
    ## second formal 'dat', and renaming ours must not break theirs
    s1 <- scorefun(dag, x, targets=targets, target.index=target.index,
                   cached.scores=cached.scores, global.sufstats=global.sufstats)
    was_in_local_maximum <- local_maximum <- s1 < s0
    trials <- escapes <- avg_trials_per_escape <- 0
    ## Both exact paths can decline and hand the work back to the older,
    ## cruder machinery, and neither degradation is otherwise visible from the
    ## result: sampler.fallbacks counts draws that reverted to the rcar() walk,
    ## whose equilibrium is NOT uniform, and escape.fallbacks counts local
    ## maxima where the class could not be enumerated and the MAXTRIALS budget
    ## was used instead.
    sampler.fallbacks <- escape.fallbacks <- 0L
    ## Whether the exhaustive enumeration has already declined for the class
    ## the search is sitting in. The MAXTRIALS trials that follow a decline
    ## re-randomise the DAG WITHIN its I-equivalence class, so they change
    ## neither the class size nor its chain components -- the I-essential
    ## graph is an invariant of the class. Re-attempting the enumeration after
    ## each of them therefore recomputes a decision that cannot have changed,
    ## and counted the same refusal MAXTRIALS + 1 times. The flag is cleared
    ## whenever a move actually leaves the class.
    escape.declined <- FALSE
    ## The same for the exact sampler: its two decline conditions -- a chain
    ## component past 64 vertices, and counts too large to stay exact -- are
    ## properties of the I-essential graph, which is an invariant of the
    ## class. So once it has declined, it declines for every within-class
    ## trial that follows, and re-asking rebuilds the I-essential graph and
    ## the Clique-Picking contexts to reach a settled answer: 90 us a call at
    ## p = 200, against the 19 us of the walk it falls back to.
    ##
    ## Unlike the escape, though, every iteration really does make a draw, and
    ## it is the walk that makes it. So the DECISION is cached while
    ## sampler.fallbacks keeps counting one per draw -- the statistic is the
    ## number of draws that were not uniform, and that is one per iteration
    ## whether or not the reason was recomputed.
    sampler.declined <- FALSE

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
        st <- .Call(C_dag_new, ncol(x))
        while (!local_maximum) {
            s0 <- s1
            ## the exact draw declines for a chain component past the 64
            ## vertices it represents as a bitmask, or for counts too large
            ## to stay exact; the decision is cached per class, the draw is
            ## still counted per iteration
            drew <- FALSE
            if (sampler == "exact" && !sampler.declined) {
                drew <- .Call(C_dag_imec_sample, st, targets)
                if (!drew) sampler.declined <- TRUE
            }
            if (!drew) {
                if (sampler == "exact")
                    sampler.fallbacks <- sampler.fallbacks + 1L
                .Call(C_dag_rcar, st, rlen, targets)
            }
            ne <- .Call(C_dag_nh, st, 3L, targets)    ## 3 = ncr
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
                escape.declined <- sampler.declined <- FALSE
            } else if (escape == "exhaustive" && !escape.declined &&
                       !is.null(mm <- {
                           m0 <- .Call(C_dag_imec_members, st, targets,
                                       escape.maxD)
                           if (is.null(m0)) {
                               escape.fallbacks <- escape.fallbacks + 1L
                               escape.declined <- TRUE
                           }
                           m0
                       })) {
                ## every member of the I-equivalence class, and the best move
                ## available from any of them; no MAXTRIALS budget involved
                ## the arcs AND both orders: the members below run through
                ## the live DAG, and if none of them improves on where it
                ## started the original has to come back exactly, order
                ## included -- C_dag_set_edges() would canonicalise it
                here <- .Call(C_dag_edgeM, st)
                here.pas <- .Call(C_dag_pasets, st)
                best.s <- s0; best <- NULL
                for (em in mm) {
                    .Call(C_dag_set_edges, st, em)
                    ne.m <- .Call(C_dag_nh, st, 3L, targets)
                    am.m <- nh.argmax.fun(ne.m$op, vidx.nodes[ne.m$u],
                                          vidx.nodes[ne.m$v],
                                          .Call(C_dag_pasets, st), global.sufstats,
                                          cached.scores, verify.band,
                                          .Call(C_dag_pastamp, st))
                    if (am.m$total > best.s) {
                        best.s <- am.m$total
                        best <- list(em=em, op=ne.m$op[am.m$index],
                                     u=ne.m$u[am.m$index], v=ne.m$v[am.m$index])
                    }
                }
                if (is.null(best)) {
                    .Call(C_dag_restore_state, st, here, here.pas)
                    local_maximum <- TRUE
                    s1 <- s0
                } else {
                    .Call(C_dag_set_edges, st, best$em)
                    .Call(C_dag_apply_move, st, best$op, best$u, best$v)
                    s1 <- best.s
                    local_maximum <- FALSE
                    escapes <- escapes + 1
                    was_in_local_maximum <- FALSE
                    trials <- 0
                    escape.declined <- sampler.declined <- FALSE
                }
            } else if (trials < MAXTRIALS) {
                s1 <- s0
                drew <- FALSE
                if (sampler == "exact" && !sampler.declined) {
                    drew <- .Call(C_dag_imec_sample, st, targets)
                    if (!drew) sampler.declined <- TRUE
                }
                if (!drew) {
                    if (sampler == "exact")
                        sampler.fallbacks <- sampler.fallbacks + 1L
                    .Call(C_dag_rcar, st, rlen, targets)
                }
                local_maximum <- FALSE
                was_in_local_maximum <- TRUE
                trials <- trials + 1
            } else
                s1 <- s0

            .debug_assertions(st, NULL, NULL, NULL, x, vnames)

            if (verbose)
                cli_progress_update()
        }
        dag <- .graphNEL_from_edgeM(vnames, .Call(C_dag_edgeM, st))
    } else {
        while (!local_maximum) {
            s0 <- s1
            rcar.out <- if (sampler == "exact" && !sampler.declined)
                isample.move(dag, targets, anc, pasets, vidx, vnames,
                             vidx.nodes, r)
            else rcar(dag, r, targets, anc, pasets, vidx)
            if (isTRUE(rcar.out$fallback)) sampler.declined <- TRUE
            if (sampler == "exact" && sampler.declined)
                sampler.fallbacks <- sampler.fallbacks + 1L
            dag <- rcar.out$dag
            anc <- rcar.out$anc
            pasets <- rcar.out$pasets
            ne <- ncr.nh(dag, anc, targets)
            sco <- score.nh(ne, dag, x, targets, target.index, cached.scores,
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
                escape.declined <- sampler.declined <- FALSE
            } else if (escape == "exhaustive" && !escape.declined &&
                       !is.null(mm <- {
                           m0 <- imec.members(dag, targets, vnames, vidx.nodes,
                                              escape.max)
                           if (is.null(m0)) {
                               escape.fallbacks <- escape.fallbacks + 1L
                               escape.declined <- TRUE
                           }
                           m0
                       })) {
                ## Examine EVERY member of the current I-equivalence class and
                ## take the best move available from any of them.  With the
                ## class size known exactly this replaces the MAXTRIALS budget
                ## by a deterministic and complete check: if nothing here beats
                ## s0, the search really is at a local maximum of the class.
                best.s <- s0; best <- NULL
                for (M in mm) {
                    ne.m <- ncr.nh(M$dag, M$anc, targets)
                    sco.m <- score.nh(ne.m, M$dag, x, targets, target.index,
                                      cached.scores, global.sufstats, M$pasets,
                                      vidx.nodes, supports.pasets, scorefun,
                                      nh.scores.fun)
                    bm <- which.max(sco.m)
                    if (sco.m[bm] > best.s) {
                        best.s <- sco.m[bm]
                        best <- list(M=M, op=ne.m$op[bm], u=ne.m$u[bm], v=ne.m$v[bm])
                    }
                }
                if (is.null(best)) {
                    local_maximum <- TRUE
                    s1 <- s0
                } else {
                    anc <- switch(best$op,
                                  add.ancestors(best$M$anc, vnames[best$u], vnames[best$v]),
                                  remove.ancestors(best$M$anc, best$M$dag,
                                                   vnames[best$u], vnames[best$v]),
                                  reverse.ancestors(best$M$anc, best$M$dag,
                                                    vnames[best$u], vnames[best$v]))
                    pasets <- move.pasets(best$M$pasets, best$op,
                                          vidx.nodes[best$u], vidx.nodes[best$v])
                    dag <- apply.move(best$M$dag, best$op, best$u, best$v, vnames)
                    s1 <- best.s
                    local_maximum <- FALSE
                    escapes <- escapes + 1
                    was_in_local_maximum <- FALSE
                    trials <- 0
                    escape.declined <- sampler.declined <- FALSE
                }
            } else if (trials < MAXTRIALS) {
                s1 <- s0
                rcar.out <- if (sampler == "exact" && !sampler.declined)
                    isample.move(dag, targets, anc, pasets, vidx, vnames,
                                 vidx.nodes, r)
                else rcar(dag, r, targets, anc, pasets, vidx)
                if (isTRUE(rcar.out$fallback)) sampler.declined <- TRUE
                if (sampler == "exact" && sampler.declined)
                    sampler.fallbacks <- sampler.fallbacks + 1L
                dag <- rcar.out$dag
                anc <- rcar.out$anc
                pasets <- rcar.out$pasets
                local_maximum <- FALSE
                was_in_local_maximum <- TRUE
                trials <- trials + 1
            } else
                s1 <- s0

            .debug_assertions(NULL, dag, anc, pasets, x, vnames)

            if (verbose)
                cli_progress_update()
        }
    }

    if (verbose) {
        algname <- if (identical(targets, list(integer(0)))) "HCMC" else "iHCMC"
        cli_progress_done("{algname} algorithm completed")
    }

    if (verbose && (sampler.fallbacks > 0L || escape.fallbacks > 0L)) {
        msg <- paste("exact sampling fell back to RCAR {sampler.fallbacks}",
                     "time{?s}; exhaustive enumeration fell back to the",
                     "MAXTRIALS budget {escape.fallbacks} time{?s}")
        cli_alert_warning(msg)
    }

    list(dag=dag, sco=s1, sampler.fallbacks=sampler.fallbacks,
         escape.fallbacks=escape.fallbacks)
}

#' @importFrom cli cli_abort
.check_input_data <- function(dat) {
    ## a population model stands in for the data; it reports ncol() and
    ## colnames() like the matrix it replaces, so nothing downstream changes
    if (.is.population(dat))
        return(.as.population(dat))
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

#' @importFrom cli cli_abort
.check_targets <- function(targets, p) {
    if (!is.list(targets)) {
        msg <- paste("The 'targets' argument must be a list of integer",
                     "vectors, each vector specifying zero or more",
                     "intervened random variables in the data.")
        cli_abort(c("x"=msg))
    } else if (length(targets) == 0) {
        targets <- list(integer(0))
        msg <- paste("The 'targets' list is empty, replaced by",
                     "'list(integer(0))', which implies all data",
                     "are observational.")
        cli_alert_warning(msg)
        return(targets)
    }

    for (i in seq_along(targets)) {
        if (!is.numeric(targets[[i]]) || any(!is.finite(targets[[i]])) ||
            any(is.na(targets[[i]]))) {
            msg <- paste("Each element of the 'targets' list must be an",
                         "vector of finite positive integer values.")
            cli_abort(c("x"=msg))
        }
        if (any(targets[[i]] < 1 | targets[[i]] > p)) {
            msg <- paste("Each target in the 'targets' list must be an",
                         "integer between 1 and", p, "(inclusive).")
            cli_abort(c("x"=msg))
        }
        if (!is.integer(targets[[i]])) {
            if (any(targets[[i]] != floor(targets[[i]]))) {
                msg <- paste("Each element of the 'targets' list must be an",
                             "integer vector.")
                cli_abort(c("x"=msg))
            }
            targets[[i]] <- as.integer(targets[[i]])
        }
    }

    targets
}

.check_escape.max <- function(escape.max) {
    if (!is.numeric(escape.max) || length(escape.max) != 1 ||
        is.nan(escape.max) ||(!is.na(escape.max) && (escape.max <= 0 ||
        escape.max != floor(escape.max)))) {
        msg <- paste("The 'escape.max' argument must be either a positive",
                     "integer scalar or NA/Inf, which imply no limit.")
        cli_abort(c("x"=msg))
    }

    as.double(escape.max)
}
