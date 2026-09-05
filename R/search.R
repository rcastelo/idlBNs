##
## ANCESTOR MATRIX MAINTENANCE (Giudici and Castelo, ML, 2003)
##
## an ancestor matrix 'anc' is a p x p logical matrix, dimnames=list(v, v),
## where v are the DAG vertex names, and anc[u, v] is TRUE iff u is an
## ancestor of v (a directed path u -> ... -> v exists). It is maintained
## incrementally across the whole search, instead of being recomputed from
## scratch on every neighborhood-generation call, and is used to test
## whether a candidate edge addition/reversal preserves acyclicity with O(1)
## (addition) or O(degree) (reversal) lookups.

## build an all-FALSE ancestor matrix for an edgeless DAG on vertex names v
init.ancestors <- function(v) {
    p <- length(v)
    matrix(FALSE, nrow=p, ncol=p, dimnames=list(v, v))
}

## incremental update of 'anc' after adding edge u -> v
add.ancestors <- function(anc, u, v) {
    anc[u, v] <- TRUE
    delta <- anc[, u]
    delta[u] <- TRUE
    D <- c(v, rownames(anc)[anc[v, ]]) ## v and its descendants
    anc[, D] <- anc[, D] | delta
    anc
}

## incremental update of 'anc' after removing edge u -> v from 'dag', where
## 'dag' is the DAG as it stood BEFORE the removal (used only to read parent
## sets and a topological order)

#' @importFrom graph nodes edgeMatrix
#' @importFrom RBGL tsort
remove.ancestors <- function(anc, dag, u, v) {
    vnodes <- nodes(dag)
    em <- edgeMatrix(dag)
    pasets <- split(vnodes[em["from", ]], factor(vnodes[em["to", ]],
                                                 levels=vnodes))
    D <- c(v, rownames(anc)[anc[v, ]]) ## v and its descendants
    ## any topological order of 'dag' remains valid after removing an edge
    to <- tsort(dag)
    D <- to[to %in% D]
    for (k in D) {
        pa.k <- pasets[[k]]
        if (k == v)
            pa.k <- setdiff(pa.k, u)
        newcol <- rep(FALSE, length(vnodes))
        names(newcol) <- vnodes
        if (length(pa.k) > 0) {
            newcol[pa.k] <- TRUE
            newcol <- newcol | apply(anc[, pa.k, drop=FALSE], 1, any)
        }
        anc[, k] <- newcol
    }
    anc
}

## incremental update of 'anc' after reversing edge u -> v into v -> u in
## 'dag', where 'dag' is the DAG as it stood BEFORE the reversal
reverse.ancestors <- function(anc, dag, u, v) {
    anc <- remove.ancestors(anc, dag, u, v)
    add.ancestors(anc, v, u)
}

##
## PARENT-SET MAINTENANCE
##
## a pasets structure is a list of length p, one element per vertex in
## dat-column order, where pasets[[i]] is an integer vector of the
## dat-column-index parents of vertex i (the same representation iBIC()/
## iBGe() build internally from a graph via .build_pasets()). It is
## maintained incrementally across the whole search, instead of being
## recomputed from scratch (match()+edgeMatrix()+split()) on every call to
## a score function. Unlike the ancestor matrix, parent-set membership is
## NOT transitive, so an edge operation between u and v only ever touches
## the 1-2 list elements for u and/or v directly -- no cascading
## recomputation of other vertices, no topological order needed. 'u'/'v'
## here are integer dat-column indices, not vertex names; callers translate
## a neighbor's (character) u/v via a name -> index lookup built once per
## search (see hillclimbing()/hcmc()).

## build an all-empty pasets list for an edgeless DAG on p vertices
init.pasets <- function(p) replicate(p, integer(0), simplify=FALSE)

## incremental update of 'pasets' after adding edge u -> v
add.pasets <- function(pasets, u, v) {
    stopifnot(!(u %in% pasets[[v]]))
    pasets[[v]] <- c(pasets[[v]], u)
    pasets
}

## incremental update of 'pasets' after removing edge u -> v
remove.pasets <- function(pasets, u, v) {
    stopifnot(u %in% pasets[[v]])
    pasets[[v]] <- setdiff(pasets[[v]], u)
    pasets
}

## incremental update of 'pasets' after reversing edge u -> v into v -> u
reverse.pasets <- function(pasets, u, v) {
    stopifnot(u %in% pasets[[v]])
    pasets[[v]] <- setdiff(pasets[[v]], u)
    pasets[[u]] <- c(pasets[[u]], v)
    pasets
}

##
## NEIGHBORHOODS (Castelo and Kocka, JMLR, 2003)
##

## NR: non-reversals neighborhood (addition and removal only)
## each returned entry is a list(graph=, op=, u=, v=) describing the move

#' @importFrom graph edgeL removeEdge addEdge nodes
nr.nh <- function(dag, anc) {
    v <- nodes(dag)
    e <- edgeL(dag)
    nr <- list()
    nr.i <- 0
    for (i in seq_along(e)) {
        a <- v[e[[i]]$edges]
        na <- setdiff(v, c(a, names(e)[i])) ## exclude i itself (no self-loops)
        for (j in seq_along(na)) { ## go through non-adjacent vertices
            if (!anc[na[j], names(e)[i]]) {
                tmp.g <- addEdge(names(e)[i], na[j], dag)
                nr.i <- nr.i + 1
                nr[[nr.i]] <- list(graph=tmp.g, op="add", u=names(e)[i], v=na[j])
            }
        }
        for (j in seq_along(a)) { ## go through adjacent vertices
            tmp.g <- removeEdge(names(e)[i], a[j], dag)
            nr.i <- nr.i + 1
            nr[[nr.i]] <- list(graph=tmp.g, op="remove", u=names(e)[i], v=a[j])
        }
    }
    nr
}

## AR: all-reversals neighborhood (NR + all-arc-reversals)

#' @importFrom graph edgeL removeEdge addEdge nodes
ar.nh <- function(dag, anc) {
    v <- nodes(dag)
    e <- edgeL(dag)
    ar <- nr.nh(dag, anc)
    ar.i <- length(ar)
    for (i in seq_along(e)) { ## reverse edges
        a <- v[e[[i]]$edges]
        for (j in seq_along(a)) { ## go through adjacent vertices
            if (!any(anc[a[-j], a[j]])) {
                tmp.g <- removeEdge(names(e)[i], a[j], dag)
                tmp.g <- addEdge(a[j], names(e)[i], tmp.g)
                ar.i <- ar.i + 1
                ar[[ar.i]] <- list(graph=tmp.g, op="reverse", u=names(e)[i], v=a[j])
            }
        }
    }
    ar
}

## NCR: non-covered arc reversals neighborhood (NR + non-covered-arc-reversals)

#' @importFrom graph edgeL removeEdge addEdge edgeMatrix nodes
ncr.nh <- function(dag, anc, utargets=integer(0)) {
    v <- nodes(dag)
    e <- edgeL(dag)
    ncr <- nr.nh(dag, anc)
    ncr.i <- length(ncr)
    em <- edgeMatrix(dag)
    pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
    for (i in seq_along(e)) { ## reverse edges
        a <- v[e[[i]]$edges]
        for (j in seq_along(a)) { ## go through adjacent vertices
            ced <- identical(sort(pasets[[names(e)[i]]]), sort(setdiff(pasets[[a[j]]], names(e)[i])))
            if (!ced || any(c(e[[i]]$edges[j], i) %in% utargets)) { ## NCR including not interventionally covered
                if (!any(anc[a[-j], a[j]])) {
                    tmp.g <- removeEdge(names(e)[i], a[j], dag)
                    tmp.g <- addEdge(a[j], names(e)[i], tmp.g)
                    ncr.i <- ncr.i + 1
                    ncr[[ncr.i]] <- list(graph=tmp.g, op="reverse", u=names(e)[i], v=a[j])
                }
            }
        }
    }
    ncr
}

##
## REPEATED COVERED ARC REVERSAL ALGORITHM (Castelo and Kocka, JMLR, 2003)
## ADAPTED TO INTERVENTIONS IN (Castelo, 2025)
##

## build a logical mask indicated what edges are "covered" in the input DAG
## utargets should be a vector of unique target vertices, which when non-empty
## restricts covered edges to those without any target vertex

#' @importFrom graph nodes edgeMatrix
cedges <- function(dag, utargets) {
    v <- nodes(dag)
    em <- edgeMatrix(dag)
    pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
    cemask <- mapply(function(pafrom, pato, from) identical(sort(pafrom), sort(setdiff(pato, from))),
                     pasets[em["from", ]], pasets[em["to", ]], v[em["from", ]])
    temask <- rep(FALSE, ncol(em))
    if (length(utargets) > 0)
        temask <- colSums(matrix(as.vector(em) %in% utargets, ncol=ncol(em))) > 0
    cemask & !temask
}

## resample helper function
resample <- function(x, ...) x[sample.int(length(x), ...)]

## RCAR: repeated covered arc reversal algorithm
## utargets should be a vector of unique target vertices
## returns a list(dag=, anc=, pasets=) since every reversal it performs,
## although always cycle-safe by construction (a covered edge cannot
## introduce a cycle), still changes true ancestor relationships and parent
## sets, and must keep both 'anc' and 'pasets' in sync for subsequent
## neighborhood generation and scoring to remain correct. 'vidx' is a named
## integer vector mapping vertex name -> dat-column index (as used by
## 'pasets'), built once per search by the caller.

#' @importFrom graph removeEdge addEdge numEdges edgeMatrix nodes
rcar <- function(dag, r, utargets, anc, pasets, vidx) {
    if (numEdges(dag) == 0)
        return(list(dag=dag, anc=anc, pasets=pasets))
    cemask <- cedges(dag, utargets)
    if (!any(cemask))
        return(list(dag=dag, anc=anc, pasets=pasets))

    tmp.g <- dag
    v <- nodes(tmp.g)
    rr <- sample(0:r, size=1)
    for (i in seq_len(rr)) {
        em <- edgeMatrix(tmp.g)
        cemask <- cedges(tmp.g, utargets)
        rndce <- resample(which(cemask), size=1)
        u <- v[em["from", rndce]]
        w <- v[em["to", rndce]]
        anc <- reverse.ancestors(anc, tmp.g, u, w) ## a covered edge cannot introduce a cycle
        pasets <- reverse.pasets(pasets, vidx[[u]], vidx[[w]])
        tmp.g <- removeEdge(u, w, tmp.g)
        tmp.g <- addEdge(w, u, tmp.g)
    }
    list(dag=tmp.g, anc=anc, pasets=pasets)
}



##
## SEARCH ALGORITHMS OTHER THAN (i)HCMC
##

#' @title Straightforward (classical) hill-climbing algorithm
#'
#' @description Learn the structure of a Bayesian network from observational
#' and interventional data using a straightforward (classical) hill-climbing
#' algorithm that at each step during the search adds, removes and reverses all
#' possible arcs.
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
#' @param scorefun (Default is [`iBIC`]) A function to calculate the goodness
#' of fit (GoF) score of a DAG on a given data set.
#'
#' @param verbose (Default TRUE) Show progress in the calculations.
#'
#' @return A list containing a [`graphNEL`][graph::graphNEL-class] object with
#' the structure of the learned DAG, and its corresponding score.
#'
#' @seealso [iBIC()], [iBGe()]
#'
#' @importFrom graph graphNEL
#' @importClassesFrom graph graphNEL
#' @importFrom cli cli_progress_step cli_progress_update
#' @importFrom stats setNames
#' @export
hillclimbing <- function(dat, targets=list(integer(0)),
                         target.index=rep(1L, nrow(dat)),  scorefun=iBIC,
                         verbose=TRUE) {

    dat <- .check_input_data(dat)
    dag <- graphNEL(colnames(dat), edgemode="directed")
    attr(dat, "sanitycheck") <- TRUE

    stopifnot(is.list(targets)) ## QC
    scorefun <- match.fun(scorefun)

    cached.scores <- list()
    for (i in seq_len(ncol(dat)))
        cached.scores[[i]] <- new.env(hash=TRUE, parent=emptyenv())

    global.sufstats <- NULL
    global.sufstats.fun <- attr(scorefun, "global.sufstats.fun")
    if (!is.null(global.sufstats.fun)) {
        if (verbose)
            cli_alert_info("Calculating global sufficient statistics")
        global.sufstats <- global.sufstats.fun(dat, targets, target.index)
    }
    scorefun.name <- NULL
    if (!is.null(attr(scorefun, "scorefun.name")))
        scorefun.name <- attr(scorefun, "scorefun.name")
    supports.pasets <- isTRUE(attr(scorefun, "supports.pasets"))

    anc <- init.ancestors(colnames(dat))
    vidx <- setNames(seq_len(ncol(dat)), colnames(dat))
    pasets <- init.pasets(ncol(dat))

    s0 <- -Inf
    s1 <- scorefun(g=dag, dat=dat, targets=targets, target.index=target.index,
                   cached.scores=cached.scores, global.sufstats=global.sufstats)

    if (verbose) {
        msg <- "Running a straightforward hill-climbing algorithm"
        if (!is.null(scorefun.name))
            msg <- paste(msg, "with the {scorefun.name} score function")
        cli_progress_bar(msg)
        cli_progress_step("Score {s1}", spinner=TRUE)
    }

    while (s1 > s0) {
        s0 <- s1
        ne <- ar.nh(dag, anc)
        s1 <- sapply(ne, function(nb, d, tgts, tgt.idx, chd.sco, gbl.sst,
                                  pas, vix, use.pas) {
                          args <- list(g=nb$graph, dat=d, targets=tgts,
                                       target.index=tgt.idx, cached.scores=chd.sco,
                                       global.sufstats=gbl.sst)
                          if (use.pas)
                              args$pasets <- switch(nb$op,
                                                    add     = add.pasets(pas, vix[[nb$u]], vix[[nb$v]]),
                                                    remove  = remove.pasets(pas, vix[[nb$u]], vix[[nb$v]]),
                                                    reverse = reverse.pasets(pas, vix[[nb$u]], vix[[nb$v]]))
                          do.call(scorefun, args)
                      },
                     dat, targets, target.index, cached.scores, global.sufstats,
                     pasets, vidx, supports.pasets)
        best <- ne[[which.max(s1)]]
        anc <- switch(best$op,
                      add     = add.ancestors(anc, best$u, best$v),
                      remove  = remove.ancestors(anc, dag, best$u, best$v),
                      reverse = reverse.ancestors(anc, dag, best$u, best$v))
        pasets <- switch(best$op,
                         add     = add.pasets(pasets, vidx[[best$u]], vidx[[best$v]]),
                         remove  = remove.pasets(pasets, vidx[[best$u]], vidx[[best$v]]),
                         reverse = reverse.pasets(pasets, vidx[[best$u]], vidx[[best$v]]))
        dag <- best$graph
        s1 <- max(s1)

        if (isTRUE(getOption("idlBNs.debug.pasets", FALSE)))
            stopifnot(identical(unname(lapply(pasets, function(x) unname(sort.int(x)))),
                                unname(lapply(.build_pasets(dag, dat), function(x) unname(sort.int(x))))))

        if (verbose)
          cli_progress_update()
    }

    if (verbose)
        cli_progress_done("straightforward hill-climbing algorithm completed")

    list(dag=dag, sco=s1)
}
