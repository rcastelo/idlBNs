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
## a neighborhood is returned as a list(op=, u=, v=) of three parallel
## integer vectors of the same length k, describing k candidate moves on the
## input DAG: op[m] is one of the OP.* codes below, and u[m], v[m] are
## indices into nodes(dag) naming the arc u -> v that the move acts on (for
## OP.REVERSE, the arc as it stands BEFORE the reversal). no candidate graph
## is built here: a search step scores the whole neighborhood from the
## (op, u, v) parent-set deltas alone (see score.nh()) and materialises a
## single graph, for the winning move, via apply.move(). building one
## graphNEL per candidate used to dominate the whole search -- addEdge() on
## a p=30 DAG costs ~60us against ~1us for the C score of one vertex, and
## every step discards all but one of its O(p^2) candidate graphs.

OP.ADD     <- 1L
OP.REMOVE  <- 2L
OP.REVERSE <- 3L

## apply a single (op, u, v) move to 'dag' and return the resulting DAG.
## 'vnames' is nodes(dag), taken as an argument so a caller inside the
## search loop passes the copy it already holds instead of re-fetching it.

#' @importFrom graph addEdge removeEdge nodes
apply.move <- function(dag, op, u, v, vnames=nodes(dag)) {
    switch(op,
           addEdge(vnames[u], vnames[v], dag),                  ## OP.ADD
           removeEdge(vnames[u], vnames[v], dag),               ## OP.REMOVE
           addEdge(vnames[v], vnames[u],                        ## OP.REVERSE
                   removeEdge(vnames[u], vnames[v], dag)))
}

## build the parent sets of 'dag' as a list of p integer vectors of vertex
## indices, from the child lists that edgeL() already returns, avoiding a
## second traversal of the graph through edgeMatrix()
.pasets.from.edgeL <- function(e, p) {
    ch <- lapply(e, `[[`, "edges")
    to <- unlist(ch, use.names=FALSE)
    if (is.null(to))
        return(replicate(p, integer(0), simplify=FALSE))
    from <- rep(seq_len(p), lengths(ch))
    unname(split(from, factor(to, levels=seq_len(p))))
}

## NR: non-reversals neighborhood (addition and removal only)

#' @importFrom graph edgeL nodes
nr.nh <- function(dag, anc) {
    p <- length(nodes(dag))
    e <- edgeL(dag)
    kmax <- p * (p - 1L) ## p-1 additions + at most p-1 removals per vertex
    op <- integer(kmax)
    uu <- integer(kmax)
    vv <- integer(kmax)
    k <- 0L
    for (i in seq_len(p)) {
        a <- e[[i]]$edges                ## children of i
        na <- seq_len(p)[-c(a, i)]       ## non-adjacent, i excluded (no self-loops)
        ## i -> na[j] closes a cycle iff na[j] is already an ancestor of i
        na <- na[!anc[na, i]]
        if (length(na) > 0L) {           ## additions
            idx <- k + seq_along(na)
            op[idx] <- OP.ADD
            uu[idx] <- i
            vv[idx] <- na
            k <- k + length(na)
        }
        if (length(a) > 0L) {            ## removals
            idx <- k + seq_along(a)
            op[idx] <- OP.REMOVE
            uu[idx] <- i
            vv[idx] <- a
            k <- k + length(a)
        }
    }
    idx <- seq_len(k)

    list(op=op[idx], u=uu[idx], v=vv[idx])
}

## which of the arcs i -> a[j] can be reversed without closing a cycle: the
## reversal a[j] -> i creates a cycle iff some other child of i is a
## descendant of a[j], i.e. iff any(anc[a[-j], a[j]]). computed for all j at
## once off the |a| x |a| submatrix of 'anc', as colSums() minus the
## diagonal, rather than one any() call per j
.reversible <- function(anc, a) {
    M <- anc[a, a, drop=FALSE]

    (colSums(M) - diag(M)) == 0
}

## append the reversals of the arcs i -> a[keep] to the move vectors held in
## the 'acc' accumulator environment
.add.reversals <- function(acc, i, a, keep) {
    if (!any(keep))
        return(invisible(NULL))
    a <- a[keep]
    acc$n <- acc$n + 1L
    acc$u[[acc$n]] <- rep(i, length(a))
    acc$v[[acc$n]] <- a
    invisible(NULL)
}

## combine an NR neighborhood with the reversals gathered in 'acc'
.bind.reversals <- function(nr, acc) {
    if (acc$n == 0L)
        return(nr)
    ru <- unlist(acc$u[seq_len(acc$n)], use.names=FALSE)
    rv <- unlist(acc$v[seq_len(acc$n)], use.names=FALSE)

    list(op=c(nr$op, rep(OP.REVERSE, length(ru))),
         u=c(nr$u, ru),
         v=c(nr$v, rv))
}

## AR: all-reversals neighborhood (NR + all-arc-reversals)

#' @importFrom graph edgeL nodes
ar.nh <- function(dag, anc) {
    p <- length(nodes(dag))
    e <- edgeL(dag)
    nr <- nr.nh(dag, anc)
    acc <- new.env(parent=emptyenv())
    acc$n <- 0L
    acc$u <- vector("list", p)
    acc$v <- vector("list", p)
    for (i in seq_len(p)) {
        a <- e[[i]]$edges
        if (length(a) == 0L)
            next
        .add.reversals(acc, i, a, .reversible(anc, a))
    }

    .bind.reversals(nr, acc)
}

## NCR: non-covered arc reversals neighborhood (NR + non-covered-arc-reversals)

#' @importFrom graph edgeL nodes
ncr.nh <- function(dag, anc, utargets=integer(0)) {
    p <- length(nodes(dag))
    e <- edgeL(dag)
    nr <- nr.nh(dag, anc)
    pasets <- .pasets.from.edgeL(e, p)
    acc <- new.env(parent=emptyenv())
    acc$n <- 0L
    acc$u <- vector("list", p)
    acc$v <- vector("list", p)
    for (i in seq_len(p)) {
        a <- e[[i]]$edges
        if (length(a) == 0L)
            next
        ## arc i -> a[j] is covered iff pa(i) == pa(a[j]) \ {i}; sort pa(i)
        ## once for all j rather than once per j
        pa.i <- sort.int(pasets[[i]])
        ced <- vapply(a,
                      function(w) identical(pa.i,
                                            sort.int(setdiff(pasets[[w]], i))),
                      logical(1))
        ## an I-covered arc is one that is covered AND has no target vertex
        ## at either endpoint, so a covered arc touching a target stays in
        ## the neighborhood
        tgt <- (i %in% utargets) | (a %in% utargets)
        .add.reversals(acc, i, a, (!ced | tgt) & .reversible(anc, a))
    }

    .bind.reversals(nr, acc)
}

## apply a move's parent-set delta to 'pasets'. 'pu'/'pv' are the move's
## endpoints already translated to the dat-column indices 'pasets' is keyed
## by (see vidx.nodes in hcmc()/hillclimbing()).

#' @importFrom cli cli_abort
move.pasets <- function(pasets, op, pu, pv) {
    if (length(op) != 1L || is.na(op))
        cli_abort(c=("move.pasets: unknown operation code {op}"))
    
    if (op == OP.ADD)
        return(add.pasets(pasets, pu, pv))

    if (op == OP.REMOVE)
        return(remove.pasets(pasets, pu, pv))

    if (op == OP.REVERSE)
        return(reverse.pasets(pasets, pu, pv))

    cli_abort(c=("move.pasets: unknown operation code {op}"))
}

## score every candidate move of a neighborhood 'ne' against the DAG it was
## generated from, returning a numeric vector of length |ne| in neighborhood
## order.
##
## when 'scorefun' carries a neighborhood scorer (attr "nh.scores.fun"), the
## whole neighborhood is scored in a single call that rescores only the one
## or two vertices whose parent set each move changes, reusing the current
## DAG's terms for the other p-1 or p-2 (see .iBIC.nh.scores() and
## src/nh_scores.c). each candidate's terms are still summed over all p
## vertices, so the totals are bit-identical to scoring every candidate
## from scratch and the search follows exactly the same trajectory.
##
## failing that, when 'scorefun' declares that it accepts a 'pasets'
## argument (attr "supports.pasets"), the move's parent-set delta is applied
## to 'pasets' and the *current* 'dag' is passed as 'g': the score functions
## read 'g' only to validate the vertex count and the cached-scores list,
## and every neighbor shares both with 'dag', so no candidate graph is
## needed. a custom scorefun with neither attribute still gets a real
## neighbor graph, built here on demand. 'vidx.nodes' maps a nodes(dag)
## position to the dat-column index that 'pasets' is keyed by.
score.nh <- function(ne, dag, dat, targets, target.index, cached.scores,
                     global.sufstats, pasets, vidx.nodes, supports.pasets,
                     scorefun, nh.scores.fun=NULL) {
    if (!is.null(nh.scores.fun))
        return(nh.scores.fun(ne$op, vidx.nodes[ne$u], vidx.nodes[ne$v],
                             pasets, global.sufstats, cached.scores))

    k <- length(ne$op)
    sco <- numeric(k)
    vnames <- if (supports.pasets) NULL else nodes(dag)
    for (m in seq_len(k)) {
        if (supports.pasets) {
            pas <- move.pasets(pasets, ne$op[m], vidx.nodes[ne$u[m]],
                               vidx.nodes[ne$v[m]])
            sco[m] <- scorefun(g=dag, dat=dat, targets=targets,
                               target.index=target.index,
                               cached.scores=cached.scores,
                               global.sufstats=global.sufstats, pasets=pas)
        } else {
            g <- apply.move(dag, ne$op[m], ne$u[m], ne$v[m], vnames)
            sco[m] <- scorefun(g=g, dat=dat, targets=targets,
                               target.index=target.index,
                               cached.scores=cached.scores,
                               global.sufstats=global.sufstats)
        }
    }

    sco
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
    nh.scores.fun <- attr(scorefun, "nh.scores.fun")

    anc <- init.ancestors(colnames(dat))
    vidx <- setNames(seq_len(ncol(dat)), colnames(dat))
    ## nodes(dag) never changes during the search, only its edges do, so the
    ## vertex names and the nodes(dag)-position -> dat-column map that
    ## translate a move's integer u/v are built once here
    vnames <- nodes(dag)
    vidx.nodes <- unname(vidx[vnames])
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
        sco <- score.nh(ne, dag, dat, targets, target.index, cached.scores,
                        global.sufstats, pasets, vidx.nodes, supports.pasets,
                        scorefun, nh.scores.fun)
        b <- which.max(sco)
        b.op <- ne$op[b]
        b.u <- ne$u[b]
        b.v <- ne$v[b]
        ## 'anc' and 'pasets' are updated before 'dag', since
        ## remove.ancestors()/reverse.ancestors() read the DAG as it stood
        ## BEFORE the move
        anc <- switch(b.op,
                      add.ancestors(anc, vnames[b.u], vnames[b.v]),
                      remove.ancestors(anc, dag, vnames[b.u], vnames[b.v]),
                      reverse.ancestors(anc, dag, vnames[b.u], vnames[b.v]))
        pasets <- move.pasets(pasets, b.op, vidx.nodes[b.u], vidx.nodes[b.v])
        dag <- apply.move(dag, b.op, b.u, b.v, vnames)
        s1 <- sco[b]

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
