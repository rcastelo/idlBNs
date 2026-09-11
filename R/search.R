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
remove.ancestors <- function(anc, dag, u, v) {
    vnodes <- nodes(dag)
    em <- edgeMatrix(dag)
    pasets <- split(vnodes[em["from", ]], factor(vnodes[em["to", ]],
                                                 levels=vnodes))
    D <- c(v, rownames(anc)[anc[v, ]]) ## v and its descendants
    ## any topological order of 'dag' remains valid after removing an edge
    to <- RBGL::tsort(dag)
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
ncr.nh <- function(dag, anc, targets=list()) {
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
        ## an I-covered arc is covered AND unseparated by every target, so a
        ## covered arc whose endpoints some target tells apart stays in the
        ## neighborhood (see .separates())
        tgt <- .separates(targets, rep(i, length(a)), a)
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
##
## 'dag' may be NULL, but only when 'nh.scores.fun' is supplied: that branch
## returns before 'dag' is read. The C search engine relies on this, because
## it holds no graphNEL while the search is running.
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
            sco[m] <- scorefun(dag, dat, targets=targets,
                               target.index=target.index,
                               cached.scores=cached.scores,
                               global.sufstats=global.sufstats, pasets=pas)
        } else {
            g <- apply.move(dag, ne$op[m], ne$u[m], ne$v[m], vnames)
            sco[m] <- scorefun(g, dat, targets=targets,
                               target.index=target.index,
                               cached.scores=cached.scores,
                               global.sufstats=global.sufstats)
        }
    }

    sco
}

## rebuild a graphNEL from an edge matrix, replaying the arcs in column
## order. graph::addEdge appends, so this reproduces the @edgeL slot of the
## graph the edge matrix came from, element for element -- which matters,
## because a vertex's child order is what nr.nh() emits its removals in.
## Used to turn the C search state back into the graphNEL that hcmc() and
## hillclimbing() return, and to materialise it for the debug assertions.

#' @importFrom graph graphNEL addEdge
.graphNEL_from_edgeM <- function(vnames, em) {
    g <- graphNEL(vnames, edgemode="directed")
    for (k in seq_len(ncol(em)))
        g <- addEdge(vnames[em["from", k]], vnames[em["to", k]], g)

    g
}

## the in-loop debug assertions, shared by both engines and by both search
## algorithms. 'st' is the C state when the C engine is running, NULL
## otherwise; 'dag'/'anc'/'pasets' are the R state.
##
## idlBNs.debug.pasets keeps exactly the semantics it had when the search was
## pure R: the incrementally maintained parent sets must equal a from-scratch
## rebuild off the graph. Note it sorts both sides, so it is deliberately
## blind to parent ORDER -- which is why the C engine's parent order is
## pinned by tests/test_c_dag.R with identical() instead.
##
## idlBNs.debug.anc and idlBNs.debug.dag are new, and close the gap that
## remove.ancestors()/reverse.ancestors() had no in-loop assertion at all:
## the first compares the ancestor matrix against a from-scratch transitive
## closure, the second runs the C structure's own full self-consistency pass.
.debug_assertions <- function(st, dag, anc, pasets, dat, vnames) {
    if (isTRUE(getOption("idlBNs.debug.pasets", FALSE))) {
        g <- if (is.null(dag))
                 .graphNEL_from_edgeM(vnames, .Call(C_dag_edgeM, st))
             else dag
        pas <- if (is.null(pasets)) .Call(C_dag_pasets, st) else pasets
        stopifnot(identical(unname(lapply(pas,
                                          function(x) unname(sort.int(x)))),
                            unname(lapply(.build_pasets(g, dat),
                                          function(x) unname(sort.int(x))))))
    }
    if (isTRUE(getOption("idlBNs.debug.anc", FALSE))) {
        g <- if (is.null(dag))
                 .graphNEL_from_edgeM(vnames, .Call(C_dag_edgeM, st))
             else dag
        a <- if (is.null(anc)) .Call(C_dag_anc, st) else unname(anc)
        stopifnot(identical(a, .anc_closure(g)))
    }
    if (isTRUE(getOption("idlBNs.debug.dag", FALSE)) && !is.null(st))
        .Call(C_dag_check, st)

    invisible(NULL)
}

## the ancestor matrix computed from scratch: the transitive closure of the
## adjacency matrix by repeated boolean squaring. shares no machinery with
## add.ancestors()/remove.ancestors(), which is what makes it a usable
## reference for them.

#' @importFrom graph nodes edgeMatrix
.anc_closure <- function(g) {
    v <- nodes(g)
    p <- length(v)
    A <- matrix(FALSE, p, p)
    em <- edgeMatrix(g)
    if (ncol(em) > 0)
        A[cbind(em["from", ], em["to", ])] <- TRUE
    R <- A
    repeat {
        N <- R | ((R %*% A) > 0)
        if (identical(N, R))
            break
        R <- N
    }

    R
}

##
## REPEATED COVERED ARC REVERSAL ALGORITHM (Castelo and Kocka, JMLR, 2003)
## ADAPTED TO INTERVENTIONS IN (Castelo, 2025)
##

## build a logical mask indicated what edges are "covered" in the input DAG
## An arc is target-protected iff some SINGLE target contains exactly one of
## its endpoints (Hauser and Buehlmann 2012), which is what .separates()
## tests. The union of the targets is not a substitute: it agrees only when
## every target is a singleton. For targets = list(integer(0), c(1L, 2L))
## both ends of 1 -> 2 lie in the union, yet no target separates them, so the
## arc is I-covered -- reversing it stays inside the I-equivalence class, and
## the union test wrongly declared it protected. That made rcar() unable to
## walk those arcs and left ncr.nh() exposing a within-class reversal.
##
## `targets` is the family, a list of integer vectors of vertex indices.
## Returns a logical vector, one per (from, to) pair given.
.separates <- function(targets, from, to) {
    if (length(targets) == 0L || length(from) == 0L)
        return(rep(FALSE, length(from)))
    sep <- rep(FALSE, length(from))
    for (I in targets)
        sep <- sep | xor(from %in% I, to %in% I)
    sep
}

## `targets` is the target family; a covered arc stays out of the covered set
## only when some target separates its endpoints

#' @importFrom graph nodes edgeMatrix
cedges <- function(dag, targets) {
    v <- nodes(dag)
    em <- edgeMatrix(dag)
    ## an edgeless DAG has no covered arcs. the early return is needed
    ## because mapply() over zero-length inputs returns list(), not
    ## logical(0), and the '&' below would then fail on a list. rcar(), the
    ## only caller, guards numEdges(dag) == 0 before it gets here, so this
    ## never fired in the search -- but the function should be total, and
    ## the regression tests call it directly.
    if (ncol(em) == 0L)
        return(logical(0))
    pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
    cemask <- mapply(function(pafrom, pato, from) identical(sort(pafrom), sort(setdiff(pato, from))),
                     pasets[em["from", ]], pasets[em["to", ]], v[em["from", ]])
    temask <- .separates(targets, em["from", ], em["to", ])
    cemask & !temask
}

## resample helper function
resample <- function(x, ...) x[sample.int(length(x), ...)]

## one draw of R's own R_unif_index(), the function sample.int() uses
## internally, exposed from C. rcar() is the only RNG consumer in the
## search, and a C port of it has to reproduce R's stream draw for draw --
## which means calling R_unif_index() rather than scaling a uniform, since
## under the default sample.kind="Rejection" the number of unif_rand() calls
## per draw is data dependent (1 to 3 for a size-1 draw). This wrapper
## exists so tests/test_rng_equivalence.R can pin the three identities the
## port relies on, on both the value and the resulting .Random.seed:
##
##   sample.int(n, 1)      == .unif_index(n) + 1
##   sample(0:r, size=1)   == .unif_index(length(0:r))
##   resample(x, size=1)   == x[.unif_index(length(x)) + 1]
##
## 'n' is the population size, and the result lies in 0:(n-1) -- 0-based, as
## C wants it. See src/rng.c.
.unif_index <- function(n) .Call(C_unif_index, as.double(n))

## RCAR: repeated covered arc reversal algorithm
## `targets` is the target family (see .separates())
## returns a list(dag=, anc=, pasets=) since every reversal it performs,
## although always cycle-safe by construction (a covered edge cannot
## introduce a cycle), still changes true ancestor relationships and parent
## sets, and must keep both 'anc' and 'pasets' in sync for subsequent
## neighborhood generation and scoring to remain correct. 'vidx' is a named
## integer vector mapping vertex name -> dat-column index (as used by
## 'pasets'), built once per search by the caller.

#' @importFrom graph removeEdge addEdge numEdges edgeMatrix nodes
rcar <- function(dag, r, targets, anc, pasets, vidx) {
    if (numEdges(dag) == 0)
        return(list(dag=dag, anc=anc, pasets=pasets))
    cemask <- cedges(dag, targets)
    if (!any(cemask))
        return(list(dag=dag, anc=anc, pasets=pasets))

    tmp.g <- dag
    v <- nodes(tmp.g)
    rr <- sample(0:r, size=1)
    for (i in seq_len(rr)) {
        em <- edgeMatrix(tmp.g)
        cemask <- cedges(tmp.g, targets)
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
#' their sum. It may then be omitted, in which case [`population`]'s own `n`
#' and `C` supply the size and how it is split between observational and
#' interventional environments.
#'
#' @param scorefun (Default is [`iBIC`]) A function to calculate the goodness
#' of fit (GoF) score of a DAG on a given data set.
#'
#' @param verbose (Default TRUE) Show progress in the calculations.
#'
#' @param engine (Default `"C"`) A character string selecting the search
#' engine: `"C"` (default) maintains the DAG, its ancestor relation and its
#' parent sets in compiled code; `"R"` uses the pure-R implementation and is
#' provided for testing and verification. Both follow the same trajectory and
#' return the same result. `"C"` requires a `scorefun` able to score a whole
#' neighbourhood at once, which [`iBIC`] and [`iBGe`] are; with any other
#' score function the `"R"` engine is used regardless.
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
hillclimbing <- function(x, targets=list(integer(0)),
                         target.index=NULL,  scorefun=iBIC,
                         verbose=TRUE, engine=c("C", "R")) {

    engine <- match.arg(engine)

    x <- .check_input_data(x)
    dag <- graphNEL(colnames(x), edgemode="directed")
    attr(x, "sanitycheck") <- TRUE

    targets <- .check_targets(targets, ncol(x))
    target.index <- .resolve.target.index(x, targets, target.index)
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

    s0 <- -Inf
    s1 <- scorefun(dag, x, targets=targets, target.index=target.index,
                   cached.scores=cached.scores, global.sufstats=global.sufstats)

    if (verbose) {
        msg <- "Running a straightforward hill-climbing algorithm"
        if (!is.null(scorefun.name))
            msg <- paste(msg, "with the {scorefun.name} score function")
        cli_progress_bar(msg)
        cli_progress_step("Score {s1}", spinner=TRUE)
    }

    ## the C engine keeps the DAG, its ancestor relation and its parent sets
    ## in compiled state behind an external pointer, and enumerates the
    ## neighbourhood there too; the R engine is the reference implementation
    ## the differential tests score it against. Scoring is identical in both:
    ## it goes through score.nh() either way.
    ##
    ## The C engine needs a scorefun that can score a whole neighbourhood in
    ## one call, because it holds no graphNEL during the search and so cannot
    ## serve score.nh()'s per-candidate fallback. iBIC() and iBGe() qualify;
    ## anything else falls back to R rather than failing.
    if (use.c) {
        st <- .Call(C_dag_new, ncol(x))
        while (s1 > s0) {
            s0 <- s1
            ne <- .Call(C_dag_nh, st, 2L, list())       ## 2 = ar
            pasets <- .Call(C_dag_pasets, st)
            ## only the winner is needed, so the O(p) exact summation is
            ## paid for the provable handful of candidates that could still
            ## be the maximum rather than for all O(p^2) of them
            am <- nh.argmax.fun(ne$op, vidx.nodes[ne$u], vidx.nodes[ne$v],
                                pasets, global.sufstats, cached.scores,
                                verify.band, .Call(C_dag_pastamp, st))
            b <- am$index
            .Call(C_dag_apply_move, st, ne$op[b], ne$u[b], ne$v[b])
            s1 <- am$total

            .debug_assertions(st, NULL, NULL, NULL, x, vnames)

            if (verbose)
                cli_progress_update()
        }
        dag <- .graphNEL_from_edgeM(vnames, .Call(C_dag_edgeM, st))
    } else {
        while (s1 > s0) {
            s0 <- s1
            ne <- ar.nh(dag, anc)
            sco <- score.nh(ne, dag, x, targets, target.index, cached.scores,
                            global.sufstats, pasets, vidx.nodes,
                            supports.pasets, scorefun, nh.scores.fun)
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
            pasets <- move.pasets(pasets, b.op, vidx.nodes[b.u],
                                  vidx.nodes[b.v])
            dag <- apply.move(dag, b.op, b.u, b.v, vnames)
            s1 <- sco[b]

            .debug_assertions(NULL, dag, anc, pasets, x, vnames)

            if (verbose)
                cli_progress_update()
        }
    }

    if (verbose)
        cli_progress_done("straightforward hill-climbing algorithm completed")

    list(dag=dag, sco=s1)
}
