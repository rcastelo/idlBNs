##
## NEIGHBORHOODS
##

## NR: non-reversals neighborhood (addition and removal only)

#' @importFrom graph edgeL removeEdge addEdge nodes
#' @importFrom RBGL dag.sp
nr.nh <- function(dag) {
  v <- nodes(dag)
  e <- edgeL(dag)
  nr <- list()
  nr.i <- 0
  for (i in seq_along(e)) {
    a <- v[e[[i]]$edges]
    na <- setdiff(v, a)
    for (j in seq_along(na)) { ## go through non-adjacent vertices
      d <- dag.sp(dag, na[j])$distance[names(e)[i]]
      if (is.nan(d) || is.infinite(d)) {
        tmp.g <- addEdge(names(e)[i], na[j], dag)
        nr.i <- nr.i + 1
        nr[[nr.i]] <- tmp.g
      }
    }
    for (j in seq_along(a)) { ## go through adjacent vertices
      tmp.g <- removeEdge(names(e)[i], a[j], dag)
      nr.i <- nr.i + 1
      nr[[nr.i]] <- tmp.g
    }
  }
  nr
}

## AR: all-reversals neighborhood (NR + all-arc-reversals)

#' @importFrom graph edgeL removeEdge addEdge nodes
#' @importFrom RBGL dag.sp
ar.nh <- function(dag) {
  v <- nodes(dag)
  e <- edgeL(dag)
  ar <- nr.nh(dag)
  ar.i <- length(ar)
  for (i in seq_along(e)) { ## reverse edges
    a <- v[e[[i]]$edges]
    for (j in seq_along(a)) { ## go through adjacent vertices
      tmp.g <- removeEdge(names(e)[i], a[j], dag)
      d <- dag.sp(tmp.g, names(e)[i])$distance[a[j]]
      if (is.nan(d) || is.infinite(d)) {
        tmp.g <- addEdge(a[j], names(e)[i], tmp.g)
        ar.i <- ar.i + 1
        ar[[ar.i]] <- tmp.g
      }
    }
  }
  ar
}

## NCR: non-covered arc reversals neighborhood (NR + non-covered-arc-reversals)

#' @importFrom graph edgeL removeEdge addEdge edgeMatrix nodes
#' @importFrom RBGL dag.sp
ncr.nh <- function(dag, utargets=integer(0)) {
  v <- nodes(dag)
  e <- edgeL(dag)
  ncr <- nr.nh(dag)
  ncr.i <- length(ncr)
  em <- edgeMatrix(dag)
  pasets <- split(v[em["from", ]], factor(v[em["to", ]], levels=v))
  for (i in seq_along(e)) { ## reverse edges
    a <- v[e[[i]]$edges]
    for (j in seq_along(a)) { ## go through adjacent vertices
      ced <- identical(sort(pasets[[names(e)[i]]]), sort(setdiff(pasets[[a[j]]], names(e)[i])))
      if (!ced || any(c(e[[i]]$edges[j], i) %in% utargets)) { ## NCR including not interventionally covered
        tmp.g <- removeEdge(names(e)[i], a[j], dag)
        d <- dag.sp(tmp.g, names(e)[i])$distance[a[j]]
        if (is.nan(d) || is.infinite(d)) {
          tmp.g <- addEdge(a[j], names(e)[i], tmp.g)
          ncr.i <- ncr.i + 1
          ncr[[ncr.i]] <- tmp.g
        }
      }
    }
  }
  ncr
}

##
## REPEATED COVERED ARC REVERSAL ALGORITHM
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

#' @importFrom graph removeEdge addEdge numEdges edgeMatrix nodes
rcar <- function(dag, r, utargets) {
  if (numEdges(dag) == 0)
    return(dag)
  cemask <- cedges(dag, utargets)
  if (!any(cemask))
    return(dag)

  tmp.g <- dag
  v <- nodes(tmp.g)
  rr <- sample(0:r, size=1)
  for (i in seq_len(rr)) {
    em <- edgeMatrix(tmp.g)
    cemask <- cedges(tmp.g, utargets)
    rndce <- resample(which(cemask), size=1)
    tmp.g <- removeEdge(v[em["from", rndce]], v[em["to", rndce]], tmp.g)
    tmp.g <- addEdge(v[em["to", rndce]], v[em["from", rndce]], tmp.g) ## a covered edge cannot introduce a cycle
  }
  tmp.g
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

  s0 <- -Inf
  s1 <- scorefun(g=dag, dat=dat, targets=targets, target.index=target.index,
                 cached.scores=cached.scores, global.sufstats=global.sufstats)

  if (verbose)
    cli_progress_step("Score {s1}")

  while (s1 > s0) {
    s0 <- s1
    ne <- ar.nh(dag)
    s1 <- sapply(ne, function(g, d, tgts, tgt.idx, chd.sco)
                       scorefun(g=g, dat=d, targets=tgts, target.index=tgt.idx,
                                cached.scores=chd.sco), dat, targets,
                                target.index, cached.scores)
    dag <- ne[[which.max(s1)]]
    s1 <- max(s1)

    if (verbose)
      cli_progress_update()
  }
  list(dag=dag, sco=s1)
}
