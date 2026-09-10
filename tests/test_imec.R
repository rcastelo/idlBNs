## 2026-09-10 tests for the exact I-equivalence class machinery of R/imec.R,
## which replaces the rcar() random walk by uniform sampling on the class.
##
## Three things are pinned here, each against an INDEPENDENT reference rather
## than against imec.R itself:
##
##   1. .iessgraph() against pcalg::dag2essgraph(), which is a different
##      implementation of the same construction (pcalg's C++
##      EssentialGraph::replaceUnprotected).
##   2. imec.list()/imec.size() against brute force over every acyclic
##      orientation of the skeleton, filtered by Hauser & Buehlmann's (2012)
##      I-Markov equivalence criterion -- same skeleton and immoralities, and
##      the same skeleton in every intervention graph.
##   3. uniformity of imec.sample() by chi-square against that brute-force class.
##
## Also pinned: reversing an I-covered arc leaves the I-essential graph
## unchanged, which is what makes rcar()'s moves score-neutral in the first
## place, and is the assumption the fallback path still relies on.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

iess   <- idlBNs:::.iessgraph
gnel   <- function(A) { nn <- as.character(seq_len(nrow(A)))
    dimnames(A) <- list(nn, nn); as(as(A * 1, "graphAM"), "graphNEL") }
rdag   <- function(p, prob) { o <- sample(p); A <- matrix(FALSE, p, p)
    for (i in 1:(p-1)) for (j in (i+1):p) if (runif(1) < prob) A[o[i], o[j]] <- TRUE; A }
## Target families. Single-vertex interventions are the easy case: whenever
## either endpoint of an arc is targeted the two have different membership
## masks, so the arc is target-protected. A target containing BOTH endpoints
## gives them equal masks and leaves the arc unprotected -- a different branch
## of .protects() / tm_protects() that single-vertex families never reach. So
## draw some targets as the two ends of an actual arc, some as random sets of
## two or three, and some as single vertices.
rtgts  <- function(p, k, A = NULL) {
    if (k == 0) return(list(integer(0)))
    arcs <- if (is.null(A)) NULL else which(A, arr.ind = TRUE)
    one <- function(i) {
        u <- runif(1)
        if (!is.null(arcs) && nrow(arcs) > 0 && u < 0.40)
            as.integer(sort(arcs[sample(nrow(arcs), 1), ]))
        else if (u < 0.70)
            as.integer(sort(sample(p, min(p, sample(2:3, 1)))))
        else
            as.integer(sample(p, 1))
    }
    c(list(integer(0)), lapply(seq_len(k), one))
}
key    <- function(A) paste(which(A), collapse = ",")
acyc   <- function(A) { a <- A; repeat { b <- a | ((a %*% a) > 0)
    if (all(b == a)) break; a <- b }; !any(diag(a)) }
skel   <- function(A) A | t(A)
immor  <- function(A) { p <- nrow(A); out <- character(0)
    for (v in seq_len(p)) { pa <- which(A[, v])
        if (length(pa) > 1) for (i in seq_along(pa)) for (j in seq_len(i - 1))
            if (!A[pa[i], pa[j]] && !A[pa[j], pa[i]])
                out <- c(out, sprintf("%d-%d-%d", min(pa[i],pa[j]), v, max(pa[i],pa[j]))) }
    sort(unique(out)) }
subI   <- function(A, I) { B <- A; if (length(I)) B[, I] <- FALSE; B }

brute <- function(A0, tg) {
    p <- nrow(A0); E <- which(skel(A0) & upper.tri(skel(A0)), arr.ind = TRUE)
    m <- nrow(E); ref.im <- immor(A0); ref.sk <- lapply(tg, function(I) skel(subI(A0, I)))
    out <- list()
    for (b in 0:(2^m - 1)) {
        A <- matrix(FALSE, p, p)
        for (e in seq_len(m)) if (bitwAnd(bitwShiftR(b, e - 1), 1L))
            A[E[e,1], E[e,2]] <- TRUE else A[E[e,2], E[e,1]] <- TRUE
        if (!acyc(A) || !identical(immor(A), ref.im)) next
        if (all(vapply(seq_along(tg), function(t)
                all(skel(subI(A, tg[[t]])) == ref.sk[[t]]), TRUE)))
            out[[key(A)]] <- A
    }
    out
}

## 1. the I-essential graph, against pcalg -------------------------------------
set.seed(1)
for (it in 1:150) {
    p <- sample(4:8, 1); A <- rdag(p, 0.3); tg <- rtgts(p, sample(0:3, 1), A)
    ref <- as(as(dag2essgraph(gnel(A), targets = tg), "graphNEL"), "matrix") > 0
    dimnames(ref) <- NULL
    stopifnot(all(iess(A, tg) == ref))
}

## 2. the class itself, against brute force ------------------------------------
set.seed(2)
for (it in 1:80) {
    p <- sample(5:7, 1); A <- rdag(p, 0.35); tg <- rtgts(p, sample(0:2, 1), A)
    L <- idlBNs:::imec.list.ref(A, tg); if (is.null(L)) next
    B <- brute(A, tg)
    stopifnot(setequal(names(B), vapply(L, key, "")),
              isTRUE(all.equal(idlBNs:::imec.size(A, tg), length(B))))
}

## 3. uniformity ---------------------------------------------------------------
set.seed(3)
found <- 0
while (found < 3) {
    p <- sample(5:7, 1); A <- rdag(p, 0.45); tg <- rtgts(p, 1, A)
    cs <- idlBNs:::imec.size(A, tg)
    if (is.na(cs) || cs < 5) next
    found <- found + 1
    draws <- replicate(6000, key(idlBNs:::imec.sample(A, tg)))
    tab <- table(factor(draws, levels = vapply(idlBNs:::imec.list(A, tg), key, "")))
    stopifnot(all(tab > 0), chisq.test(tab)$p.value > 1e-4)
}

## 4. an I-covered arc reversal does not leave the I-equivalence class ---------
##
## I-covered means covered AND unseparated: some single target must contain
## exactly one endpoint for the arc to be target-protected. Skipping every
## arc with an endpoint in the UNION of the targets, as this once did, threw
## away precisely the arcs a multi-vertex target makes I-covered -- both ends
## of 1 -> 2 lie in the union of list(integer(0), c(1L, 2L)), yet no target
## separates them. Those are the arcs rcar() walks, so they are the ones this
## invariant has to hold for. cedges(), which the search actually consults,
## is checked against the same definition.
separates <- function(tg, a, b)
    any(vapply(tg, function(I) xor(a %in% I, b %in% I), NA))
set.seed(4)
nicov <- 0L
for (it in 1:100) {
    p <- sample(5:8, 1); A <- rdag(p, 0.35); tg <- rtgts(p, sample(1:2, 1), A)
    E <- which(A, arr.ind = TRUE)
    ref <- iess(A, tg)
    ## the search's own covered-arc mask, in edgeMatrix order
    g <- gnel(A); em <- graph::edgeMatrix(g)
    mR <- unname(idlBNs:::cedges(g, tg))
    for (e in seq_len(nrow(E))) {
        a <- E[e, 1]; b <- E[e, 2]
        icov <- !separates(tg, a, b) &&
                setequal(which(A[, b]), union(a, which(A[, a])))
        ## cedges() must agree, arc for arc
        col <- which(as.integer(em["from", ]) == a & as.integer(em["to", ]) == b)
        stopifnot(length(col) == 1L, identical(mR[col], icov))
        if (!icov) next
        nicov <- nicov + 1L
        B <- A; B[a, b] <- FALSE; B[b, a] <- TRUE
        stopifnot(all(iess(B, tg) == ref))
    }
}
stopifnot(nicov > 0L)      ## anti-vacuity

cat("test_imec.R: OK\n")
