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
rtgts  <- function(p, k) if (k == 0) list(integer(0)) else
    c(list(integer(0)), lapply(sample(p, k), function(v) as.integer(v)))
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
    p <- sample(4:8, 1); A <- rdag(p, 0.3); tg <- rtgts(p, sample(0:3, 1))
    ref <- as(as(dag2essgraph(gnel(A), targets = tg), "graphNEL"), "matrix") > 0
    dimnames(ref) <- NULL
    stopifnot(all(iess(A, tg) == ref))
}

## 2. the class itself, against brute force ------------------------------------
set.seed(2)
for (it in 1:80) {
    p <- sample(5:7, 1); A <- rdag(p, 0.35); tg <- rtgts(p, sample(0:2, 1))
    L <- idlBNs:::imec.list(A, tg); if (is.null(L)) next
    B <- brute(A, tg)
    stopifnot(setequal(names(B), vapply(L, key, "")),
              isTRUE(all.equal(idlBNs:::imec.size(A, tg), length(B))))
}

## 3. uniformity ---------------------------------------------------------------
set.seed(3)
found <- 0
while (found < 3) {
    p <- sample(5:7, 1); A <- rdag(p, 0.45); tg <- rtgts(p, 1)
    dec <- idlBNs:::imec.decompose(A, tg)
    cs <- idlBNs:::imec.size(A, tg, dec = dec)
    if (is.na(cs) || cs < 5) next
    found <- found + 1
    draws <- replicate(6000, key(idlBNs:::imec.sample(A, tg, dec = dec)))
    tab <- table(factor(draws, levels = vapply(idlBNs:::imec.list(A, tg, dec = dec), key, "")))
    stopifnot(all(tab > 0), chisq.test(tab)$p.value > 1e-4)
}

## 4. an I-covered arc reversal does not leave the I-equivalence class ---------
set.seed(4)
for (it in 1:100) {
    p <- sample(5:8, 1); A <- rdag(p, 0.35); tg <- rtgts(p, sample(1:2, 1))
    ut <- sort(unique(unlist(tg))); E <- which(A, arr.ind = TRUE)
    ref <- iess(A, tg)
    for (e in seq_len(nrow(E))) {
        a <- E[e, 1]; b <- E[e, 2]
        if (a %in% ut || b %in% ut) next
        if (!setequal(which(A[, b]), union(a, which(A[, a])))) next   # not covered
        B <- A; B[a, b] <- FALSE; B[b, a] <- TRUE
        stopifnot(all(iess(B, tg) == ref))
    }
}

cat("test_imec.R: OK\n")
