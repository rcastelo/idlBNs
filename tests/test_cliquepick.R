## 2026-09-10 tests for the Clique-Picking implementation in src/cliquepick.c,
## which counts and uniformly samples the acyclic moral orientations of a
## chordal graph in polynomial time (Wienoebst, Bannach & Liskiewicz, JMLR 24,
## 2023) and replaces the enumerate-and-deduplicate routine that was
## Theta(m! * #AMO).
##
## Everything is checked against independent references:
##   1. the count, against brute-force enumeration over vertex orderings;
##   2. the count, against two closed forms -- a complete graph on m has m!
##      acyclic moral orientations and a tree on m has m;
##   3. the sampler, by chi-square against the enumerated class.
##
## Test graphs are the undirected chain components of essential graphs of random
## DAGs: chordal by construction, and exactly the input the search feeds it.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

cpcount <- function(U) .Call(idlBNs:::C_cp_amo_count, U)
cpsamp  <- function(U) .Call(idlBNs:::C_cp_amo_sample, U)

perms <- function(v) if (length(v) <= 1) list(v) else
    unlist(lapply(seq_along(v), function(i)
        lapply(perms(v[-i]), function(r) c(v[i], r))), recursive = FALSE)

ord2amo <- function(U, ord) {
    m <- nrow(U); pos <- integer(m); pos[ord] <- seq_len(m)
    ei <- which(U & upper.tri(U), arr.ind = TRUE); A <- matrix(FALSE, m, m)
    if (nrow(ei))
        A[cbind(ifelse(pos[ei[,1]] < pos[ei[,2]], ei[,1], ei[,2]),
                ifelse(pos[ei[,1]] < pos[ei[,2]], ei[,2], ei[,1]))] <- TRUE
    A
}
amokey <- function(A) paste0("o", paste(which(A), collapse = ","))

## brute force: every vertex ordering, keep the moral ones, de-duplicate
brute <- function(U) {
    m <- nrow(U)
    seen <- new.env(parent = emptyenv()); n <- 0L
    for (perm in perms(seq_len(m))) {
        A <- ord2amo(U, perm)
        bad <- FALSE
        for (v in seq_len(m)) { pa <- which(A[, v])
            if (length(pa) > 1) for (a in seq_along(pa)) for (b in seq_len(a - 1))
                if (!A[pa[a], pa[b]] && !A[pa[b], pa[a]]) { bad <- TRUE; break } }
        if (bad) next
        k <- amokey(A)
        if (is.null(seen[[k]])) { assign(k, TRUE, envir = seen); n <- n + 1L }
    }
    list(n = n, keys = ls(seen))
}

## chain components of essential graphs of random DAGs
pool <- function(p, prob, howmany, lo = 2L, hi = 7L, seed) {
    set.seed(seed); out <- list()
    while (length(out) < howmany) {
        o <- sample(p); A <- matrix(FALSE, p, p)
        for (i in 1:(p-1)) for (j in (i+1):p) if (runif(1) < prob) A[o[i], o[j]] <- TRUE
        nn <- as.character(seq_len(p)); dimnames(A) <- list(nn, nn)
        am <- as(as(dag2essgraph(as(as(A * 1, "graphAM"), "graphNEL"),
                                 targets = list(integer(0))), "graphNEL"), "matrix") > 0
        un <- am & t(am); dimnames(un) <- NULL
        lab <- integer(p); nl <- 0
        for (v in seq_len(p)) {
            if (lab[v]) next
            nl <- nl + 1; st <- v; lab[v] <- nl
            while (length(st)) { w <- st[1]; st <- st[-1]
                nb <- which(un[w, ] & lab == 0); lab[nb] <- nl; st <- c(st, nb) }
        }
        for (k in seq_len(nl)) { v <- which(lab == k)
            if (length(v) >= lo && length(v) <= hi)
                out[[length(out) + 1L]] <- un[v, v, drop = FALSE] }
    }
    out
}

## 1. the count, against brute force ------------------------------------------
P <- pool(9, 0.45, 250, seed = 1)
for (U in P) {
    b <- brute(U)
    stopifnot(isTRUE(all.equal(as.numeric(b$n), cpcount(U))))
}

## 2. closed forms -------------------------------------------------------------
set.seed(2)
for (m in 2:9) {
    Kc <- matrix(TRUE, m, m); diag(Kc) <- FALSE
    stopifnot(isTRUE(all.equal(cpcount(Kc), factorial(m))))
    Tr <- matrix(FALSE, m, m)
    for (i in 2:m) { j <- sample(i - 1, 1); Tr[i, j] <- TRUE; Tr[j, i] <- TRUE }
    stopifnot(isTRUE(all.equal(cpcount(Tr), as.numeric(m))))
}
## no size limit: a complete graph on 20 vertices
K20 <- matrix(TRUE, 20, 20); diag(K20) <- FALSE
stopifnot(isTRUE(all.equal(cpcount(K20), factorial(20))))

## 3. uniformity ---------------------------------------------------------------
set.seed(3)
tested <- 0
for (U in pool(9, 0.5, 60, lo = 3L, hi = 6L, seed = 4)) {
    n <- cpcount(U)
    if (n < 4 || n > 60) next
    tested <- tested + 1
    if (tested > 8) break
    keys <- brute(U)$keys
    draws <- replicate(400 * n, amokey(ord2amo(U, cpsamp(U))))
    tab <- table(factor(draws, levels = keys))
    stopifnot(all(tab > 0), chisq.test(tab)$p.value > 1e-4)
}
stopifnot(tested > 0)

cat("test_cliquepick.R: OK\n")
