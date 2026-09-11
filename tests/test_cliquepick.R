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
    ## Same draws, same keys, computed without rebuilding a matrix each time.
    ## ord2amo() recomputed the edge index and allocated an m x m matrix per
    ## draw, and that R overhead -- not the sampling -- was 4.7 s of this
    ## file's 6.3. The orientation of an edge is decided entirely by which
    ## endpoint comes first in the order, so the key is a vector expression
    ## over a precomputed edge list. Nothing about the TEST changes: the same
    ## 400 draws per class member, over the same eight components, compared
    ## against the same brute-force key set. Trimming the draws instead would
    ## have cost real power -- at 150 per member a sampler biased 3% toward
    ## one clique goes undetected, where at 400 it is caught.
    fastkey <- local({
        m <- nrow(U); ei <- which(U & upper.tri(U), arr.ind = TRUE)
        ia <- ei[, 1]; ib <- ei[, 2]
        function(ord) {
            pos <- integer(m); pos[ord] <- seq_len(m)
            first <- pos[ia] < pos[ib]
            idx <- ifelse(first, ia + (ib - 1L) * m, ib + (ia - 1L) * m)
            paste0("o", paste(sort(idx), collapse = ","))
        }
    })
    stopifnot(identical(fastkey(seq_len(nrow(U))),          # same as ord2amo
                        amokey(ord2amo(U, seq_len(nrow(U))))))
    draws <- replicate(400 * n, fastkey(cpsamp(U)))
    tab <- table(factor(draws, levels = keys))
    stopifnot(all(tab > 0), chisq.test(tab)$p.value > 1e-4)
}
stopifnot(tested > 0)

cat("test_cliquepick.R: OK\n")

## Listing is proportional to the output, not quadratic in it. cp_unrank()
## decodes a rank by counting allowed completions and skipping whole prefix
## subtrees; walking to the rank-th permutation instead made a listing cost
## O(c^2), which for a complete component of 8 vertices was 7.1 s for its
## 40320 members (0.18 ms each, against 0.003 ms at m = 6). The guard below
## is deliberately loose -- it is a complexity check, not a timing pin.
local({
    Km <- function(m) { U <- matrix(TRUE, m, m); diag(U) <- FALSE; U }
    tm <- function(m) { U <- Km(m)
        min(replicate(3, system.time(
            .Call(idlBNs:::C_cp_amo_list, U, Inf))[["elapsed"]])) }
    ## |class| grows 9-fold from m = 8 to m = 9 (8! to 9!). Linear listing
    ## costs about 9x, quadratic about 81x; m = 6 and 7 are below the timer's
    ## resolution, so the two measurable sizes are the ones compared.
    t8 <- tm(8L); t9 <- tm(9L)
    stopifnot(t8 > 0, t9 < 25 * t8)
})
cat("test_cliquepick.R: listing is output-proportional\n")

## The clique is chosen with an EXACT uniform integer, and the counts must be
## exact for that to mean anything.
##
## unif_rand() lives on a ~2^-32 grid under Mersenne-Twister, so the old
## `unif_rand() * total` spread at most 2^32 distinct values over [0, total):
## at total = 2^40 the top 256 units of the range were unreachable, and any
## branch narrower than total * 2^-32 could hold no grid point at all,
## however positive its weight. A 13-vertex clique already has 13! > 2^32.
## R_unif_index() draws bits and rejects instead, so every branch gets its
## share exactly. Above 2^53 -- or a clique past 18!, where cp_fac() stops
## being exact -- the weights themselves round, and the samplers decline
## rather than return a draw they cannot vouch for. COUNTING still answers.
local({
    K <- function(m) { A <- matrix(TRUE, m, m); diag(A) <- FALSE; A }
    ## exact side of both boundaries: draws, and the count is m!
    for (m in c(13L, 18L)) {
        stopifnot(isTRUE(all.equal(.Call(idlBNs:::C_cp_amo_count, K(m)),
                                   factorial(m))),
                  !is.null(.Call(idlBNs:::C_cp_amo_sample, K(m))))
    }
    ## past it: declines, but counting is unaffected -- an approximate double
    ## is the expected answer for a class of 10^18 members
    for (m in c(20L, 25L)) {
        stopifnot(isTRUE(all.equal(.Call(idlBNs:::C_cp_amo_count, K(m)),
                                   factorial(m))),
                  is.null(.Call(idlBNs:::C_cp_amo_sample, K(m))),
                  is.null(.Call(idlBNs:::C_cp_amo_list, K(m), Inf)))
    }
    ## and the draw is uniform where a clique choice is actually made: this
    ## component has several maximal cliques, so cp_sample() weights them
    U <- matrix(FALSE, 6, 6)
    for (e in list(c(1,2), c(2,3), c(1,3), c(3,4), c(4,5), c(4,6), c(5,6)))
        U[e[1], e[2]] <- U[e[2], e[1]] <- TRUE
    n <- .Call(idlBNs:::C_cp_amo_count, U)
    set.seed(1)
    ## 120000 draws, keyed by a base-7 encoding of the order rather than by
    ## paste(collapse = ","): the same draws and the same categories, without
    ## the per-draw string build. The count is not negotiable -- against a
    ## sampler biased 3% toward one clique, 120000 draws catch it and 30000
    ## do not -- so the saving has to come from the cost per draw.
    w <- 7^(seq_len(6L) - 1L)
    tab <- table(replicate(120000,
                 sum(.Call(idlBNs:::C_cp_amo_sample, U) * w)))
    stopifnot(length(tab) == n,
              stats::chisq.test(as.vector(tab))$p.value > 1e-4)
})
cat("test_cliquepick.R: exact weighted choice, and declines when counts round\n")
