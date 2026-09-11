## 2026-09-10 differential test for the C port of the exact I-MEC machinery in
## src/imec.c against the R reference in R/imec.R.
##
## Sections, the first three in increasing strength:
##
##   1. |[D]_I| and the enumeration agree, as SETS.
##   2. The enumeration is duplicate free and exactly of size |[D]_I|, so it is
##      a bijection onto the class. Order is NOT pinned: the R reference walks
##      permutations of each chain component while the C engine unranks the
##      Clique-Picking decomposition, and the sampler draws its own clique
##      rather than indexing into the enumeration.
##   3. The sampler agrees bit for bit AND leaves .Random.seed in the same
##      place. The stream is the load-bearing part, but the invariant is
##      EQUIVALENCE, not a draw count: Clique-Picking spends one unif_rand()
##      on the clique choice (none when the node has a single clique), then
##      R_unif_index() per Fisher-Yates step -- 1 to 3 unif_rand() calls each
##      under sample.kind = "Rejection" -- retried for as long as the
##      forbidden-prefix test rejects the permutation, and recurses into every
##      subproblem. The total varies with the component AND with the seed:
##      over 400 seeds a complete component of 5 vertices consumed between 4
##      and 13 unif_rand() calls, an 8-vertex path between 2 and 7. What is
##      pinned is that both engines reach the same cp_sample() on the same
##      components in the same order, so they consume the same draws and end
##      in the same place. (The identity sample.int(n, 1) == R_unif_index(n)
##      + 1 belongs to rcar(); see tests/test_rng_equivalence.R.)
##   4. Both engines run end to end and agree on score -- bitwise -- DAG and
##      final seed, over every sampler/escape combination.
##   5. C_dag_set_edges rejects a malformed edge list without mutating state.
##   6. Exact sampling has no whole-graph size limit; the 64-bit vertex sets
##      are per chain component.
##   7. A declined draw consumes no randomness at all.

suppressPackageStartupMessages({
  library(graph)
  library(idlBNs)
})

dag_new    <- function(p) .Call(idlBNs:::C_dag_new, as.integer(p))
dag_set    <- function(st, em) .Call(idlBNs:::C_dag_set_edges, st, em)
dag_edgeM  <- function(st) .Call(idlBNs:::C_dag_edgeM, st)
imec_size  <- function(st, tg, m) .Call(idlBNs:::C_dag_imec_size, st, tg)
imec_memb  <- function(st, tg, m, mx = Inf) .Call(idlBNs:::C_dag_imec_members, st, tg, as.double(mx))
imec_samp  <- function(st, tg, m) .Call(idlBNs:::C_dag_imec_sample, st, tg)

adj2em <- function(A) { E <- which(A, arr.ind = TRUE)
    matrix(as.integer(c(rbind(E[, 1], E[, 2]))), nrow = 2,
           dimnames = list(c("from", "to"), NULL)) }
em2adj <- function(em, p) { A <- matrix(FALSE, p, p)
    if (ncol(em)) A[cbind(em["from", ], em["to", ])] <- TRUE; A }
key    <- function(A) paste(which(A), collapse = ",")
## Two learned DAGs compared by adjacency matrix. Sorting the flattened edge
## matrix, as this once did, discards the pairing of each arc's tail with its
## head: 1->2, 3->4 and 1->4, 3->2 both flatten and sort to (1,2,3,4), so two
## genuinely different DAGs compared equal.
dagadj <- function(g) { m <- as(g, "matrix"); m[] <- m != 0; m }
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

set.seed(7)
for (it in 1:150) {
    p <- sample(4:8, 1); A <- rdag(p, 0.4); tg <- rtgts(p, sample(0:2, 1), A)
    st <- dag_new(p); dag_set(st, adj2em(A))
    stopifnot(identical(em2adj(dag_edgeM(st), p), A))     # set_edges round trip

    ## 1 & 2: size and enumeration
    sz.R <- idlBNs:::imec.size(A, tg); sz.C <- imec_size(st, tg, 8L)
    stopifnot(isTRUE(all.equal(sz.R, sz.C)))
    L.R <- idlBNs:::imec.list.ref(A, tg); L.C <- imec_memb(st, tg, 8L)
    if (is.null(L.R)) { stopifnot(is.null(L.C)); next }
    k.R <- vapply(L.R, key, ""); k.C <- vapply(L.C, function(em) key(em2adj(em, p)), "")
    stopifnot(length(k.C) == length(k.R),
              !anyDuplicated(k.C), !anyDuplicated(k.R),
              identical(sort(k.C), sort(k.R)),
              isTRUE(all.equal(length(k.C), sz.C)))
    ## listing must leave the state untouched
    stopifnot(identical(em2adj(dag_edgeM(st), p), A))

    ## 3: the sampler, values and RNG stream
    set.seed(1000 + it)
    B.R <- idlBNs:::imec.sample(A, tg)
    seed.R <- .Random.seed
    set.seed(1000 + it)
    ok <- imec_samp(st, tg, 8L)
    seed.C <- .Random.seed
    stopifnot(isTRUE(ok), identical(seed.R, seed.C))
    stopifnot(identical(em2adj(dag_edgeM(st), p), B.R))
    dag_set(st, adj2em(A))
}

## Counting, sampling AND enumeration go through Clique-Picking and have no
## component-size limit: a complete DAG on 10 vertices is one chain component of 10 vertices
## whose class has 10! members, and both must handle it. Only the ENUMERATION
## used by the exhaustive escape is still capped, and it must decline rather
## than try.
set.seed(21)
p <- 10; A <- matrix(FALSE, p, p)
for (i in 1:(p-1)) for (j in (i+1):p) A[i, j] <- TRUE       # complete DAG
st <- dag_new(p); invisible(dag_set(st, adj2em(A)))
stopifnot(isTRUE(all.equal(imec_size(st, list(integer(0)), 4L), factorial(10))),
          is.null(imec_memb(st, list(integer(0)), 4L, 100)), # over budget: declines
          identical(imec_samp(st, list(integer(0)), 4L), TRUE))
## and what it sampled is a member: same skeleton, no immoralities
B <- em2adj(dag_edgeM(st), p)
stopifnot(all((B | t(B)) == (A | t(A))))
for (v in 1:p) { pa <- which(B[, v])
    if (length(pa) > 1) for (a in seq_along(pa)) for (b in seq_len(a - 1))
        stopifnot(B[pa[a], pa[b]] || B[pa[b], pa[a]]) }

cat("test_c_imec.R: OK\n")

## ---------------------------------------------------------------------------
## end to end: the two engines must agree on the exact-sampling path too, both
## in the DAG they return and in where they leave the RNG stream.  This is the
## same contract test_search_engines.R pins for the rcar() path.
suppressPackageStartupMessages(library(pcalg))
set.seed(4242)
for (it in 1:6) {
    p <- 8
    d <- r.gauss.pardag(p, 0.3, top.sort = TRUE, normalize = TRUE)
    tg <- c(list(integer(0)), lapply(sample(p, 2), function(v) as.integer(v)))
    nb <- c(120L, 40L, 40L)
    dat <- do.call(rbind, lapply(seq_along(tg), function(i)
        rmvnorm.ivent(nb[i], d, target = tg[[i]],
                      target.value = rep(2, length(tg[[i]])))))
    colnames(dat) <- as.character(seq_len(p))
    ti <- rep(seq_along(nb), nb)
    for (cf in list(c("exact", "trials"), c("exact", "exhaustive"),
                    c("rcar", "exhaustive"))) {
        set.seed(99 + it)
        hC <- hcmc(dat, 20, tg, ti, verbose = FALSE, engine = "C",
                   sampler = cf[1], escape = cf[2])
        seedC <- .Random.seed
        set.seed(99 + it)
        hR <- hcmc(dat, 20, tg, ti, verbose = FALSE, engine = "R",
                   sampler = cf[1], escape = cf[2])
        ## BITWISE, not all.equal: the two engines run the same scoring code
        ## on the same parent sets, so the only way the last bits can differ
        ## is if they present those sets in different orders. They did --
        ## sample_and_apply() re-orients only the changed arcs, leaving a
        ## history-dependent pa[]/ch[] order, while the R path rebuilds both
        ## ascending from the sampled matrix -- and the scores drifted apart
        ## in 35 of 40 runs until idl_dag_canonical_order() was added.
        stopifnot(identical(seedC, .Random.seed),
                  identical(hC$sco, hR$sco),
                  identical(dagadj(hC$dag), dagadj(hR$dag)))
    }
}
cat("test_c_imec.R: engines agree on the exact path\n")

## 5: C_dag_set_edges rejects a malformed edge list WITHOUT mutating the state.
## The validation used to happen on the way in, so a cyclic list raised only
## after the old arcs had been removed and the acyclic prefix inserted, and a
## repeated arc was not detected at all -- link_edge() has no duplicate guard,
## so it was pushed twice into pa[v] and ch[u] and counted twice in nedges,
## silently desynchronising the vectors from the adjacency bitsets.
mkem <- function(...) matrix(as.integer(c(...)), nrow = 2,
                             dimnames = list(c("from", "to"), NULL))
st   <- dag_new(5L)
good <- mkem(1,2, 2,3, 1,3, 4,5)
invisible(dag_set(st, good))
before <- dag_edgeM(st)
malformed <- list(
    cyclic     = mkem(1,2, 2,3, 3,1),      # cycle among three
    two.cycle  = mkem(1,2, 2,1),           # cycle of length two
    duplicate  = mkem(1,2, 1,2, 2,3),      # arc given twice
    dup.only   = mkem(3,4, 3,4),
    out.range  = mkem(1,2, 6,1),           # endpoint past p
    zero.index = mkem(0,1),                # 0-based index
    self.loop  = mkem(2,2),
    three.rows = matrix(1:6, nrow = 3),    # wrong shape, read as pairs before
    not.matrix = 1:4,                      # no dim at all
    doubles    = matrix(c(1,2,2,3), nrow = 2))
for (nm in names(malformed)) {
    err <- tryCatch({ dag_set(st, malformed[[nm]]); NULL },
                    error = function(e) conditionMessage(e))
    stopifnot(!is.null(err),                        # rejected
              identical(dag_edgeM(st), before))     # and nothing changed
}
invisible(dag_set(st, good))                        # still usable afterwards
stopifnot(identical(dag_edgeM(st), before))
cat("test_c_imec.R: malformed edge lists are rejected before any mutation\n")

## 6: the exact sampler has NO whole-graph size limit -- the 64-bit vertex
## sets are per chain component. R/imec.R used to guard with `p > 64L`, which
## made the R engine fall back to rcar() for any DAG on more than 64 vertices,
## however small its components: at p = 80 that was 288 fallbacks, a different
## DAG, a different score and a different RNG state from the C engine, i.e. a
## silent breach of the engines' documented equivalence.
A <- matrix(FALSE, 65L, 65L)                     # no undirected component
stopifnot(isTRUE(all.equal(idlBNs:::imec.size(A), 1)),
          !is.null(idlBNs:::imec.sample(A)))
A <- matrix(FALSE, 200L, 200L)                   # 100 components of two
for (i in seq(1L, 199L, by = 2L)) A[i, i + 1L] <- TRUE
stopifnot(isTRUE(all.equal(idlBNs:::imec.size(A), 2^100)),
          !is.null(idlBNs:::imec.sample(A)))

set.seed(21); p <- 70L; n <- 2000L
o <- sample(p); B <- matrix(0, p, p)
for (i in 1:(p-1)) for (j in (i+1):p)
    if (runif(1) < 4/p) B[o[i], o[j]] <- runif(1, .5, 1.5)
X <- matrix(0, n, p)
for (v in o) X[, v] <- X %*% B[, v] + rnorm(n)
colnames(X) <- as.character(seq_len(p))
set.seed(5); hC <- hcmc(X, verbose = FALSE, engine = "C", sampler = "exact")
seedC <- .Random.seed
set.seed(5); hR <- hcmc(X, verbose = FALSE, engine = "R", sampler = "exact")
stopifnot(identical(seedC, .Random.seed),                  # same random stream
          identical(hC$sco, hR$sco),                       # bitwise, see above
          identical(dagadj(hC$dag), dagadj(hR$dag)),
          hC$sampler.fallbacks == 0L, hR$sampler.fallbacks == 0L)
cat("test_c_imec.R: exact sampling is unrestricted in p, engines still agree\n")

## 7: a declined exact draw must consume no randomness, and the Get/PutRNGstate
## bracket must cover only the draws. cp_build() allocates with R_alloc(),
## which longjmps on failure; inside the bracket that would leave R's
## generator advanced while .Random.seed stayed behind, so a caller trapping
## the error would resume on a stale stream. Building every context first also
## makes a decline free: it used to leave the stream advanced by the
## components already sampled before a later one failed to build, even though
## the caller was told to fall back and would sample again.
local({
    q <- 70L                                   # one chain component of 70 > 64
    em <- matrix(as.integer(rbind(1:(q-1), 2:q)), nrow = 2,
                 dimnames = list(c("from", "to"), NULL))
    U <- matrix(FALSE, q, q); for (i in 1:(q-1)) U[i, i+1] <- U[i+1, i] <- TRUE
    st <- dag_new(q); invisible(dag_set(st, em))
    set.seed(1); before <- .Random.seed
    stopifnot(identical(.Call(idlBNs:::C_dag_imec_sample, st, list(integer(0))), FALSE),
              identical(.Random.seed, before))          # C engine: nothing drawn
    set.seed(1); before <- .Random.seed
    stopifnot(is.null(.Call(idlBNs:::C_cp_amo_sample, U)),
              identical(.Random.seed, before))          # R engine: nothing drawn
    A <- matrix(FALSE, q, q); for (i in 1:(q-1)) A[i, i+1] <- TRUE
    stopifnot(is.null(idlBNs:::imec.sample(A)))
    V <- matrix(FALSE, 6L, 6L); for (i in 1:5) V[i, i+1] <- V[i+1, i] <- TRUE
    set.seed(1); before <- .Random.seed
    invisible(.Call(idlBNs:::C_cp_amo_sample, V))
    stopifnot(!identical(.Random.seed, before))         # a real draw still advances
})
cat("test_c_imec.R: a declined draw consumes no randomness\n")

## cp_sample() reports its one failure ("no allowed permutation after 100000
## trials") as a 0 rather than raising, so the caller closes the RNG bracket
## before reporting it and nothing between GetRNGstate() and PutRNGstate()
## can longjmp. That path cannot be reached from R without injecting a fault,
## so it is not exercised here; it was verified out of tree by patching
## cp_sample() to return 0 on its second call and checking that the error
## surfaces AND .Random.seed is committed. To reproduce: add
##     { static int n = 0; if (getenv("IDLBNS_FAIL_DRAW") && ++n == 2)
##           return 0; }
## at the top of cp_sample(), rebuild, and draw from a graph with two chain
## components with IDLBNS_FAIL_DRAW set.

## 8: the exhaustive escape asks the enumeration ONCE per escape episode.
## A decline is followed by up to MAXTRIALS within-class trials, and those
## re-randomise the DAG inside its I-equivalence class -- whose I-essential
## graph, hence class size and chain components, is an invariant of the class.
## So the enumeration cannot answer differently, and re-attempting it after
## every trial recomputed a settled decision and counted the same refusal
## MAXTRIALS + 1 times (6 at the default, 21 at MAXTRIALS = 20).
local({
    set.seed(5); p <- 20L; n <- 1500L
    o <- sample(p); B <- matrix(0, p, p)
    for (i in 1:(p-1)) for (j in (i+1):p)
        if (runif(1) < 0.15) B[o[i], o[j]] <- runif(1, .5, 1.5)
    X <- matrix(0, n, p)
    for (v in o) X[, v] <- X %*% B[, v] + rnorm(n)
    colnames(X) <- as.character(seq_len(p))
    ## escape.max = 1 makes the enumeration decline at every local maximum,
    ## so the count is exactly the number of escape episodes. From this seed
    ## the search has one, whatever the trial budget: the count is 1 for both
    ## budgets, where before the fix it was MAXTRIALS + 1, i.e. 6 and 21.
    ## (The count is NOT generally independent of MAXTRIALS -- a different
    ## budget is a different trajectory, so a different number of episodes --
    ## which is why the seed is fixed here.)
    for (eng in c("C", "R")) {
        set.seed(9)
        f5  <- hcmc(X, verbose = FALSE, engine = eng, escape = "exhaustive",
                    escape.max = 1, MAXTRIALS = 5)
        set.seed(9)
        f20 <- hcmc(X, verbose = FALSE, engine = eng, escape = "exhaustive",
                    escape.max = 1, MAXTRIALS = 20)
        stopifnot(identical(f5$escape.fallbacks, 1L),
                  identical(f20$escape.fallbacks, 1L))
    }
    ## and an escape that always declines must be indistinguishable from
    ## escape = "trials": the enumeration consumes no randomness and, once
    ## refused, has no side effect on the state
    for (eng in c("C", "R")) {
        set.seed(3)
        a <- hcmc(X, verbose = FALSE, engine = eng, escape = "exhaustive",
                  escape.max = 1)
        sa <- .Random.seed
        set.seed(3)
        b <- hcmc(X, verbose = FALSE, engine = eng, escape = "trials")
        stopifnot(identical(a$sco, b$sco), identical(sa, .Random.seed),
                  identical(dagadj(a$dag), dagadj(b$dag)))
    }
})
## and the same for the exact sampler: its decline conditions are properties
## of the I-essential graph, hence of the class, so once it has declined it
## declines for every within-class trial. The DECISION is cached, but the
## count is not: every iteration still makes a draw, and it is the walk that
## makes it, so sampler.fallbacks must be one per draw exactly as before.
local({
    set.seed(2); p <- 80L; n <- 4000L
    X <- matrix(0, n, p); X[, 1] <- rnorm(n)
    for (v in 2:p) X[, v] <- 0.9 * X[, v-1] + rnorm(n) * 0.5
    colnames(X) <- as.character(seq_len(p))
    ## a chain has no immoralities, so its essential graph is one undirected
    ## component of p > 64 vertices and the exact sampler must decline
    set.seed(4); hC <- hcmc(X, verbose = FALSE, engine = "C", sampler = "exact")
    set.seed(4); hR <- hcmc(X, verbose = FALSE, engine = "R", sampler = "exact")
    stopifnot(hC$sampler.fallbacks > 0L,
              identical(hC$sampler.fallbacks, hR$sampler.fallbacks),
              identical(hC$sco, hR$sco),
              identical(dagadj(hC$dag), dagadj(hR$dag)))
})
cat("test_c_imec.R: the exhaustive escape is attempted once per episode\n")
