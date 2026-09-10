## 2026-09-10 differential test for the C port of the exact I-MEC machinery in
## src/imec.c against the R reference in R/imec.R.
##
## Three levels, in increasing strength:
##
##   1. |[D]_I| and the enumeration agree, as SETS.
##   2. The enumeration is duplicate free and exactly of size |[D]_I|, so it is
##      a bijection onto the class. Order is NOT pinned: the R reference walks
##      permutations of each chain component while the C engine unranks the
##      Clique-Picking decomposition, and the sampler draws its own clique
##      rather than indexing into the enumeration.
##   3. The sampler agrees bit for bit AND leaves .Random.seed in the same
##      place. As with rcar(), the stream is the load-bearing part: one
##      R_unif_index() per chain component, in component order, and R's
##      sample.int(n, 1) is (int) R_unif_index(n) + 1.

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
rdag   <- function(p, prob) { o <- sample(p); A <- matrix(FALSE, p, p)
    for (i in 1:(p-1)) for (j in (i+1):p) if (runif(1) < prob) A[o[i], o[j]] <- TRUE; A }
rtgts  <- function(p, k) if (k == 0) list(integer(0)) else
    c(list(integer(0)), lapply(sample(p, k), function(v) as.integer(v)))

set.seed(7)
for (it in 1:150) {
    p <- sample(4:8, 1); A <- rdag(p, 0.4); tg <- rtgts(p, sample(0:2, 1))
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
        stopifnot(identical(seedC, .Random.seed),
                  isTRUE(all.equal(hC$sco, hR$sco)),
                  identical(sort(as.vector(graph::edgeMatrix(hC$dag))),
                            sort(as.vector(graph::edgeMatrix(hR$dag)))))
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
dag_set(st, good)
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
dag_set(st, good)                                   # still usable afterwards
stopifnot(identical(dag_edgeM(st), before))
cat("test_c_imec.R: malformed edge lists are rejected before any mutation\n")
