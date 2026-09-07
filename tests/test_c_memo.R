## 2026-09-07 tests for the addition memo.
##
## WHY IT EXISTS. Even with the compiled hash cache, the lookup is the
## dominant per-step cost at large p: measured 26.6 ns per candidate at
## p = 400 against a 2.1 ns floor for a bare array read. And additions are
## about 99% of the candidates -- at p = 400 with 800 arcs, 154,000 of
## 155,764 -- because removals and reversals number one per arc. So the memo
## covers additions only, keyed on (u, v).
##
## WHY IT IS SAFE, which is the part that needs testing. An addition's delta
## is a function of pa(v) alone, through both s(v, pa(v) + u) and base[v], so
## a single compare against the DAG's pav_stamp[v] decides validity. And a
## memo hit SKIPS A CACHE LOOKUP THAT WOULD HAVE BEEN A HIT -- the key is
## the one the earlier computation used, the cache never evicts, so the
## lookup could not have missed. Skipping a hit stores nothing and changes
## nothing, which is why the memo does not perturb the cache's contents and
## therefore does not perturb which order-variant of a parent set is stored.
## That last point is the whole hazard, and section 2 tests it directly.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

sc_new    <- function(p) .Call(idlBNs:::C_sccache_new, as.integer(p))
sc_dump   <- function(x) .Call(idlBNs:::C_sccache_dump, x)
sc_stats  <- function(x) .Call(idlBNs:::C_sccache_stats, x)
dag_new   <- function(p) .Call(idlBNs:::C_dag_new, as.integer(p))
dag_move  <- function(st, op, u, v)
    invisible(.Call(idlBNs:::C_dag_apply_move, st, as.integer(op),
                    as.integer(u), as.integer(v)))
dag_pas   <- function(st) .Call(idlBNs:::C_dag_pasets, st)
dag_stamp <- function(st) .Call(idlBNs:::C_dag_pastamp, st)
dag_nh    <- function(st, kind, ut = integer(0))
    .Call(idlBNs:::C_dag_nh, st, as.integer(kind), as.integer(ut))

################################################################################
## 1. the stamps themselves: pav_stamp[v] must change exactly when pa(v)
## does, and nothing else must change
################################################################################

p <- 6
st <- dag_new(p)
s0 <- dag_stamp(st)
stopifnot(identical(s0, rep(1L, p)))        ## starts at 1, never 0

dag_move(st, 1L, 1L, 3L)                    ## add 1 -> 3: only pa(3) changes
s1 <- dag_stamp(st)
stopifnot(s1[3] != s0[3], identical(s1[-3], s0[-3]))

dag_move(st, 1L, 2L, 3L)                    ## add 2 -> 3: pa(3) again
s2 <- dag_stamp(st)
stopifnot(s2[3] != s1[3], identical(s2[-3], s1[-3]))

dag_move(st, 2L, 1L, 3L)                    ## remove 1 -> 3
s3 <- dag_stamp(st)
stopifnot(s3[3] != s2[3], identical(s3[-3], s2[-3]))

dag_move(st, 3L, 2L, 3L)                    ## reverse 2 -> 3: BOTH change
s4 <- dag_stamp(st)
stopifnot(s4[2] != s3[2], s4[3] != s3[3],
          identical(s4[-c(2, 3)], s3[-c(2, 3)]))
cat("pav_stamp changes exactly when a vertex's parent set does\n")

################################################################################
## 2. THE CENTRAL PROPERTY: the memo must not change the score cache's
## contents. Two identical searches, one with the memo enabled and one
## without, must end with byte-identical cache contents -- because a memo hit
## only ever skips a lookup that would have been a hit.
################################################################################

mkdata <- function(p, k = 2L, n = 150L, seed = 1L) {
  set.seed(seed)
  Mg <- r.gauss.pardag(p, 0.3, top.sort = TRUE, normalize = TRUE)
  I <- c(list(integer(0)), sample(p, size = k, replace = FALSE))
  nb <- rep(floor(n / (k + 1)), k)
  nb <- c(n - sum(nb), nb)
  dat <- do.call(rbind,
                 lapply(seq_along(I),
                        function(v) rmvnorm.ivent(nb[v], Mg, target = I[[v]],
                                                  target.value = rep(2, length(I[[v]])))))

  list(dat = dat, targets = I, target.index = rep(seq_along(nb), nb))
}

norm_hash <- function(x)
    lapply(sc_dump(x), function(v) v[order(names(v))])

## drive a search by hand, so the memo can be turned on and off while
## everything else -- the DAG, the move sequence, the cache -- stays identical
manual_search <- function(x, sf, p, use.memo, nsteps = 40L) {
  gs <- attr(sf, "global.sufstats.fun")(x$dat, x$targets, x$target.index)
  cs <- sc_new(p)
  st <- dag_new(p)
  amf <- attr(sf, "nh.argmax.fun")
  ut <- as.integer(sort(unique(unlist(x$targets))))
  for (k in seq_len(nsteps)) {
    ne <- dag_nh(st, 3L, ut)
    if (length(ne$op) == 0L)
      break
    pas <- dag_pas(st)
    am <- amf(ne$op, ne$u, ne$v, pas, gs, cs, FALSE,
              if (use.memo) dag_stamp(st) else NULL)
    dag_move(st, ne$op[am$index], ne$u[am$index], ne$v[am$index])
  }

  list(cache = norm_hash(cs), stats = sc_stats(cs))
}

ncase <- 0L
for (p in c(6, 10, 14))
  for (seed in 1:2) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe)) {
      a <- manual_search(x, sf, p, use.memo = TRUE)
      b <- manual_search(x, sf, p, use.memo = FALSE)
      ## the cache contents must be byte-identical -- same keys, same values
      stopifnot(identical(a$cache, b$cache))
      ## the memo really was used
      stopifnot(a$stats[["memo.hits"]] > 0)
      stopifnot(b$stats[["memo.hits"]] == 0)
      ## and it saved real lookups
      stopifnot(a$stats[["hits"]] + a$stats[["misses"]] <
                b$stats[["hits"]] + b$stats[["misses"]])
      ## while the ENTRIES stored are the same count
      stopifnot(a$stats[["entries"]] == b$stats[["entries"]])
      ncase <- ncase + 1L
    }
  }
stopifnot(ncase >= 12L)
cat(sprintf("memo does not perturb the cache: %d searches, contents byte-identical\n",
            ncase))

################################################################################
## 3. full searches: both engines still agree, bit for bit, with the RNG
################################################################################

nsearch <- 0L
for (p in c(8, 14, 20))
  for (seed in 1:2) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe))
      for (alg in c("hcmc", "hillclimbing")) {
        f <- match.fun(alg)
        set.seed(77L)
        rC <- f(x$dat, targets = x$targets, target.index = x$target.index,
                scorefun = sf, verbose = FALSE, engine = "C")
        sC <- .Random.seed
        set.seed(77L)
        rR <- f(x$dat, targets = x$targets, target.index = x$target.index,
                scorefun = sf, verbose = FALSE, engine = "R")
        stopifnot(identical(rC$sco, rR$sco),
                  identical(unname(edgeMatrix(rC$dag)),
                            unname(edgeMatrix(rR$dag))),
                  identical(sC, .Random.seed))
        nsearch <- nsearch + 1L
      }
  }
cat(sprintf("searches: %d cases bit-identical between engines with the memo on\n",
            nsearch))

################################################################################
## 4. with the band verified at every step, which also re-checks the bound
## now that it is computed from the delta alone rather than from s
################################################################################

old <- getOption("idlBNs.debug.band", FALSE)
on.exit(options(idlBNs.debug.band = old), add = TRUE)
options(idlBNs.debug.band = TRUE)
x <- mkdata(12, seed = 5L)
for (sf in list(iBIC, iBGe))
  for (alg in c("hcmc", "hillclimbing")) {
    f <- match.fun(alg)
    set.seed(5L)
    rC <- f(x$dat, targets = x$targets, target.index = x$target.index,
            scorefun = sf, verbose = FALSE, engine = "C")
    set.seed(5L)
    rR <- f(x$dat, targets = x$targets, target.index = x$target.index,
            scorefun = sf, verbose = FALSE, engine = "R")
    stopifnot(identical(rC$sco, rR$sco),
              identical(unname(edgeMatrix(rC$dag)),
                        unname(edgeMatrix(rR$dag))))
  }
options(idlBNs.debug.band = old)
cat("band verified at every step with the memo enabled\n")

################################################################################
## 5. the memo stays off unless BOTH a compiled cache and the stamps are
## supplied, since without stamps a stale entry is indistinguishable from a
## fresh one
################################################################################

p <- 8
x <- mkdata(p, seed = 3L)
gs <- idlBNs:::.iBIC.global.sufstats(x$dat, x$targets, x$target.index)
st <- dag_new(p)
dag_move(st, 1L, 1L, 2L)
dag_move(st, 1L, 2L, 3L)
ne <- dag_nh(st, 2L)
pas <- dag_pas(st)

## no stamps -> memo off
cs <- sc_new(p)
invisible(idlBNs:::.iBIC.nh.argmax(ne$op, ne$u, ne$v, pas, gs, cs, FALSE, NULL))
stopifnot(sc_stats(cs)[["memo.hits"]] == 0,
          sc_stats(cs)[["memo.misses"]] == 0)
## environment cache + stamps -> memo off too (nowhere to keep it)
csE <- lapply(seq_len(p), function(i) new.env(hash = TRUE, parent = emptyenv()))
a <- idlBNs:::.iBIC.nh.argmax(ne$op, ne$u, ne$v, pas, gs, csE, FALSE,
                              dag_stamp(st))
## and the answer is the same either way
b <- idlBNs:::.iBIC.nh.argmax(ne$op, ne$u, ne$v, pas, gs, cs, FALSE,
                              dag_stamp(st))
d <- idlBNs:::.iBIC.nh.argmax(ne$op, ne$u, ne$v, pas, gs, NULL, FALSE, NULL)
stopifnot(identical(a$index, b$index), identical(a$total, b$total),
          identical(a$index, d$index), identical(a$total, d$total))
## a wrong-length stamp vector is rejected
stopifnot(inherits(tryCatch(idlBNs:::.iBIC.nh.argmax(ne$op, ne$u, ne$v, pas,
                                                     gs, cs, FALSE,
                                                     rep(1L, p - 1L)),
                            error = function(e) e), "error"))
cat("memo requires both a compiled cache and stamps; answer unchanged either way\n")

cat("all addition memo tests passed\n")
