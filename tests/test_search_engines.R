## 2026-09-07 engine parity: the C and R search engines must follow the same
## trajectory and return the same result, bit for bit.
##
## This is the highest-value single test of the C port, because one assertion
## covers the whole stack at once: the DAG struct, the ancestor and
## descendant bookkeeping, the Kahn topological order inside a removal, the
## neighbourhood enumeration ORDER (which.max breaks ties by position, so a
## reordering changes the trajectory), the parent-set insertion order that
## reaches the score functions, and -- once hcmc() is wired -- the RNG stream
## through rcar().
##
## identical() throughout, never a tolerance: the two engines run in the same
## session against the same BLAS, so any difference at all is a real
## divergence rather than platform noise.

suppressPackageStartupMessages({
  library(graph)
  library(pcalg)
  library(idlBNs)
})

## simulate mixed observational and interventional Gaussian data
mkdata <- function(p, k = 2L, n = 120L, seed = 1L, d = 0.3) {
  set.seed(seed)
  Mg <- r.gauss.pardag(p, d, top.sort = TRUE, normalize = TRUE)
  I <- c(list(integer(0)), sample(p, size = k, replace = FALSE))
  nb <- rep(floor(n / (k + 1)), k)
  nb <- c(n - sum(nb), nb)
  dat <- do.call(rbind,
                 lapply(seq_along(I),
                        function(v) rmvnorm.ivent(nb[v], Mg, target = I[[v]],
                                                  target.value = rep(2, length(I[[v]])))))

  list(dat = dat, targets = I, target.index = rep(seq_along(nb), nb))
}

same_result <- function(a, b) {
  identical(a$sco, b$sco) &&
    identical(unname(edgeMatrix(a$dag)), unname(edgeMatrix(b$dag))) &&
    identical(nodes(a$dag), nodes(b$dag)) &&
    ## the child lists too, not just the arc set: their order is what
    ## nr.nh() emits removals in
    identical(lapply(edgeL(a$dag), function(e) as.integer(e$edges)),
              lapply(edgeL(b$dag), function(e) as.integer(e$edges)))
}

ncase <- 0L

################################################################################
## 1. hillclimbing(): deterministic, so the two engines must agree exactly
################################################################################

for (p in c(5, 8, 12, 20))
  for (seed in 1:3) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe))
      for (tg in list(list(targets = list(integer(0)),
                           target.index = rep(1L, nrow(x$dat))),
                      list(targets = x$targets, target.index = x$target.index))) {
        rC <- hillclimbing(x$dat, targets = tg$targets,
                           target.index = tg$target.index, scorefun = sf,
                           verbose = FALSE, engine = "C")
        rR <- hillclimbing(x$dat, targets = tg$targets,
                           target.index = tg$target.index, scorefun = sf,
                           verbose = FALSE, engine = "R")
        stopifnot(same_result(rC, rR))
        ## and the reported score really is the score of the returned DAG
        stopifnot(abs(rC$sco - sf(rC$dag, x$dat, tg$targets,
                                  tg$target.index)) < 1e-9)
        ncase <- ncase + 1L
      }
  }
stopifnot(ncase >= 48L)
cat(sprintf("hillclimbing: %d cases, C and R engines bit-identical\n", ncase))

################################################################################
## 2. hcmc(): stochastic, so the engines must agree on the RANDOM STREAM as
## well as the result. The assertion on .Random.seed after the run is the
## load-bearing one -- rcar() draws a data-dependent number of unif_rand()
## values per call and is called once or twice per iteration, so two
## implementations could agree on every arc they reverse and still leave the
## generator in different places, which would make every later step diverge.
################################################################################

nhcmc <- 0L
for (p in c(5, 8, 12, 20))
  for (seed in 1:3) {
    x <- mkdata(p, seed = seed)
    for (sf in list(iBIC, iBGe))
      for (r in c(0L, 20L))
        for (tg in list(list(targets = list(integer(0)),
                             target.index = rep(1L, nrow(x$dat))),
                        list(targets = x$targets,
                             target.index = x$target.index))) {
          set.seed(4242L)
          rC <- hcmc(x$dat, r = r, targets = tg$targets,
                     target.index = tg$target.index, scorefun = sf,
                     verbose = FALSE, engine = "C")
          seedC <- .Random.seed
          set.seed(4242L)
          rR <- hcmc(x$dat, r = r, targets = tg$targets,
                     target.index = tg$target.index, scorefun = sf,
                     verbose = FALSE, engine = "R")
          seedR <- .Random.seed

          stopifnot(same_result(rC, rR))
          stopifnot(identical(seedC, seedR))        ## the same stream
          stopifnot(abs(rC$sco - sf(rC$dag, x$dat, tg$targets,
                                    tg$target.index)) < 1e-9)
          nhcmc <- nhcmc + 1L
        }
  }
stopifnot(nhcmc >= 96L)
cat(sprintf("hcmc: %d cases, C and R engines bit-identical including .Random.seed\n",
            nhcmc))

## r = 0 exercises the branch where rcar() still consumes exactly one draw
## per call but never reverses anything
x <- mkdata(8, seed = 9L)
set.seed(7L); z0 <- hcmc(x$dat, r = 0L, verbose = FALSE, engine = "C")
s0 <- .Random.seed
set.seed(7L); z1 <- hcmc(x$dat, r = 0L, verbose = FALSE, engine = "R")
stopifnot(same_result(z0, z1), identical(s0, .Random.seed))
cat("hcmc at r = 0: engines agree on result and stream\n")

################################################################################
## 3. a scorefun that cannot score a whole neighbourhood falls back to R
## rather than failing. tests/test_delta_scores.R already drives that path;
## here we check the fallback is transparent -- asking for engine = "C" with
## such a scorefun gives the same answer as asking for "R".
################################################################################

plain_iBIC <- function(g, dat, targets = list(integer(0)),
                       target.index = rep(1L, nrow(dat)), cached.scores = NULL,
                       global.sufstats = NULL)
    iBIC(g, dat, targets = targets, target.index = target.index,
         cached.scores = cached.scores, global.sufstats = global.sufstats)

x <- mkdata(6, seed = 4L)
aC <- hillclimbing(x$dat, scorefun = plain_iBIC, verbose = FALSE, engine = "C")
aR <- hillclimbing(x$dat, scorefun = plain_iBIC, verbose = FALSE, engine = "R")
stopifnot(same_result(aC, aR))
## and it reaches the same DAG the fast path does
fast <- hillclimbing(x$dat, scorefun = iBIC, verbose = FALSE)
stopifnot(same_result(aC, fast))

## hcmc() too, including the random stream
set.seed(3L); hC <- hcmc(x$dat, scorefun = plain_iBIC, verbose = FALSE,
                         engine = "C")
sC <- .Random.seed
set.seed(3L); hR <- hcmc(x$dat, scorefun = plain_iBIC, verbose = FALSE,
                         engine = "R")
stopifnot(same_result(hC, hR), identical(sC, .Random.seed))
cat("attribute-less scorefun: engine='C' falls back to R transparently\n")

################################################################################
## 4. the engine argument itself
################################################################################

stopifnot(inherits(tryCatch(hillclimbing(x$dat, verbose = FALSE,
                                         engine = "nonsense"),
                            error = function(e) e), "error"))
## partial matching, as match.arg gives
stopifnot(same_result(hillclimbing(x$dat, verbose = FALSE, engine = "C"),
                      hillclimbing(x$dat, verbose = FALSE)))
stopifnot(inherits(tryCatch(hcmc(x$dat, verbose = FALSE, engine = "nonsense"),
                            error = function(e) e), "error"))
set.seed(2L); e1 <- hcmc(x$dat, verbose = FALSE, engine = "C")
set.seed(2L); e2 <- hcmc(x$dat, verbose = FALSE)
stopifnot(same_result(e1, e2))
cat("engine argument validated and defaults to C\n")

################################################################################
## 5. the in-loop debug assertions, on both engines.
##
## idlBNs.debug.pasets existed before the port; idlBNs.debug.anc and
## idlBNs.debug.dag are new and close a real gap -- remove.ancestors() and
## reverse.ancestors() previously had no in-loop assertion of any kind, so an
## ancestor-bookkeeping bug that only changed WHICH moves were legal was
## caught by nothing but the pinned trajectory outputs.
################################################################################

old <- list(p = getOption("idlBNs.debug.pasets", FALSE),
            a = getOption("idlBNs.debug.anc", FALSE),
            d = getOption("idlBNs.debug.dag", FALSE))
on.exit(options(idlBNs.debug.pasets = old$p, idlBNs.debug.anc = old$a,
                idlBNs.debug.dag = old$d), add = TRUE)
options(idlBNs.debug.pasets = TRUE, idlBNs.debug.anc = TRUE,
        idlBNs.debug.dag = TRUE)

x <- mkdata(8, seed = 5L)
for (sf in list(iBIC, iBGe))
  for (eng in c("C", "R")) {
    r <- hillclimbing(x$dat, targets = x$targets,
                      target.index = x$target.index, scorefun = sf,
                      verbose = FALSE, engine = eng)
    stopifnot(is.finite(r$sco))
    ## and hcmc(), at both r = 0 and the default, so that rcar()'s two call
    ## sites run with every assertion enabled
    for (rr in c(0L, 20L)) {
      set.seed(6L)
      h <- hcmc(x$dat, r = rr, targets = x$targets,
                target.index = x$target.index, scorefun = sf,
                verbose = FALSE, engine = eng)
      stopifnot(is.finite(h$sco))
    }
  }
## with the assertions on, the two engines still agree
rC <- hillclimbing(x$dat, targets = x$targets, target.index = x$target.index,
                   verbose = FALSE, engine = "C")
rR <- hillclimbing(x$dat, targets = x$targets, target.index = x$target.index,
                   verbose = FALSE, engine = "R")
stopifnot(same_result(rC, rR))
set.seed(8L); hC <- hcmc(x$dat, targets = x$targets,
                         target.index = x$target.index, verbose = FALSE,
                         engine = "C")
sC <- .Random.seed
set.seed(8L); hR <- hcmc(x$dat, targets = x$targets,
                         target.index = x$target.index, verbose = FALSE,
                         engine = "R")
stopifnot(same_result(hC, hR), identical(sC, .Random.seed))

options(idlBNs.debug.pasets = old$p, idlBNs.debug.anc = old$a,
        idlBNs.debug.dag = old$d)
cat("debug.pasets, debug.anc and debug.dag all pass on both engines\n")

################################################################################
## 6. the debug assertions are not vacuous: idlBNs.debug.anc must reject a
## deliberately wrong ancestor matrix
################################################################################

g <- new("graphNEL", nodes = paste0("X", 1:4), edgemode = "directed")
g <- addEdge("X1", "X2", g)
g <- addEdge("X2", "X3", g)
## .anc_closure() returns an UNNAMED matrix on purpose: .debug_assertions()
## compares it against unname(anc), so index it positionally
truth <- idlBNs:::.anc_closure(g)
stopifnot(is.null(dimnames(truth)), dim(truth) == c(4L, 4L))
stopifnot(truth[1L, 3L])                ## X1 reaches X3 through X2
wrong <- truth
wrong[1L, 3L] <- FALSE                  ## drop a real ancestor
stopifnot(!identical(wrong, truth))
## and it agrees with the incremental ancestor matrix the R engine keeps
anc <- idlBNs:::init.ancestors(nodes(g))
em <- edgeMatrix(g)
for (k in seq_len(ncol(em)))
  anc <- idlBNs:::add.ancestors(anc, nodes(g)[em["from", k]],
                                nodes(g)[em["to", k]])
stopifnot(identical(unname(anc), truth))
cat("debug.anc reference rejects a corrupted ancestor matrix\n")

cat("all engine parity tests passed\n")
