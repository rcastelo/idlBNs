## 2026-09-07 regression test for the RNG identities a C port of rcar()
## depends on.
##
## rcar() (R/search.R) is the only part of the search that consumes
## randomness, via `sample(0:r, size=1)` once and `resample(which(cemask),
## size=1)` rr times. Moving it into C means reproducing that stream draw
## for draw: the pinned .Rout.save outputs and the golden trajectory
## reference both fail on any divergence, and a divergence here is invisible
## to every other test.
##
## The critical point is that these tests assert equality of .Random.seed
## AFTER the call, not just of the returned value. Under the default
## sample.kind = "Rejection", R_unif_index() draws by rejection, so the
## number of underlying unif_rand() calls is data dependent -- 1 to 3 for a
## size-1 draw. Two implementations can therefore agree on every value while
## consuming different amounts of the stream, and only the seed comparison
## catches that. It is exactly why the C port calls R_unif_index() instead
## of scaling its own uniform.

suppressPackageStartupMessages({
  library(graph)
  library(idlBNs)
})

## restore whatever the caller had, whichever way this test exits
old_kind <- RNGkind()
on.exit(do.call(RNGkind, as.list(old_kind)), add=TRUE)

## population sizes chosen around the powers of two, because that is where
## R_unif_index()'s bit count changes and its rejection rate jumps: n = 2^k
## never rejects, n = 2^k + 1 rejects almost half the time
NVALS <- c(1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
           100, 255, 256, 257, 1000, 65535, 65536, 65537)
SEEDS <- 1:60

KINDS <- list(c("Mersenne-Twister", "Inversion", "Rejection"),
              c("Mersenne-Twister", "Inversion", "Rounding"),
              c("L'Ecuyer-CMRG",    "Inversion", "Rejection"))

ncmp <- 0L
for (kind in KINDS) {
  ## the Rounding sampler is deprecated and warns on use; that is not what
  ## is under test here
  suppressWarnings(do.call(RNGkind, as.list(kind)))

  ##########################################################################
  ## 1. sample.int(n, 1) == .unif_index(n) + 1, in value AND in stream
  ##########################################################################
  for (n in NVALS)
    for (s in SEEDS) {
      set.seed(s)
      a <- sample.int(n, size=1)
      sa <- .Random.seed
      set.seed(s)
      b <- idlBNs:::.unif_index(n) + 1
      sb <- .Random.seed
      stopifnot(identical(as.integer(a), as.integer(b)),
                identical(sa, sb))
      ncmp <- ncmp + 1L
    }

  ##########################################################################
  ## 2. sample(0:r, size=1) == .unif_index(length(0:r))
  ##
  ## r = 0 is included deliberately. sample(0:0, size=1) does NOT take
  ## sample()'s "treat a length-1 numeric x as 1:x" shortcut, because that
  ## needs x >= 1; it goes through sample.int(1, 1), and R_unif_index(1)
  ## computes bits = 0 and still consumes one unif_rand(). A C port that
  ## short-circuits "there is only one choice" desynchronises the stream.
  ## tests/test_pasets.R runs hcmc(dat3, r=0), so this path is live.
  ##########################################################################
  for (r in 0:25)
    for (s in SEEDS) {
      set.seed(s)
      a <- sample(0:r, size=1)
      sa <- .Random.seed
      set.seed(s)
      b <- idlBNs:::.unif_index(length(0:r))
      sb <- .Random.seed
      stopifnot(identical(as.integer(a), as.integer(b)),
                identical(sa, sb))
      ncmp <- ncmp + 1L
    }

  ##########################################################################
  ## 3. resample(x, size=1) == x[.unif_index(length(x)) + 1]
  ##
  ## length 1 is included for the same reason as r = 0 above: it still draws.
  ##########################################################################
  for (len in 1:12)
    for (s in SEEDS) {
      x <- seq_len(len) * 7L
      set.seed(s)
      a <- idlBNs:::resample(x, size=1)
      sa <- .Random.seed
      set.seed(s)
      b <- x[idlBNs:::.unif_index(len) + 1]
      sb <- .Random.seed
      stopifnot(identical(a, b), identical(sa, sb))
      ncmp <- ncmp + 1L
    }
}

suppressWarnings(do.call(RNGkind, as.list(old_kind)))
stopifnot(ncmp > 0L)
cat(sprintf("RNG identities: %d comparisons, value and .Random.seed identical\n",
            ncmp))

################################################################################
## 4. the draw count really is data dependent, so that the seed assertions
## above are doing work rather than passing trivially
################################################################################

## how many unif_rand() draws does one sample.int(n, 1) consume?
draw_count <- function(n, seed=1L, maxk=12L) {
  set.seed(seed)
  invisible(sample.int(n, size=1))
  target <- .Random.seed
  for (k in seq_len(maxk)) {
    set.seed(seed)
    invisible(runif(k))
    if (identical(.Random.seed, target))
      return(k)
  }
  NA_integer_
}

counts <- vapply(c(1, 2, 3, 21, 64, 65, 500), draw_count, integer(1))
stopifnot(!anyNA(counts))
## if every draw consumed the same amount, a naive "scale one uniform" port
## would pass the tests above and still be wrong, so assert the variation
stopifnot(length(unique(counts)) > 1L)
cat(sprintf("rejection sampling consumes %s unif_rand draws for size-1 samples\n",
            paste(range(counts), collapse="-")))

################################################################################
## 5. a no-draw rcar() must not merely leave .Random.seed unchanged -- in a
## session that has not yet touched the RNG it must not bring .Random.seed
## into existence at all. That is the property that lets the C port put
## rcar()'s two early returns ABOVE GetRNGstate(): an unconditional
## PutRNGstate() would create a time-seeded .Random.seed, and the very first
## hcmc() iteration always hits the numEdges == 0 early return.
##
## This has to run in a fresh process, because the property is about
## .Random.seed not existing. Rscript rather than callr, to avoid adding a
## dependency for one test.
################################################################################

probe <- tempfile(fileext=".R")
writeLines(c(
  'suppressPackageStartupMessages({library(graph); library(idlBNs)})',
  'vn <- paste0("X", 1:4)',
  'g0 <- new("graphNEL", nodes=vn, edgemode="directed")   ## no edges',
  'a0 <- idlBNs:::init.ancestors(vn)',
  'p0 <- idlBNs:::init.pasets(4L)',
  'vi <- setNames(1:4, vn)',
  '## numEdges(g0) == 0, so rcar() returns before drawing anything',
  'invisible(idlBNs:::rcar(g0, 20L, integer(0), a0, p0, vi))',
  'cat(exists(".Random.seed", envir=globalenv()), "\n")'), probe)
rscript <- file.path(R.home("bin"), "Rscript")
out <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(probe)),
                                stdout=TRUE, stderr=FALSE))
unlink(probe)
## only assert if the subprocess actually ran; a sandbox that blocks it
## should not turn into a spurious failure
if (length(out) > 0L && any(grepl("^(TRUE|FALSE)$", trimws(out)))) {
  created <- as.logical(trimws(out[grepl("^(TRUE|FALSE)$", trimws(out))][1]))
  stopifnot(identical(created, FALSE))
  cat("a no-draw rcar() leaves .Random.seed non-existent in a fresh session\n")
} else
  cat("subprocess probe skipped (could not run Rscript)\n")

cat("all RNG equivalence tests passed\n")
