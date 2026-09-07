## Compare the current tree against the golden trajectory reference.
## Run at every implementation gate. Exits non-zero on any divergence.
source(file.path("regression", "golden.R"))

if (!file.exists(GOLDEN_FILE))
    stop(GOLDEN_FILE, " not found -- run regression/make_reference.R first")

ref <- read.csv(GOLDEN_FILE, stringsAsFactors=FALSE, colClasses="character")
cat("Re-running ", nrow(ref), " golden runs\n", sep="")
now <- golden_run()

key <- function(d) paste(d$p, d$score, d$alg, d$seed, sep="/")
ref <- ref[match(key(now), key(ref)), ]
stopifnot(!any(is.na(ref$sco)))

bad <- 0L
for (i in seq_len(nrow(now))) {
    dsco <- !identical(now$sco[i], ref$sco[i])
    dedg <- !identical(now$edges[i], ref$edges[i])
    if (dsco || dedg) {
        bad <- bad + 1L
        cat(sprintf("\nDIVERGED %s\n", key(now)[i]))
        if (dsco) cat(sprintf("  score : golden %s\n          now    %s\n",
                              ref$sco[i], now$sco[i]))
        if (dedg) cat(sprintf("  edges : golden %s arcs\n          now    %s arcs\n",
                              ref$nedges[i], now$nedges[i]))
    }
}

cat(sprintf("\n%d/%d runs bit-identical to the golden reference\n",
            nrow(now) - bad, nrow(now)))
if (bad > 0L)
    quit(status=1L)
