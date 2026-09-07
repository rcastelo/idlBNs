## Create the golden trajectory reference. Run ONCE, on a known-good tree,
## BEFORE any of the porting work starts -- it cannot be reconstructed later.
source(file.path("regression", "golden.R"))

if (file.exists(GOLDEN_FILE))
    stop("refusing to overwrite an existing ", GOLDEN_FILE,
         " -- delete it deliberately if you really mean to re-baseline")

cat("Generating golden reference (", length(GOLDEN_P) * length(GOLDEN_SCORE) *
    length(GOLDEN_ALG) * length(GOLDEN_SEED), " runs)\n", sep="")
ref <- golden_run()
write.csv(ref, GOLDEN_FILE, row.names=FALSE)
cat("\nWrote ", GOLDEN_FILE, " (", nrow(ref), " runs)\n", sep="")
cat("idlBNs version: ", as.character(packageVersion("idlBNs")), "\n", sep="")
cat("BLAS: ", sessionInfo()$BLAS, "\n", sep="")
