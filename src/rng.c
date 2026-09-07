#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

/*
 * RNG SUPPORT FOR THE C PORT OF rcar()
 *
 * rcar() is the only part of the search that consumes randomness, through
 * two expressions in R/search.R:
 *
 *     rr    <- sample(0:r, size = 1)
 *     rndce <- resample(which(cemask), size = 1)   # x[sample.int(n, 1)]
 *
 * A C port has to reproduce that stream draw for draw, or every trajectory
 * moves and the pinned .Rout.save outputs under tests/ fail. The way to do
 * that is to CALL R's own R_unif_index() rather than to reimplement it:
 *
 *   - Under the default sample.kind = "Rejection", R_unif_index() draws
 *     ceil(log2(n)) bits by taking floor(bits/16)+1 unif_rand() values and
 *     retries while the result is >= n. The number of unif_rand() calls is
 *     therefore DATA DEPENDENT -- measured 1 to 3 for size-1 draws -- so
 *     "scale one uniform" desynchronises the stream immediately.
 *   - It also dispatches internally on R_sample_kind(), so calling it stays
 *     correct if the user has set sample.kind = "Rounding".
 *
 * The identities the port relies on, all pinned by
 * tests/test_rng_equivalence.R:
 *
 *     sample.int(n, 1)       ==  (int) R_unif_index((double) n) + 1
 *     sample(0:r, size = 1)  ==  (int) R_unif_index((double) length(0:r))
 *     resample(x, size = 1)  ==  x[(int) R_unif_index((double) length(x)) + 1]
 *
 * Note that length(0:r) is passed in from R rather than computed as r + 1
 * in C: 0:r for a non-integer r is 0:floor(r), and for negative r it counts
 * down, so leaving that coercion on the R side keeps the corner cases where
 * they already behave correctly.
 */

/*
 * C_unif_index
 *
 * Test-only wrapper exposing one R_unif_index() draw to R, so that the
 * three identities above can be asserted on both the returned value and
 * the resulting .Random.seed. Not used by the search itself.
 *
 * n_R  REALSXP  scalar: the exclusive upper bound, i.e. the population size
 *
 * Returns a length-1 REALSXP in [0, n-1], or 0 when n <= 0 (matching
 * R_unif_index(), which returns 0 without drawing in that case).
 */
SEXP
C_unif_index(SEXP n_R) {
    if (TYPEOF(n_R) != REALSXP || XLENGTH(n_R) != 1)
        error("C_unif_index: 'n' must be a numeric scalar");

    double n = REAL(n_R)[0];

    if (!R_FINITE(n))
        error("C_unif_index: 'n' must be finite");

    /* GetRNGstate()/PutRNGstate() bracket the draw, which is what keeps
       .Random.seed in step with the C-side generator. In the real rcar()
       port this bracket goes around the WHOLE call, with the two early
       returns above GetRNGstate() so that a no-draw call leaves the R
       generator untouched -- and, in a session that has not used the RNG
       yet, does not bring .Random.seed into existence at all. */
    GetRNGstate();
    double v = R_unif_index(n);
    PutRNGstate();

    return Rf_ScalarReal(v);
}
