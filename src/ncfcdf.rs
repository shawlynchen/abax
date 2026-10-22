use crate::{betainc, fcdf, gammaln};

/// Noncentral F cumulative distribution function.
///
/// `upper`: if true, return the upper-tail probability P(X > x).
pub fn ncfcdf(x: f64, nu1: f64, nu2: f64, delta: f64, upper: bool) -> f64 {
    let cvgtol = f64::EPSILON.powf(0.75);

    // NaN handling
    if x.is_nan() || nu1.is_nan() || nu2.is_nan() || delta.is_nan() {
        return f64::NAN;
    }

    // Illegal parameters
    let k = nu1 <= 0.0 || nu2 <= 0.0 || delta < 0.0;
    if k {
        return f64::NAN;
    }

    // x == Inf
    if x.is_infinite() && x > 0.0 {
        return if upper { 0.0 } else { 1.0 };
    }

    // Central distribution
    if delta == 0.0 {
        return fcdf(x, nu1, nu2, upper);
    }

    // Upper tail: p(x <= 0) = 1
    if upper && x <= 0.0 {
        return 1.0;
    }

    // Lower tail: p(x <= 0) = 0
    if !upper && x <= 0.0 {
        return 0.0;
    }

    // Pre-divide nu1, nu2, delta
    let nu1 = nu1 / 2.0;
    let nu2 = nu2 / 2.0;
    let delta = delta / 2.0;

    // Value passed to the Beta distribution function
    let tmp = nu1 * x / (nu2 + nu1 * x);
    let logtmp = tmp.ln();
    let nu2const = nu2 * (1.0 - tmp).ln() - gammaln(nu2);

    // Start from j0 = floor(delta)
    let j0 = delta.floor();

    // Poisson pdf and beta cdf at the starting point
    let bcdf0 = betainc(j0 + nu1, nu2, tmp, upper);
    let ppdf0 = (-delta + j0 * delta.ln() - gammaln(j0 + 1.0)).exp();

    // Downward loop
    let mut y = ppdf0 * bcdf0;
    let mut ppdf = ppdf0;
    let mut bcdf = bcdf0;
    let mut olddy = 0.0;
    let mut j = j0 - 1.0;
    while j >= 0.0 {
        ppdf = ppdf * (j + 1.0) / delta;
        if upper {
            bcdf = betainc(j + nu1, nu2, tmp, true);
        } else {
            let db = ((j + nu1) * logtmp
                + nu2const
                + gammaln(j + nu1 + nu2)
                - gammaln(j + nu1 + 1.0))
            .exp();
            bcdf += db;
        }
        let deltay = ppdf * bcdf;
        y += deltay;

        // Convergence: change small and not increasing
        if !(deltay > y * cvgtol || deltay.abs() > olddy) {
            break;
        }
        olddy = deltay.abs();
        j -= 1.0;
    }

    // Upward loop
    ppdf = ppdf0;
    bcdf = bcdf0;
    olddy = 0.0;
    j = j0 + 1.0;
    let mut jj = 0;
    for _ in 0..5000 {
        jj += 1;
        ppdf = ppdf * delta / j;
        if upper {
            bcdf = betainc(j + nu1, nu2, tmp, true);
        } else {
            let sub = ((j + nu1 - 1.0) * logtmp
                + nu2const
                + gammaln(j + nu1 + nu2 - 1.0)
                - gammaln(j + nu1))
            .exp();
            bcdf -= sub;
        }
        let deltay = ppdf * bcdf;
        y += deltay;

        // Convergence: change small and not increasing
        if !(deltay > y * cvgtol || deltay.abs() > olddy) {
            break;
        }
        olddy = deltay.abs();
        j += 1.0;
    }

    if jj == 5000 {
        // No convergence — MATLAB issues a warning; here we just return.
    }

    y
}