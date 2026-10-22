use crate::{finv, ncfcdf, ncfpdf};

/// Inverse of the noncentral F cumulative distribution function.
///
/// `p`: probability in (0, 1) for which to invert the cdf
/// `nu1`: numerator degrees of freedom
/// `nu2`: denominator degrees of freedom
/// `delta`: noncentrality parameter
///
/// Returns `f64::NAN` for out-of-range parameters or probabilities.
pub fn ncfinv(p: f64, nu1: f64, nu2: f64, delta: f64) -> f64 {
    let crit = f64::EPSILON.sqrt();
    let count_limit = 100;

    // NaN handling
    if p.is_nan() || nu1.is_nan() || nu2.is_nan() || delta.is_nan() {
        return f64::NAN;
    }

    // Central distribution: delegate to finv
    if delta == 0.0 {
        return finv(p, nu1, nu2);
    }

    let ok_params = nu1 > 0.0 && nu2 > 0.0 && delta >= 0.0;

    // Boundary probabilities
    if !ok_params {
        return f64::NAN;
    }
    if p == 0.0 {
        return 0.0;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    if p < 0.0 || p > 1.0 {
        // p outside [0,1] -> NaN
        return f64::NAN;
    }

    // Start at the mean if it exists (nu2 > 2)
    let mut y = nu2 * (nu1 + delta) / (nu1 * (nu2 - 2.0).max(1.0));

    // Newton's method
    let mut count = 0;
    let mut oldh: f64 = 0.0;
    let mut f_val = ncfcdf(y, nu1, nu2, delta, false);

    while count < count_limit {
        count += 1;

        let fpdf_val = ncfpdf(y, nu1, nu2, delta);
        let mut h = if fpdf_val == 0.0 {
            // Avoid division by zero: no usable Newton step.
            // This can happen in extreme tails; fall back to a small step.
            if f_val > p { -y * 0.5 } else { y * 0.5 }
        } else {
            (f_val - p) / fpdf_val
        };

        // If iterations appear to be oscillating, damp them out
        if h.signum() == -oldh.signum() && oldh != 0.0 {
            h = h.signum() * h.abs().min(oldh.abs()) / 2.0;
        }

        // Avoid stepping too far
        let newy = (y / 5.0).max((5.0 * y).min(y - h));

        // Back off if the step gives a worse result
        let mut newy = newy;
        let mut newf = ncfcdf(newy, nu1, nu2, delta, false);
        loop {
            let worse = (newf - p).abs() > (f_val - p).abs() * (1.0 + crit)
                && (y - newy).abs() > crit * y.abs();
            if !worse {
                break;
            }
            newy = 0.5 * (newy + y);
            newf = ncfcdf(newy, nu1, nu2, delta, false);
        }

        // Convergence check
        h = y - newy;
        let converged = h.abs() <= crit * y.abs();

        // Advance state
        y = newy;
        f_val = newf;
        oldh = h;

        if converged {
            break;
        }
    }

    if count == count_limit {
        // MATLAB issues a warning here; we silently return the best estimate.
    }

    y
}