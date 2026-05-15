use crate::{betainc, gammaln, normcdf, tcdf};

/// Noncentral T cumulative distribution function (CDF).
///
/// Returns the cumulative distribution function for the noncentral T distribution
/// with `nu` degrees of freedom and noncentrality parameter `delta` at the value `x`.
///
/// # Mathematical Definition
/// The noncentral T distribution is the distribution of the random variable:
/// <math display="block">
///   <mi>T</mi>
///   <mo>=</mo>
///   <mfrac>
///     <mrow>
///       <mi>Z</mi>
///       <mo>+</mo>
///       <mi>δ</mi>
///     </mrow>
///     <msqrt>
///       <mi>V</mi>
///       <mo>/</mo>
///       <mi>ν</mi>
///     </msqrt>
///   </mfrac>
/// </math>
/// where <math><mi>Z</mi></math> is a standard normal distribution and <math><mi>V</mi></math> is a chi-squared distribution 
/// with <math><mi>ν</mi></math> degrees of freedom.
///
/// # Domain
/// - `nu > 0`
/// - Returns `NaN` if any input is `NaN` or `nu <= 0`.
///
/// # Examples
/// ```
/// use abax::nctcdf;
///
/// // delta = 0 reduces to the central Student's T distribution
/// let p = nctcdf(0.0, 5.0, 0.0, false);
/// assert!((p - 0.5).abs() < 1e-15);
/// ```
pub fn nctcdf(x: f64, nu: f64, delta: f64, upper: bool) -> f64 {
    // ------------------------------------------------------------------
    // Validation
    // ------------------------------------------------------------------

    if x.is_nan() || nu.is_nan() || delta.is_nan() {
        return f64::NAN;
    }

    if nu <= 0.0 || delta.is_infinite() {
        return f64::NAN;
    }

    // Central t special case
    if delta == 0.0 {
        return tcdf(x, nu, upper);
    }

    // Infinite x
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }

    if x == f64::NEG_INFINITY {
        return if upper { 1.0 } else { 0.0 };
    }

    // Large-nu normal approximation
    // Johnson & Kotz eq 26.7.10
    if nu > 2e6 {
        let s = 1.0 - 1.0 / (4.0 * nu);
        let d = (1.0 + x * x / (2.0 * nu)).sqrt();
        return normcdf(x * s, delta, d, upper);
    }

    // Symmetry relation:
    // F(-x;nu,delta) = 1 - F(x;nu,-delta)
    if x < 0.0 {
        return nctcdf(-x, nu, -delta, !upper);
    }

    // ------------------------------------------------------------------
    // Main computation (x >= 0)
    // ------------------------------------------------------------------

    let xsq = x * x;
    let denom = nu + xsq;

    let p_beta = xsq / denom;
    let q_beta = nu / denom;

    let dsq = delta * delta;
    let signd = delta.signum();

    // Base probability P(T < 0)
    let mut p = if upper {
        if x == 0.0 {
            normcdf(-delta, 0.0, 1.0, true)
        } else {
            0.0
        }
    } else {
        normcdf(-delta, 0.0, 1.0, false)
    };

    if x == 0.0 {
        return p;
    }

    // ------------------------------------------------------------------
    // Start near peak of Poisson weights
    // ------------------------------------------------------------------

    let jj0 = 2.0 * (dsq / 2.0).floor();

    // E terms
    //
    // exp(0.5*j*log(0.5*dsq) - dsq/2 - gammaln(j/2+1))
    //
    // Handle extremely tiny dsq robustly.
    let log_half_dsq = if dsq > 0.0 {
        (0.5 * dsq).ln()
    } else {
        f64::NEG_INFINITY
    };

    let e1_0 = (0.5 * jj0 * log_half_dsq
        - dsq / 2.0
        - gammaln(jj0 / 2.0 + 1.0))
        .exp();

    let e2_0 = signd
        * (0.5 * (jj0 + 1.0) * log_half_dsq
            - dsq / 2.0
            - gammaln((jj0 + 1.0) / 2.0 + 1.0))
        .exp();

    // ------------------------------------------------------------------
    // Initial B terms
    // ------------------------------------------------------------------

    let use_p = p_beta < 0.5;

    let (b1_0, b2_0) = if upper {
        if use_p {
            (
                betainc(
                    p_beta,
                    (jj0 + 1.0) / 2.0,
                    nu / 2.0,
                    false,
                ),
                betainc(
                    p_beta,
                    (jj0 + 2.0) / 2.0,
                    nu / 2.0,
                    false,
                ),
            )
        } else {
            (
                betainc(
                    q_beta,
                    nu / 2.0,
                    (jj0 + 1.0) / 2.0,
                    true,
                ),
                betainc(
                    q_beta,
                    nu / 2.0,
                    (jj0 + 2.0) / 2.0,
                    true,
                ),
            )
        }
    } else {
        if use_p {
            (
                betainc(
                    p_beta,
                    (jj0 + 1.0) / 2.0,
                    nu / 2.0,
                    true,
                ),
                betainc(
                    p_beta,
                    (jj0 + 2.0) / 2.0,
                    nu / 2.0,
                    true,
                ),
            )
        } else {
            (
                betainc(
                    q_beta,
                    nu / 2.0,
                    (jj0 + 1.0) / 2.0,
                    false,
                ),
                betainc(
                    q_beta,
                    nu / 2.0,
                    (jj0 + 2.0) / 2.0,
                    false,
                ),
            )
        }
    };

    // ------------------------------------------------------------------
    // Initial recurrence ratios
    // ------------------------------------------------------------------

    let r1_0 = (
        gammaln((jj0 + 1.0) / 2.0 + nu / 2.0)
            - gammaln((jj0 + 3.0) / 2.0)
            - gammaln(nu / 2.0)
            + ((jj0 + 1.0) / 2.0) * p_beta.ln()
            + (nu / 2.0) * q_beta.ln()
    )
        .exp();

    let r2_0 = (
        gammaln((jj0 + 2.0) / 2.0 + nu / 2.0)
            - gammaln((jj0 + 4.0) / 2.0)
            - gammaln(nu / 2.0)
            + ((jj0 + 2.0) / 2.0) * p_beta.ln()
            + (nu / 2.0) * q_beta.ln()
    )
        .exp();

    let mut subtotal = 0.0;

    // ------------------------------------------------------------------
    // Upward sweep
    // ------------------------------------------------------------------

    let mut jj = jj0;

    let mut e1 = e1_0;
    let mut e2 = e2_0;

    let mut b1 = b1_0;
    let mut b2 = b2_0;

    let mut r1 = r1_0;
    let mut r2 = r2_0;

    loop {
        let twoterms = e1 * b1 + e2 * b2;

        subtotal += twoterms;

        // Convergence
        if twoterms.abs() <= subtotal.abs() * 1e-15 {
            break;
        }

        jj += 2.0;

        e1 *= dsq / jj;
        e2 *= dsq / (jj + 1.0);

        if upper {
            b1 = betainc(
                p_beta,
                (jj + 1.0) / 2.0,
                nu / 2.0,
                false,
            );

            b2 = betainc(
                p_beta,
                (jj + 2.0) / 2.0,
                nu / 2.0,
                false,
            );
        } else {
            b1 -= r1;
            b2 -= r2;

            r1 *= p_beta * (jj + nu - 1.0) / (jj + 1.0);

            r2 *= p_beta * (jj + nu) / (jj + 2.0);
        }

        if jj > 100000.0 {
            break;
        }
    }

    // ------------------------------------------------------------------
    // Downward sweep
    // ------------------------------------------------------------------

    jj = jj0;

    e1 = e1_0;
    e2 = e2_0;

    b1 = b1_0;
    b2 = b2_0;

    r1 = r1_0;
    r2 = r2_0;

    while jj > 0.0 {
        let jj_cur = jj;

        e1 *= jj_cur / dsq;
        e2 *= (jj_cur + 1.0) / dsq;

        r1 *= (jj_cur + 1.0)
            / ((jj_cur + nu - 1.0) * p_beta);

        r2 *= (jj_cur + 2.0)
            / ((jj_cur + nu) * p_beta);

        if upper {
            b1 = betainc(
                p_beta,
                (jj_cur - 1.0) / 2.0,
                nu / 2.0,
                false,
            );

            b2 = betainc(
                p_beta,
                jj_cur / 2.0,
                nu / 2.0,
                false,
            );
        } else {
            b1 += r1;
            b2 += r2;
        }

        let twoterms = e1 * b1 + e2 * b2;

        subtotal += twoterms;

        if twoterms.abs() <= subtotal.abs() * 1e-15 {
            break;
        }

        jj -= 2.0;
    }

    // ------------------------------------------------------------------
    // Final probability
    // ------------------------------------------------------------------

    if upper {
        let p_gt_0 = normcdf(-delta, 0.0, 1.0, true);

        p = p_gt_0 - subtotal / 2.0;
    } else {
        p += subtotal / 2.0;
    }

    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nctcdf_central() {
        let x = 1.5;
        let nu = 10.0;
        // Central T (delta=0)
        let p_nct = nctcdf(x, nu, 0.0, false);
        let p_t = tcdf(x, nu, false);
        assert!((p_nct - p_t).abs() < 1e-15);
    }

    #[test]
    fn test_nctcdf_known_values() {
        let tol = 1e-10;
        println!("{:.17}", nctcdf(2.0, 5.0, 1.0, false));
        println!("{:.17}", nctcdf(2.0, 10.0, 2.0, false));
        // Comparison with standard statistical tables/software (e.g., R nctcdf)
        // nu=5, delta=1, x=2
        assert!((nctcdf(2.0, 5.0, 1.0, false) - 7.7807466261621487e-1).abs() < tol);
        // nu=10, delta=2, x=2
        assert!((nctcdf(2.0, 10.0, 2.0, false) - 4.8097315281790721e-1).abs() < tol);
    }

    #[test]
    fn test_nctcdf_symmetry() {
        let x = 1.0;
        let nu = 8.0;
        let delta = 1.5;
        // F(x, v, delta) = 1 - F(-x, v, -delta)
        let p1 = nctcdf(x, nu, delta, false);
        let p2 = nctcdf(-x, nu, -delta, true);
        assert!((p1 - p2).abs() < 1e-14);
    }

    #[test]
    fn test_nctcdf_limits() {
        assert_eq!(nctcdf(f64::INFINITY, 5.0, 1.0, false), 1.0);
        assert_eq!(nctcdf(f64::NEG_INFINITY, 5.0, 1.0, false), 0.0);
        assert!(nctcdf(1.0, 0.0, 1.0, false).is_nan());
    }
}
