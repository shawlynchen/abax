use crate::consts::{M_2PI, M_SQRT2, ML_NEGINF, ML_POSINF};
use crate::erfcinv;
use crate::rutils::*;
use std::f64::consts::SQRT_2;

/// Inverse of the normal cumulative distribution function (quantile function).
///
/// Given a probability `p`, a mean `mu`, and a standard deviation `sigma`,
/// this function returns the value `x` such that the probability of a
/// normal random variable being less than or equal to `x` is `p`.
///
/// See also [`qnorm`](./fn.qnorm.html)
///
/// # Mathematical Definition
/// The quantile function is expressed via the inverse error function. For <math><mi>p</mi><mo>∈</mo><mo>(</mo><mn>0</mn><mo>,</mo><mn>1</mn><mo>)</mo></math>:
/// <math display="block"><mi>x</mi><mo>=</mo><mi>μ</mi><mo>+</mo><mi>σ</mi><msqrt><mn>2</mn></msqrt><msup><mi>erf</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>(</mo><mn>2</mn><mi>p</mi><mo>-</mo><mn>1</mn><mo>)</mo></math>
///
/// For a standard normal distribution (<math><mi>μ</mi><mo>=</mo><mn>0</mn><mo>,</mo><mi>σ</mi><mo>=</mo><mn>1</mn></math>), this is the Probit function:
/// <math display="block"><mi>z</mi><mo>=</mo><msup><mi>Φ</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>(</mo><mi>p</mi><mo>)</mo><mo>=</mo><mo>-</mo><msqrt><mn>2</mn></msqrt><msup><mi>erfc</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>(</mo><mn>2</mn><mi>p</mi><mo>)</mo></math>
///
/// # Examples
/// ```
/// use abax::{norminv, qnorm};
///
/// // Median of standard normal is 0
/// assert!((norminv(0.5, 0.0, 1.0) - 0.0).abs() < 1e-15);
/// assert!((qnorm(0.5, 0.0, 1.0, true, false) - 0.0).abs() < 1e-15);
/// // ~1 standard deviation for p=0.8413
/// assert!((norminv(0.8413447460685429, 0.0, 1.0) - 1.0).abs() < 1e-15);
/// assert!((qnorm(0.8413447460685429, 0.0, 1.0, true, false) - 1.0).abs() < 1e-15);
/// ```
pub fn norminv(p: f64, mu: f64, sigma: f64) -> f64 {
    if p.is_nan() || mu.is_nan() || sigma.is_nan() || !(0.0..=1.0).contains(&p) || sigma <= 0.0 {
        return f64::NAN;
    }

    if p == 0.0 {
        return f64::NEG_INFINITY;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }

    // Calculate the standard normal quantile (z-score)
    // Relationship: normcdf(z) = 0.5 * erfc(-z / sqrt(2))
    // Setting p = 0.5 * erfc(-z / sqrt(2)) leads to:
    let z = -SQRT_2 * erfcinv(2.0 * p);

    mu + sigma * z
}

/// See [`norminv`](./fn.norminv.html) for details.
pub fn qnorm(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    qnorm5(p, mu, sigma, lower_tail, log_p)
}

fn qnorm5(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    //double p_, q, r, val;

    if ISNAN(p) || ISNAN(mu) || ISNAN(sigma) {
        return f64::NAN;
    }

    if let Some(q) = R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF, lower_tail, log_p) {
        return q;
    }

    if sigma < 0.0 {
        return f64::NAN;
    }
    if sigma == 0.0 {
        return mu;
    }

    let p_ = R_DT_qIv(p, lower_tail, log_p); /* real lower_tail prob. p */
    let q = p_ - 0.5;

    let mut val: f64;
    if f64::abs(q) <= 0.425 {
        let r = 0.180625 - q * q;
        val = q
            * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r
                + 67265.770927008700853)
                * r
                + 45921.953931549871457)
                * r
                + 13731.693765509461125)
                * r
                + 1971.5909503065514427)
                * r
                + 133.14166789178437745)
                * r
                + 3.387132872796366608)
            / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r
                + 39307.89580009271061)
                * r
                + 21213.794301586595867)
                * r
                + 5394.1960214247511077)
                * r
                + 687.1870074920579083)
                * r
                + 42.313330701600911252)
                * r
                + 1.);
    } else {
        let lp: f64 = if log_p && ((lower_tail && q <= 0.0) || (!lower_tail && q > 0.0)) {
            p
        } else {
            f64::ln(if q > 0.0 {
                R_DT_CIv(p, lower_tail, log_p)
            } else {
                p_
            })
        };

        let mut r: f64 = f64::sqrt(-lp);

        if r <= 5.0 {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r
                + 0.24178072517745061177)
                * r
                + 1.27045825245236838258)
                * r
                + 3.64784832476320460504)
                * r
                + 5.7694972214606914055)
                * r
                + 4.6303378461565452959)
                * r
                + 1.42343711074968357734)
                / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r
                    + 0.0151986665636164571966)
                    * r
                    + 0.14810397642748007459)
                    * r
                    + 0.68976733498510000455)
                    * r
                    + 1.6763848301838038494)
                    * r
                    + 2.05319162663775882187)
                    * r
                    + 1.);
        } else if r <= 27.0 {
            r += -5.0;
            val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r
                + 0.0012426609473880784386)
                * r
                + 0.026532189526576123093)
                * r
                + 0.29656057182850489123)
                * r
                + 1.7848265399172913358)
                * r
                + 5.4637849111641143699)
                * r
                + 6.6579046435011037772)
                / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r
                    + 1.8463183175100546818e-5)
                    * r
                    + 7.868691311456132591e-4)
                    * r
                    + 0.0148753612908506148525)
                    * r
                    + 0.13692988092273580531)
                    * r
                    + 0.59983220655588793769)
                    * r
                    + 1.);
        } else {
            if r >= 6.4e8 {
                val = r * M_SQRT2;
            } else {
                let s2: f64 = -2.0 * lp;
                let mut x2 = s2 - f64::ln(M_2PI * s2);

                if r < 36000.0 {
                    x2 = s2 - f64::ln(M_2PI * x2) - 2.0 / (2.0 + x2);
                    if r < 840.0 {
                        x2 = s2 - f64::ln(M_2PI * x2)
                            + 2.0 * f64::ln_1p(-(1.0 - 1.0 / (4.0 + x2)) / (2.0 + x2));
                        if r < 109.0 {
                            x2 = s2 - f64::ln(M_2PI * x2)
                                + 2.0
                                    * f64::ln_1p(
                                        -(1.0 - (1.0 - 5.0 / (6.0 + x2)) / (4.0 + x2)) / (2.0 + x2),
                                    );
                            if r < 55.0 {
                                x2 = s2 - f64::ln(M_2PI * x2)
                                    + 2.0
                                        * f64::ln_1p(
                                            -(1.0
                                                - (1.0 - (5.0 - 9.0 / (8.0 + x2)) / (6.0 + x2))
                                                    / (4.0 + x2))
                                                / (2.0 + x2),
                                        );
                            }
                        }
                    }
                }
                val = f64::sqrt(x2);
            }
        }
        if q < 0.0 {
            val = -val;
        }
    }

    return mu + sigma * val;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_norminv_standard() {
        let tol = 1e-14;
        assert!((norminv(0.5, 0.0, 1.0) - 0.0).abs() < tol);
        assert!((qnorm(0.5, 0.0, 1.0, true, false) - 0.0).abs() < tol);
        assert!((norminv(0.975, 0.0, 1.0) - 1.959963984540054).abs() < tol);
        assert!((qnorm(0.975, 0.0, 1.0, true, false) - 1.959963984540054).abs() < tol);
        assert!((norminv(0.158655253931457, 0.0, 1.0) + 1.0).abs() < tol);
        assert!((qnorm(0.158655253931457, 0.0, 1.0, true, false) + 1.0).abs() < tol);
    }

    #[test]
    fn test_norminv_boundaries_and_invalid() {
        assert_eq!(norminv(0.0, 0.0, 1.0), f64::NEG_INFINITY);
        assert_eq!(qnorm(0.0, 0.0, 1.0, true, false), f64::NEG_INFINITY);
        assert_eq!(norminv(1.0, 0.0, 1.0), f64::INFINITY);
        assert_eq!(qnorm(1.0, 0.0, 1.0, true, false), f64::INFINITY);
        assert!(norminv(-0.1, 0.0, 1.0).is_nan());
        assert!(qnorm(-0.1, 0.0, 1.0, true, false).is_nan());
        assert!(norminv(1.1, 0.0, 1.0).is_nan());
        assert!(qnorm(1.1, 0.0, 1.0, true, false).is_nan());
        assert!(norminv(0.5, 0.0, -1.0).is_nan());
        assert!(qnorm(0.5, 0.0, -1.0, true, false).is_nan());
        assert!(norminv(0.5, 0.0, 0.0).is_nan());
        assert!(qnorm(0.5, 0.0, 0.0, true, false) == 0.0); // special case: different between qnorm and norminv.
    }
}
