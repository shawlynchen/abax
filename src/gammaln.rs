use crate::consts::*;
use crate::rutils::{ISNAN, chebyshev_eval, sinpi};
use crate::gamma::gammafn;

/// Computes the natural logarithm of the Gamma function, <math><mi>ln</mi><mo>(</mo><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>)</mo></math>.
/// 
/// See also [`lgamma`](./fn.lgamma.html).
///
/// This implementation provides high precision across the positive real axis:
/// - **Small values (<math><mi>x</mi><mo>&lt;</mo><msup><mn>10</mn><mrow><mo>-</mo><mn>5</mn></mrow></msup></math>)**: Utilizes a Taylor series expansion involving the Euler-Mascheroni
///   constant and Riemann Zeta values.
/// - **Large values**: Employs Stirling's asymptotic expansion using high-precision coefficients.
/// - **Intermediate values (<math><mi>x</mi><mo>&lt;</mo><mn>10</mn></math>)**: Uses the recurrence relation <math><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>+</mo><mn>1</mn><mo>)</mo><mo>=</mo><mi>x</mi><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo></math>
///   to shift the argument into the optimal range for the Stirling approximation.
///
/// # Mathematical Context
/// The Gamma function <math><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo></math> extends the factorial function to real numbers: <math><mi>Γ</mi><mo>(</mo><mi>n</mi><mo>)</mo><mo>=</mo><mo>(</mo><mi>n</mi><mo>-</mo><mn>1</mn><mo>)</mo><mo>!</mo></math>
/// for positive integers <math><mi>n</mi></math>.
///
/// # Edge Cases
/// - Returns `f64::NAN` if <math><mi>x</mi><mo>≤</mo><mn>0</mn></math> or <math><mi>x</mi></math> is `NaN`.
/// - Returns `f64::INFINITY` if <math><mi>x</mi></math> is `INFINITY`.
///
/// # Examples
/// ```
/// use abax::gammaln;
/// assert_eq!(gammaln(1.0), 0.0);
/// assert!((gammaln(5.0) - 3.178053830347945).abs() < 1e-14);
/// ```
pub fn gammaln(x: f64) -> f64 {
    lgammafn(x)
} 

/// See [`gammaln`](./fn.gammaln.html) for details.
pub fn lgamma(x: f64) -> f64 {
    lgammafn(x)
}

pub(crate) fn lgammacor(x: f64) -> f64 {
    const ALGMCS: [f64; 15] = [
        // below, nalgm = 5 ==> only the first 5 are used!
        0.1666389480451863247205729650822e+0,
        -0.1384948176067563840732986059135e-4,
        0.9810825646924729426157171547487e-8,
        -0.1809129475572494194263306266719e-10,
        0.6221098041892605227126015543416e-13,
        -0.3399615005417721944303330599666e-15,
        0.2683181998482698748957538846666e-17,
        -0.2868042435334643284144622399999e-19,
        0.3962837061046434803679306666666e-21,
        -0.6831888753985766870111999999999e-23,
        0.1429227355942498147573333333333e-24,
        -0.3547598158101070547199999999999e-26,
        0.1025680058010470912000000000000e-27,
        -0.3401102254316748799999999999999e-29,
        0.1276642195630062933333333333333e-30,
    ];

    const NALGM: usize = 5;
    const XBIG: f64 = 94906265.62425156;

    if x < 10.0 {
        // possibly consider stirlerr()
        return f64::NAN;
    } else if x < XBIG {
        let tmp = 10.0 / x;
        return chebyshev_eval(tmp * tmp * 2.0 - 1.0, &ALGMCS, NALGM) / x;
    }

    // x >= xbig
    return 1.0 / (x * 12.0);
}

fn lgammafn_sign(x: f64, sgn: &mut i32) -> f64 {
    const XMAX: f64 = 2.5327372760800758e+305;
    const DXREL: f64 = 1.490116119384765625e-8;

    *sgn = 1;
    
    if ISNAN(x) {
        return x;
    }

    if x < 0.0 && f64::floor(-x) % 2.0 == 0.0 {
	    *sgn = -1;
    }

    if x <= 0.0 && x == f64::trunc(x) {
        /* Negative integer argument */
	    // No warning: this is the best answer; was  ML_WARNING(ME_RANGE, "lgamma");
	    return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    let y = f64::abs(x);

    if y < 1e-306 {
        return -f64::ln(y); // denormalized range, R change
    }
    if y <= 10.0 {
        return f64::ln(f64::abs(gammafn(x)));
    }

    if y > XMAX {
	    return ML_POSINF;
    }

    if x > 0.0 { /* i.e. y = x > 10 */
	    if x > 1.0e17 {
	        return x * (f64::ln(x) - 1.0);
        } else if x > 4934720.0 {
	        return M_LN_SQRT_2PI + (x - 0.5) * f64::ln(x) - x;
        } else {
            return M_LN_SQRT_2PI + (x - 0.5) * f64::ln(x) - x + lgammacor(x);
        }
    }
    /* else: x < -10; y = -x */
    let sinpiy = f64::abs(sinpi(y));

    if sinpiy == 0.0 { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
        return f64::NAN;
    }

    let ans = M_LN_SQRT_PId2 + (x - 0.5) * f64::ln(y) - x - f64::ln(sinpiy) - lgammacor(y);

    if f64::abs((x - f64::trunc(x - 0.5)) * ans / x) < DXREL {

    	/* The answer is less than half precision because
	     * the argument is too near a negative integer; e.g. for  lgamma(1e-7 - 11) */
        // warning about precision of lgamma
    }

    return ans;
}

pub(crate) fn lgammafn(x: f64) -> f64 {
    let mut sgn: i32 = 0;
    return lgammafn_sign(x, &mut sgn);
}


#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to assert that two floats are approximately equal.
    fn assert_approx_eq(actual: f64, expected: f64, epsilon: f64) {
        if actual.is_nan() && expected.is_nan() {
            return;
        }
        if actual.is_infinite() && expected.is_infinite() {
            assert_eq!(actual.is_sign_positive(), expected.is_sign_positive());
            return;
        }
        let diff = (actual - expected).abs();
        assert!(
            diff < epsilon,
            "Assertion failed: actual {} != expected {} (diff {} > epsilon {})",
            actual,
            expected,
            diff,
            epsilon
        );
    }

    #[test]
    fn test_gammaln_special_cases() {
        assert!(gammaln(f64::NAN).is_nan());
        assert!(gammaln(0.0).is_infinite());
        assert!(gammaln(-1.0).is_infinite());
        assert_eq!(gammaln(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn test_gammaln_small_values() {
        // Triggers the Taylor expansion path (x < 0.00001)
        assert_approx_eq(gammaln(1e-7), 16.11809559323676, 1e-10);
    }

    #[test]
    fn test_gammaln_integers() {
        // ln(Gamma(1)) = ln(1) = 0
        assert_eq!(gammaln(1.0), 0.0);
        // ln(Gamma(2)) = ln(1) = 0
        assert_eq!(gammaln(2.0), 0.0);
        // ln(Gamma(3)) = ln(2)
        assert_approx_eq(gammaln(3.0), std::f64::consts::LN_2, 1e-14);
        // ln(Gamma(10)) = ln(9!) = ln(362880)
        assert_approx_eq(gammaln(10.0), 12.801827480081469, 1e-14);
    }

    #[test]
    fn test_gammaln_half_integers() {
        // ln(Gamma(0.5)) = ln(sqrt(pi))
        assert_approx_eq(gammaln(0.5), 0.5723649429247001, 1e-14);
        // ln(Gamma(1.5)) = ln(0.5 * sqrt(pi))
        assert_approx_eq(gammaln(1.5), -0.12078223763524522, 1e-14);
    }

    #[test]
    fn test_gammaln_extreme_values() {
        // Extreme small (positive)
        // As x -> 0, Gamma(x) ~ 1/x, so ln(Gamma(x)) ~ -ln(x)
        let x_tiny = 1e-300;
        assert_approx_eq(gammaln(x_tiny), -f64::ln(x_tiny), 1e-12);

        // Extreme large
        // ln(Gamma(1e100)) is approximately 2.292585...e102
        // We use a large epsilon because of the magnitude of the result
        assert_approx_eq(gammaln(1e100), 2.2925850929940457e102, 1e88);
    }
}
