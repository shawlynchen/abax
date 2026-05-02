use crate::consts::BERNOULLI_EVEN;

/// Computes the digamma function, <math><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>)</mo></math>,
/// the logarithmic derivative of the Gamma function:
/// <math><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mfrac><mrow><mi>d</mi></mrow><mrow><mi>d</mi><mi>x</mi></mrow></mfrac><mi>ln</mi><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo></math>.
///
/// # Numerical stability
/// This implementation combines three stable strategies:
/// - Reflection for negative values: <math><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mi>ψ</mi><mo>(</mo><mn>1</mn><mo>-</mo><mi>x</mi><mo>)</mo><mo>-</mo><mi>π</mi><mi>cot</mi><mo>(</mo><mi>π</mi><mi>x</mi><mo>)</mo></math>.
/// - Recurrence for small positive values: <math><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>=</mo><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>+</mo><mn>1</mn><mo>)</mo><mo>-</mo><mfrac><mn>1</mn><mi>x</mi></mfrac></math>.
/// - Asymptotic expansion for large values:
///   <math><mi>ψ</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>≈</mo><mi>ln</mi><mi>x</mi><mo>-</mo><mfrac><mn>1</mn><mrow><mn>2</mn><mi>x</mi></mrow></mfrac><mo>-</mo><munderover><mo>∑</mo><mrow><mi>n</mi><mo>=</mo><mn>1</mn></mrow><mi>∞</mi></munderover><mfrac><msub><mi>B</mi><mrow><mn>2</mn><mi>n</mi></mrow></msub><mrow><mn>2</mn><mi>n</mi><msup><mi>x</mi><mrow><mn>2</mn><mi>n</mi></mrow></msup></mrow></mfrac></math>.
///
/// # Special cases
/// - Returns `f64::NAN` for `NaN`.
/// - Returns `f64::NAN` at non-positive integer poles.
/// - Returns `f64::INFINITY` for `+∞`.
///
/// # Examples
/// ```
/// use abax::digamma;
///
/// assert!((digamma(1.0) + 0.5772156649015329).abs() < 1e-14);
/// assert!((digamma(0.5) + 1.9635100260214235).abs() < 1e-14);
/// ```
pub fn digamma(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            f64::NAN
        };
    }

    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }

    if x < 0.5 {
        let pix = std::f64::consts::PI * x;
        return digamma(1.0 - x) - std::f64::consts::PI * (pix.cos() / pix.sin());
    }

    let mut xx = x;
    let mut acc = 0.0;

    while xx < 8.0 {
        acc -= 1.0 / xx;
        xx += 1.0;
    }

    let inv = 1.0 / xx;
    let inv2 = inv * inv;

    let mut series = xx.ln() - 0.5 * inv;
    let mut p = inv2;

    for (k, &b2k) in BERNOULLI_EVEN.iter().enumerate() {
        let denom = 2.0 * (k as f64 + 1.0);
        series -= b2k * p / denom;
        p *= inv2;
        if (b2k * p / denom).abs() < 1e-18 {
            break;
        }
    }

    acc + series
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_approx_eq(actual: f64, expected: f64, eps: f64) {
        let d = (actual - expected).abs();
        assert!(
            d < eps,
            "actual={} expected={} diff={} eps={}",
            actual,
            expected,
            d,
            eps
        );
    }

    #[test]
    fn test_special_cases() {
        assert!(digamma(f64::NAN).is_nan());
        assert_eq!(digamma(f64::INFINITY), f64::INFINITY);
        assert!(digamma(-1.0).is_nan());
        assert!(digamma(0.0).is_nan());
    }

    #[test]
    fn test_known_values() {
        assert_approx_eq(digamma(1.0), -0.5772156649015329, 1e-14);
        assert_approx_eq(digamma(0.5), -1.9635100260214235, 1e-14);
        assert_approx_eq(digamma(5.0), 1.5061176684318003, 1e-14);
    }

    #[test]
    fn test_recurrence() {
        let x = 2.75;
        let lhs = digamma(x + 1.0);
        let rhs = digamma(x) + 1.0 / x;
        assert_approx_eq(lhs, rhs, 1e-14);
    }

    #[test]
    fn test_reflection_noninteger_negative() {
        assert_approx_eq(digamma(-0.5), 0.03648997397857652, 1e-14);
    }
}
