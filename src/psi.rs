use crate::{digamma, gamma, tetragamma, trigamma};

/// Computes the polygamma function `psi(k, x)` (MATLAB-compatible order convention),
/// where `k = 0` is the digamma function and `k >= 1` are higher derivatives.
///
/// Mathematically:
/// <math><msup><mi>ψ</mi><mo>(</mo><mi>k</mi><mo>)</mo></msup><mo>(</mo><mi>x</mi><mo>)</mo></math>.
///
/// # Numerical stability
/// - Uses dedicated implementations for `k = 0, 1, 2`.
/// - For `k >= 3`, uses recurrence shifting to move `x` away from poles and a rapidly convergent
///   positive-series representation.
///
/// # Special cases
/// - Returns `NaN` for `NaN` inputs.
/// - Returns `NaN` at non-positive integer poles.
/// - Returns `+∞` for `k = 0, x = +∞`.
/// - Returns `0.0` for `k >= 1, x = +∞`.
///
/// # Examples
/// ```
/// use abax::psi;
///
/// assert!((psi(0, 1.0) + 0.5772156649015329).abs() < 1e-14);
/// assert!((psi(1, 1.0) - 1.6449340668482264).abs() < 1e-14);
/// assert!((psi(2, 1.0) + 2.4041138063191885).abs() < 1e-13);
/// ```
pub fn psi(k: usize, x: f64) -> f64 {
    match k {
        0 => digamma(x),
        1 => trigamma(x),
        2 => tetragamma(x),
        _ => psi_high_order(k, x),
    }
}

fn psi_high_order(k: usize, x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() { 0.0 } else { f64::NAN };
    }
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }

    let m = k as i32;
    let mut xx = x;
    let mut acc = 0.0;
    let fact = gamma(k as f64 + 1.0);
    if !fact.is_finite() {
        return f64::NAN;
    }

    // Recurrence: ψ^(m)(x) = ψ^(m)(x+1) - (-1)^m m! / x^(m+1)
    while xx < 8.0 {
        let denom = xx.powi(m + 1);
        let step = if k.is_multiple_of(2) {
            -fact / denom
        } else {
            fact / denom
        };
        acc += step;
        xx += 1.0;

        if xx <= 0.0 && xx == xx.floor() {
            return f64::NAN;
        }
    }

    // ψ^(m)(x) = (-1)^(m+1) m! * Σ_{n=0..∞} 1/(x+n)^(m+1), for m>=1
    let mut sum = 0.0;
    let p = (k + 1) as i32;
    for n in 0..200_000usize {
        let t = 1.0 / (xx + n as f64).powi(p);
        sum += t;
        if t < 1e-18 {
            break;
        }
    }

    let sign = if (k + 1).is_multiple_of(2) { 1.0 } else { -1.0 };
    acc + sign * fact * sum
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
        assert!(psi(3, f64::NAN).is_nan());
        assert_eq!(psi(0, f64::INFINITY), f64::INFINITY);
        assert_eq!(psi(3, f64::INFINITY), 0.0);
        assert!(psi(4, 0.0).is_nan());
        assert!(psi(4, -2.0).is_nan());
    }

    #[test]
    fn test_low_order_dispatch() {
        assert_approx_eq(psi(0, 1.0), digamma(1.0), 1e-15);
        assert_approx_eq(psi(1, 1.0), trigamma(1.0), 1e-15);
        assert_approx_eq(psi(2, 1.0), tetragamma(1.0), 1e-14);
    }

    #[test]
    fn test_known_high_order_values() {
        // ψ^(3)(1) = π^4 / 15
        assert_approx_eq(psi(3, 1.0), std::f64::consts::PI.powi(4) / 15.0, 1e-12);
        // ψ^(4)(1) = -24 * ζ(5)
        assert_approx_eq(psi(4, 1.0), -24.88626612344088, 1e-11);
    }

    #[test]
    fn test_recurrence_high_order() {
        let x = 2.75;
        let lhs = psi(3, x + 1.0);
        let rhs = psi(3, x) - 6.0 / x.powi(4);
        assert_approx_eq(lhs, rhs, 1e-12);
    }
}
