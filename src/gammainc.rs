use crate::gammaln;

/// Computes the regularized incomplete gamma function with optional scaling.
///
/// This function returns either lower- or upper-tail regularized incomplete gamma values:
/// - `lower = true` returns the lower regularized incomplete gamma
///   <math><mi>P</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>)</mo></math>.
/// - `lower = false` returns the upper regularized incomplete gamma
///   <math><mi>Q</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>)</mo></math>.
///
/// The returned values are regularized, i.e. divided by
/// <math><mi>Γ</mi><mo>(</mo><mi>a</mi><mo>)</mo></math>.
///
/// If `scaled = true`, the returned tail value is explicitly scaled as:
/// - lower tail: <math><mi>P</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>)</mo><mi>Γ</mi><mo>(</mo><mi>a</mi><mo>+</mo><mn>1</mn><mo>)</mo><msup><mi>e</mi><mi>x</mi></msup><msup><mi>x</mi><mrow><mo>-</mo><mi>a</mi></mrow></msup></math>
/// - upper tail: <math><mi>Q</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>)</mo><mi>Γ</mi><mo>(</mo><mi>a</mi><mo>+</mo><mn>1</mn><mo>)</mo><msup><mi>e</mi><mi>x</mi></msup><msup><mi>x</mi><mrow><mo>-</mo><mi>a</mi></mrow></msup></math>
///
/// where <math><mi>P</mi></math> and <math><mi>Q</mi></math> are the regularized lower and upper tails.
///
/// # Domain
/// - `a` must be positive for finite results.
/// - `x` must be non-negative.
/// - Invalid domain inputs return `NaN`.
pub fn gammainc(x: f64, a: f64, lower: bool, scaled: bool) -> f64 {
    const EPS: f64 = 1e-15;
    const FPMIN: f64 = 1e-300;
    const MAX_IT: usize = 200;

    if x.is_nan() || a.is_nan() || a <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        let base = if lower { 0.0 } else { 1.0 };
        return if scaled {
            if lower { 0.0 } else { f64::INFINITY }
        } else {
            base
        };
    }
    if x.is_infinite() {
        let base = if lower { 1.0 } else { 0.0 };
        return if scaled { 0.0 } else { base };
    }

    let gln = gammaln(a);
    let ax = (a * x.ln() - x - gln).exp();

    let p = if x < a + 1.0 {
        let mut sum = 1.0 / a;
        let mut del = sum;
        let mut ap = a;
        for _ in 0..MAX_IT {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * EPS {
                break;
            }
        }
        sum * ax
    } else {
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / FPMIN;
        let mut d = 1.0 / b;
        let mut h = d;

        for i in 1..=MAX_IT {
            let fi = i as f64;
            let an = -fi * (fi - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < FPMIN {
                d = FPMIN;
            }
            c = b + an / c;
            if c.abs() < FPMIN {
                c = FPMIN;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < EPS {
                break;
            }
        }
        1.0 - ax * h
    };

    let q = 1.0 - p;
    let mut out = if lower { p } else { q };

    if scaled {
        out *= a * (x - a * x.ln() + gln).exp();
    }

    out.clamp(0.0, f64::INFINITY)
}

/// Inverse of the regularized incomplete gamma function.
///
/// Solves for `x >= 0` in either:
/// - `P(a, x) = y` when `lower = true`
/// - `Q(a, x) = y` when `lower = false`
///
/// # Domain
/// - `a > 0`
/// - `0 <= y <= 1`
/// - Invalid inputs return `NaN`.
pub fn gammaincinv(y: f64, a: f64, lower: bool) -> f64 {
    if y.is_nan() || a.is_nan() || a <= 0.0 || !(0.0..=1.0).contains(&y) {
        return f64::NAN;
    }

    if lower {
        if y == 0.0 {
            return 0.0;
        }
        if y == 1.0 {
            return f64::INFINITY;
        }
    } else {
        if y == 1.0 {
            return 0.0;
        }
        if y == 0.0 {
            return f64::INFINITY;
        }
    }

    let target_p = if lower { y } else { 1.0 - y };

    // Bracket solution x in [lo, hi] with monotonicity of P(a, x).
    let mut lo = 0.0;
    let mut hi = a.max(1.0);
    while gammainc(hi, a, true, false) < target_p {
        hi *= 2.0;
        if !hi.is_finite() || hi > 1e308 {
            return f64::INFINITY;
        }
    }

    // Initial guess via interpolation in bracket.
    let mut x = lo + (hi - lo) * target_p;
    if x <= 0.0 {
        x = 0.5 * hi;
    }

    // Safeguarded Newton iterations with bisection fallback.
    let gln = gammaln(a);
    for _ in 0..80 {
        let p = gammainc(x, a, true, false);
        let f = p - target_p;

        if f.abs() <= 2e-14 * target_p.max(1.0 - target_p).max(1e-16) {
            break;
        }

        if f > 0.0 {
            hi = x;
        } else {
            lo = x;
        }

        let pdf = ((a - 1.0) * x.ln() - x - gln).exp();
        let mut x_new = if pdf.is_finite() && pdf > 0.0 {
            x - f / pdf
        } else {
            f64::NAN
        };

        if !x_new.is_finite() || x_new <= lo || x_new >= hi {
            x_new = 0.5 * (lo + hi);
        }

        if (x_new - x).abs() <= 2e-14 * x_new.abs().max(1.0) {
            x = x_new;
            break;
        }
        x = x_new;
    }

    x.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gamma;

    #[test]
    fn test_a1_identities() {
        let x = 2.5;
        let p = gammainc(x, 1.0, true, false);
        let q = gammainc(x, 1.0, false, false);
        assert!((p - (1.0 - (-x).exp())).abs() < 1e-14);
        assert!((q - (-x).exp()).abs() < 1e-14);
        assert!((p + q - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_a_half_erf_relation() {
        let x = 1.0;
        let p = gammainc(x, 0.5, true, false);
        let erf1 = 0.8427007929497149;
        assert!((p - erf1).abs() < 2e-14);
    }

    #[test]
    fn test_scaled_consistency() {
        let x = 15.0;
        let a = 4.5;
        let unscaled = gammainc(x, a, false, false);
        let scaled = gammainc(x, a, false, true);
        let expected = unscaled * gamma(a + 1.0) * x.exp() * x.powf(-a);
        assert!((scaled - expected).abs() / expected.abs() < 1e-12);
    }

    #[test]
    fn test_gammaincinv_roundtrip_lower() {
        let a = 2.5;
        for &y in &[1e-12, 1e-9, 1e-6, 0.01, 0.2, 0.5, 0.8, 0.99, 1.0 - 1e-10] {
            let x = gammaincinv(y, a, true);
            let y2 = gammainc(x, a, true, false);
            assert!((y2 - y).abs() <= 5e-12_f64.max(5e-12 * y.abs()));
        }
    }

    #[test]
    fn test_gammaincinv_roundtrip_upper() {
        let a = 5.0;
        for &y in &[1e-12, 1e-9, 1e-6, 0.01, 0.2, 0.5, 0.8, 0.99, 1.0 - 1e-10] {
            let x = gammaincinv(y, a, false);
            let y2 = gammainc(x, a, false, false);
            assert!((y2 - y).abs() <= 5e-12_f64.max(5e-12 * y.abs()));
        }
    }

    #[test]
    fn test_gammaincinv_domain_and_edges() {
        assert!(gammaincinv(0.5, 0.0, true).is_nan());
        assert!(gammaincinv(-0.1, 1.0, true).is_nan());
        assert_eq!(gammaincinv(0.0, 2.0, true), 0.0);
        assert!(gammaincinv(1.0, 2.0, true).is_infinite());
        assert_eq!(gammaincinv(1.0, 2.0, false), 0.0);
        assert!(gammaincinv(0.0, 2.0, false).is_infinite());
    }
}
