use crate::rutils::{ISNAN, R_FINITE};
use crate::consts::{ML_NEGINF, ML_POSINF, M_LN_SQRT_2PI};
use crate::gammaln::{lgammacor, lgammafn, lgamma};
use crate::gamma::gammafn;

/// Natural logarithm of the Beta function.
///
/// See also [`lbeta`](./fn.lbeta.html).
/// 
/// The Beta function is defined as:
/// 
/// <math display="block" xmlns="http://www.w3.org/1998/Math/MathML">
///   <mi>B</mi><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo>
///   <mo>=</mo>
///   <mfrac>
///     <mrow><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>)</mo><mi>Γ</mi><mo>(</mo><mi>w</mi><mo>)</mo></mrow>
///     <mrow><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>+</mo><mi>w</mi><mo>)</mo></mrow>
///   </mfrac>
/// </math>
///
/// This function computes <math xmlns="http://www.w3.org/1998/Math/MathML">
///   <mi>ln</mi><mo>⁡</mo><mi>B</mi><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo>
/// </math> using the logarithmic Gamma function
/// to maintain numerical stability and avoid overflow/underflow.
/// 
/// <math display="block" xmlns="http://www.w3.org/1998/Math/MathML">
///   <mi>ln</mi><mo></mo><mi>B</mi><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo>
///   <mo>=</mo>
///   <mi>ln</mi><mo>⁡</mo><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>)</mo>
///   <mo>+</mo>
///   <mi>ln</mi><mo>⁡</mo><mi>Γ</mi><mo>(</mo><mi>w</mi><mo>)</mo>
///   <mo>−</mo>
///   <mi>ln</mi><mo>⁡</mo><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>+</mo><mi>w</mi><mo>)</mo>
/// </math>
pub fn betaln(z: f64, w: f64) -> f64 {
    lbeta(z, w)
}

/// See [`betaln`](./fn.betaln.html) for details.
pub fn lbeta(a: f64, b: f64) -> f64 {
    if ISNAN(a) || ISNAN(b) {
	    return f64::NAN;
    }

    let p = f64::min(a, b);
    let q = f64::max(a, b);

    /* both arguments must be >= 0 */
    if p < 0.0 {
        return f64::NAN;
    } else if p == 0.0 {
	    return ML_POSINF;
    } else if !R_FINITE(q) { /* q == +Inf */
	    return ML_NEGINF;
    }

    if p >= 10.0 {
	    /* p and q are big. */
	    let corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	    return f64::ln(q) * - 0.5 + M_LN_SQRT_2PI + corr + (p - 0.5) * f64::ln(p / (p + q)) + q * f64::ln_1p(-p / (p + q));
    } else if q >= 10.0 {
	    /* p is small, but q is big. */
	    let corr = lgammacor(q) - lgammacor(p + q);
	    return lgammafn(p) + corr + p - p * f64::ln(p + q) + (q - 0.5) * f64::ln_1p(-p / (p + q));
    } else {
	    /* p and q are small: p <= q < 10. */
	    if p < 1e-306 {
            return lgamma(p) + (lgamma(q) - lgamma(p+q));
        } else {
            return f64::ln(gammafn(p) * (gammafn(q) / gammafn(p + q)));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_betaln_values() {
        // B(1, 1) = 1, ln(1) = 0
        assert!((betaln(1.0, 1.0)).abs() < 1e-14);
        // B(2, 1) = 1/2, ln(0.5) = -0.6931471805599453
        assert!((betaln(2.0, 1.0) - (-0.6931471805599453)).abs() < 1e-14);
        // Symmetry: B(z, w) = B(w, z)
        assert!((betaln(2.5, 3.5) - betaln(3.5, 2.5)).abs() < 1e-14);
    }

    #[test]
    fn test_betaln_domain() {
        assert!(betaln(0.0, 1.0).is_nan());
        assert!(betaln(1.0, -1.0).is_nan());
        assert!(betaln(f64::NAN, 1.0).is_nan());
    }
}
