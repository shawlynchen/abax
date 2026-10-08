use crate::gamma::gammafn;
use crate::lbeta;
use crate::rutils::{ISNAN, R_FINITE};
use crate::consts::ML_POSINF;

/// Beta function <math><mrow><mi>B</mi><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo><mo>=</mo><mfrac><mrow><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>)</mo><mi>Γ</mi><mo>(</mo><mi>w</mi><mo>)</mo></mrow><mrow><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>+</mo><mi>w</mi><mo>)</mo></mrow></mfrac></mrow></math>.
///
/// This implementation evaluates the function using the natural logarithm of the Gamma function
/// to avoid intermediate overflow for large arguments:
/// 
/// <math>
///   <mrow><mi>B</mi><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo><mo>=</mo><mi>exp</mi><mo>(</mo><mi>ln</mi><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>)</mo><mo>+</mo><mi>ln</mi><mi>Γ</mi><mo>(</mo><mi>w</mi><mo>)</mo><mo>-</mo><mi>ln</mi><mi>Γ</mi><mo>(</mo><mi>z</mi><mo>+</mo><mi>w</mi><mo>)</mo><mo>)</mo></mrow>
/// </math>
///
/// # Examples
/// ```
/// use abax::beta;
///
/// let result = beta(2.0, 3.0);
/// assert!((result - 1.0/12.0).abs() < 1e-15);
/// ```
pub fn beta(a: f64, b: f64) -> f64 {
    const _XMIN: f64  = -170.5674972726612;
    const XMAX: f64  = 171.61447887182298;
    const _LNSML: f64 = -708.39641853226412; // log(DBL_MIN) = log(2 ^ -1022) = -1022 * log(2)

    /* NaNs propagated correctly */
    if ISNAN(a) || ISNAN(b) {
        return f64::NAN;
    }

    if a < 0.0 || b < 0.0 {
        return f64::NAN;
    } else if a == 0.0 || b == 0.0 {
	    return ML_POSINF;
    } else if !R_FINITE(a) || !R_FINITE(b) {
	    return 0.0;
    }

    if a + b < XMAX {
        //	return gammafn(a) * gammafn(b) / gammafn(a+b);
        /* All the terms are positive, and all can be large for large
        or small arguments.  They are never much less than one.
        gammafn(x) can still overflow for x ~ 1e-308,
        but the result would too.
        */
        let mut r1 = gammafn(a);
        let mut r2 = gammafn(b);
        let mut r3 = f64::recip(gammafn(a + b));
        if r1 > r2 { std::mem::swap(&mut r1, &mut r2); }
        if r2 > r3 {
            std::mem::swap(&mut r2, &mut r3);
            if r1 > r2 { std::mem::swap(&mut r1, &mut r2); }
        }

	    return r1 * r2 * r3;
    } else {
	    let val = lbeta(a, b);
        // underflow to 0 is not harmful per se;  exp(-999) also gives no warning
    	return f64::exp(val);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_beta_values() {
        let tol = 1e-14;
        // B(1, 1) = 1
        assert!((beta(1.0, 1.0) - 1.0).abs() < tol);
        // B(2, 1) = 1/2
        assert!((beta(2.0, 1.0) - 0.5).abs() < tol);
        // B(3, 2) = 1/12
        assert!((beta(3.0, 2.0) - 1.0 / 12.0).abs() < tol);
        // B(0.5, 0.5) = pi
        assert!((beta(0.5, 0.5) - PI).abs() < tol);
    }

    #[test]
    fn test_beta_symmetry() {
        assert_eq!(beta(2.5, 4.2), beta(4.2, 2.5));
    }

    #[test]
    fn test_beta_domain() {
        assert!(beta(0.0, 1.0) == f64::INFINITY);
        assert!(beta(1.0, 0.0) == f64::INFINITY);
        assert!(beta(-0.5, 0.5).is_nan());
        assert!(beta(0.5, -0.5).is_nan());
        assert!(beta(f64::NAN, 1.0).is_nan());
    }
}
