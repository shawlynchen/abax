use crate::gamma;
use std::f64::consts::PI;

/// Compute the Bessel function of the first kind <math><msub><mi>J</mi><mi>v</mi></msub><mo>(</mo><mi>x</mi><mo>)</mo></math> for non-negative real input.
///
/// The Bessel function of the first kind is defined as a solution to the
/// Bessel differential equation:
/// <math display="block">
///   <msup><mi>x</mi><mn>2</mn></msup><mfrac><mrow><msup><mi>d</mi><mn>2</mn></msup><mi>y</mi></mrow><mrow><mi>d</mi><msup><mi>x</mi><mn>2</mn></msup></mrow></mfrac>
///   <mo>+</mo><mi>x</mi><mfrac><mrow><mi>dy</mi></mrow><mrow><mi>dx</mi></mrow></mfrac>
///   <mo>+</mo><mo>(</mo><msup><mi>x</mi><mn>2</mn></msup><mo>-</mo><msup><mi>ν</mi><mn>2</mn></msup><mo>)</mo><mi>y</mi><mo>=</mo><mn>0</mn>
/// </math>
///
/// This implementation uses a combination of strategies:
/// - **Power Series**: Used for small arguments (<math><mi>x</mi><mo>&lt;</mo><mn>45</mn></math>).
/// - **Asymptotic Expansion**: Used for large arguments (<math><mi>x</mi><mo>≥</mo><mn>45</mn></math>).
///
/// # Arguments
/// * `nu` - Order of the Bessel function (must be non-negative).
/// * `x` - Real argument (must be non-negative).
/// * `scale` - If `true`, returns the exponentially scaled function <math><msub><mi>J</mi><mi>ν</mi></msub><mo>(</mo><mi>x</mi><mo>)</mo><mo>⋅</mo><msup><mi>e</mi><mrow><mo>-</mo><mi>x</mi></mrow></msup></math>
///   to prevent numerical overflow for large arguments.
///
/// # Domain
/// - Returns `NaN` if `nu < 0` or `x < 0`.
///
/// # Examples
/// ```
/// use abax::besselj;
///
/// // J_0(1.0) is approximately 0.7651976865
/// let val = besselj(0.0, 1.0, false);
/// assert!((val - 0.7651976865579666).abs() < 1e-15);
///
/// // J_1(1.0) is approximately 0.4400505857
/// let val = besselj(1.0, 1.0, false);
/// assert!((val - 0.4400505857449335).abs() < 1e-15);
/// ```
pub fn besselj(nu: f64, x: f64, scale: bool) -> f64 {
    if nu < 0.0 || x < 0.0 || nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }

    // Handle edge-case at x = 0
    if x == 0.0 {
        return if nu == 0.0 { 1.0 } else { 0.0 };
    }

    // Direct crossover threshold based on John Harrison's method
    if x < 45.0 {
        besselj_small(nu, x, scale)
    } else {
        besselj_large(nu, x, scale)
    }
}

/// Evaluation for small/moderate arguments using the power series expansion
fn besselj_small(nu: f64, x: f64, scale: bool) -> f64 {
    let mut sum = 0.0;
    let mut term = 1.0;
    
    // Initial leading term coefficient: (x / 2)^nu / Gamma(nu + 1)
    let leading_pow = (x / 2.0).powf(nu);
    let gamma_val = gamma(nu + 1.0);
    
    if gamma_val.is_infinite() || leading_pow.is_infinite() {
        return 0.0;
    }
    
    let mut init_coef = leading_pow / gamma_val;
    if scale {
        init_coef *= (-x).exp(); // Inject scaling directly to prevent underflow loops
    }

    let x_half_sq = (x / 2.0) * (x / 2.0);
    
    for m in 0..200 {
        if m > 0 {
            // Factorial-step recurrence update
            term *= -x_half_sq / (m as f64 * (nu + m as f64));
        }
        
        let old_sum = sum;
        sum += term;
        
        // Break early when precision reaches machine epsilon limits
        if (sum - old_sum).abs() < f64::EPSILON * sum.abs() {
            break;
        }
    }

    init_coef * sum
}

/// Evaluation for large arguments via asymptotic phase & modulus expansion
fn besselj_large(nu: f64, x: f64, scale: bool) -> f64 {
    let (alpha, beta) = if nu < 1e-9 {
        // Optimized asymptotic coefficients for order nu = 0
        let inv_x2 = 1.0 / (x * x);
        let a = 1.0 / (8.0 * x) - (25.0 / 384.0) * inv_x2 / x + (1073.0 / 5120.0) * (inv_x2 * inv_x2) / x;
        let b = 1.0 - (1.0 / 16.0) * inv_x2 + (53.0 / 512.0) * (inv_x2 * inv_x2);
        (a, b)
    } else if (nu - 1.0).abs() < 1e-9 {
        // Optimized asymptotic coefficients for order nu = 1
        let inv_x2 = 1.0 / (x * x);
        let a = -3.0 / (8.0 * x) + (21.0 / 128.0) * inv_x2 / x - (1899.0 / 5120.0) * (inv_x2 * inv_x2) / x;
        let b = 1.0 + (3.0 / 16.0) * inv_x2 - (99.0 / 512.0) * (inv_x2 * inv_x2);
        (a, b)
    } else {
        // General calculation for arbitrary positive real orders
        let mu = 4.0 * nu * nu;
        
        let p = 1.0 - (mu - 1.0) * (mu - 9.0) / (2.0 * (8.0 * x).powi(2));
        let q = (mu - 1.0) / (8.0 * x) - (mu - 1.0) * (mu - 9.0) * (mu - 25.0) / (6.0 * (8.0 * x).powi(3));
        
        let b = (p * p + q * q).sqrt();
        let a = q.atan2(p);
        (a, b)
    };

    // Calculate oscillatory phase argument
    let phase = x - (nu / 2.0 + 0.25) * PI - alpha;
    
    // Construct amplitude envelope
    let mut amplitude = (2.0 / (PI * x)).sqrt() * beta;
    if scale {
        amplitude *= (-x).exp(); 
    }

    amplitude * phase.cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_besselj_boundaries() {
        assert_eq!(besselj(0.0, 0.0, false), 1.0);
        assert_eq!(besselj(1.0, 0.0, false), 0.0);
        assert_eq!(besselj(2.5, 0.0, false), 0.0);
        assert!(besselj(-1.0, 1.0, false).is_nan());
        assert!(besselj(1.0, -1.0, false).is_nan());
    }

    #[test]
    fn test_besselj_small_x() {
        let tol = 1e-15;
        // J_0(1.0) ref: 0.7651976865579666
        assert!((besselj(0.0, 1.0, false) - 0.7651976865579666).abs() < tol);
        // J_1(1.0) ref: 0.4400505857449335
        assert!((besselj(1.0, 1.0, false) - 0.4400505857449335).abs() < tol);
    }

    #[test]
    fn test_besselj_half_integer() {
        let tol = 1e-15;
        // J_0.5(x) = sqrt(2 / (pi * x)) * sin(x)
        let x = 1.0;
        let expected = (2.0 / (PI * x)).sqrt() * x.sin();
        assert!((besselj(0.5, x, false) - expected).abs() < tol);
    }

    #[test]
    fn test_besselj_large_x() {
        let tol = 1e-10; // Asymptotic approximation threshold
        let x = 50.0;
        let val = besselj(0.0, x, false);
        assert!((val - (0.05581232766925208)).abs() < tol);
        
        let val1 = besselj(1.0, x, false);
        assert!((val1 - (-9.75118281251751e-2)).abs() < tol);
    }

    #[test]
    fn test_besselj_scaling() {
        let x = 1.0;
        let val = besselj(0.0, x, false);
        let val_scaled = besselj(0.0, x, true);
        assert!((val_scaled - val * (-x).exp()).abs() < 1e-15);
    }
}

