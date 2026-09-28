use crate::consts::{M_LN_SQRT_2PI, ML_POSINF, ML_NEGINF, M_PI};
use crate::rutils::{ISNAN, chebyshev_eval, sinpi};
use crate::gammaln::lgammacor;
use crate::stirlerr::stirlerr;


///Calculates the Gamma function defined by
///
///<math display="block" xmlns="http://www.w3.org/1998/Math/MathML">
///  <mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo>
///  <mo>=</mo>
///  <msubsup><mo>∫</mo><mn>0</mn><mo>∞</mo></msubsup>
///  <msup><mi>t</mi><mrow><mi>x</mi><mo>−</mo><mn>1</mn></mrow></msup>
///  <msup><mi>e</mi><mrow><mo>−</mo><mi>t</mi></mrow></msup>
///  <mi>d</mi><mi>t</mi>
///</math>
///
///for <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>x</mi><mo>></mo><mn>0</mn></math>.
///
///It also satisfies the recurrence
///
///<math display="block" xmlns="http://www.w3.org/1998/Math/MathML">
///  <mi>Γ</mi><mo>(</mo><mi>x</mi><mo>+</mo><mn>1</mn><mo>)</mo>
///  <mo>=</mo>
///  <mi>x</mi><mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo>
///</math>
///
///which allows the implementation to reduce many arguments to a convenient range.
///
///For negative non-integer <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>x</mi></math>, the function is computed using Euler's reflection formula:
///
///<math display="block" xmlns="http://www.w3.org/1998/Math/MathML">
///  <mi>Γ</mi><mo>(</mo><mi>x</mi><mo>)</mo>
///  <mo>=</mo>
///  <mfrac>
///    <mi>π</mi>
///    <mrow>
///      <mi>sin</mi><mo>(</mo><mi>π</mi><mi>x</mi><mo>)</mo>
///      <mi>Γ</mi><mo>(</mo><mn>1</mn><mo>−</mo><mi>x</mi><mo>)</mo>
///    </mrow>
///  </mfrac>
///</math>
///
///Thus, the negative argument is converted to a positive argument.
///
///At zero and negative integers, Gamma has poles, so the function returns the corresponding infinite or non-finite result.
///
/// # Examples
/// ```
/// use abax::gamma;
/// let result = gamma(5.0);
/// assert_eq!(result, 24.0); // 4!
/// ```
pub fn gamma(x: f64) -> f64 {
    gammafn(x)
}

#[allow(dead_code)]
fn gammalims() -> [f64; 2] {
    let xmin = -170.5674972726612;
    let xmax = 171.61447887182298;
    [xmin, xmax]
}

pub(crate) fn gammafn(x: f64) -> f64 {
    const GAMCS: [f64; 42] = [
    	 0.8571195590989331421920062399942e-2,
    	 0.4415381324841006757191315771652e-2,
    	 0.5685043681599363378632664588789e-1,
    	-0.4219835396418560501012500186624e-2,
    	 0.1326808181212460220584006796352e-2,
    	-0.1893024529798880432523947023886e-3,
    	 0.3606925327441245256578082217225e-4,
    	-0.6056761904460864218485548290365e-5,
    	 0.1055829546302283344731823509093e-5,
    	-0.1811967365542384048291855891166e-6,
    	 0.3117724964715322277790254593169e-7,
    	-0.5354219639019687140874081024347e-8,
    	 0.9193275519859588946887786825940e-9,
    	-0.1577941280288339761767423273953e-9,
    	 0.2707980622934954543266540433089e-10,
    	-0.4646818653825730144081661058933e-11,
    	 0.7973350192007419656460767175359e-12,
    	-0.1368078209830916025799499172309e-12,
    	 0.2347319486563800657233471771688e-13,
    	-0.4027432614949066932766570534699e-14,
    	 0.6910051747372100912138336975257e-15,
    	-0.1185584500221992907052387126192e-15,
    	 0.2034148542496373955201026051932e-16,
    	-0.3490054341717405849274012949108e-17,
    	 0.5987993856485305567135051066026e-18,
    	-0.1027378057872228074490069778431e-18,
    	 0.1762702816060529824942759660748e-19,
    	-0.3024320653735306260958772112042e-20,
    	 0.5188914660218397839717833550506e-21,
    	-0.8902770842456576692449251601066e-22,
    	 0.1527474068493342602274596891306e-22,
    	-0.2620731256187362900257328332799e-23,
    	 0.4496464047830538670331046570666e-24,
    	-0.7714712731336877911703901525333e-25,
    	 0.1323635453126044036486572714666e-25,
    	-0.2270999412942928816702313813333e-26,
    	 0.3896418998003991449320816639999e-27,
    	-0.6685198115125953327792127999999e-28,
    	 0.1146998663140024384347613866666e-28,
    	-0.1967938586345134677295103999999e-29,
    	 0.3376448816585338090334890666666e-30,
    	-0.5793070335782135784625493333333e-31,
    ];

    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
    * (xmin, xmax) are non-trivial, see ./gammalims.c
    * xsml = exp(.01)*DBL_MIN
    * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
    */
    const NGAM: usize = 22;
    const XMIN: f64 = -170.5674972726612;
    const XMAX: f64 =  171.61447887182298;
    const XSML: f64 = 2.2474362225598545e-308;
    const DXREL: f64 = 1.490116119384765696e-8;

    if ISNAN(x) {
        return x;
    }

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if x == 0.0 || (x < 0.0 && x == f64::round(x)) {
        return f64::NAN;
    }

    let y = f64::abs(x);
    let mut value: f64;
    if y <= 10.0 {
        /* Compute gamma(x) for -10 <= x <= 10
        * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
        * first of all. */

    	let mut n: f64 = f64::floor(x);
	    let y = x - n as f64; /* y in [ 0, 1 ) */
	    n = n - 1.0;
	    value = chebyshev_eval(y * 2.0 - 1.0, &GAMCS, NGAM) + 0.9375;
	    if n == 0.0 {
            return value;
        }

	    if n < 0.0 {
	        /* compute gamma(x) for -10 <= x < 1 */

	        /* exact 0 or "-n" checked already above */

	        /* The answer is less than half precision */
	        /* because x too near a negative integer. */
	        if x < -0.5 && f64::abs((x - (x - 0.5).floor()) / x) < DXREL {
                // warning about precision issue of gammafn
	        }

	        /* The argument is so close to 0 that the result would overflow. */
	        if y < XSML {
                // warning about range issue of gammafn
		        if x > 0.0 {
                    return ML_POSINF;
                } else {
                    return ML_NEGINF;
                }
	        }

	        n = -n;

            let mut i = 0.0;
            while i < n {
                value /= x + i;
                i += 1.0;
            }

	        return value;
	    } else {
	        /* gamma(x) for 2 <= x <= 10 */
            let mut i = 1.0;
            while i <= n {
                value *= y + i;
                i += 1.0;
            }
	        return value;
	    }
    } else {
	    /* gamma(x) for	 y = |x| > 10. */

	    if x > XMAX {
            /* Overflow */
	        return ML_POSINF;
	    }

	    if x < XMIN {
            /* Underflow */
	        return 0.0;
	    }

	    if y <= 50.0 && y == f64::floor(y) { /* compute (n - 1)! */
	        value = 1.0;
            let mut i = 2.0;
            while i < y {
                value *= i;
                i += 1.0;
            }
	    } else { /* normal case */
	        value = f64::exp((y - 0.5) * f64::ln(y) - y + M_LN_SQRT_2PI +
            if 2.0 * y == f64::floor(2.0 * y) {stirlerr(y)} else {lgammacor(y)});
        }

        if x > 0.0 {
            return value;
        }
	    // else:  x < 0, not an integer :

        if f64::abs((x - f64::floor(x - 0.5))/x) < DXREL {
            /* The answer is less than half precision because */
            /* the argument is too near a negative integer. */
        }

        let sinpiy = sinpi(y);
        if sinpiy == 0.0 {
            /* Negative integer arg - overflow */
            return ML_POSINF;
        }

        return -M_PI / (y * sinpiy * value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-14;


    #[test]
    fn test_exact_integers() {
        // Gamma(n) = (n-1)!
        assert_eq!(gamma(1.0), 1.0);
        assert_eq!(gamma(2.0), 1.0);
        assert_eq!(gamma(3.0), 2.0);
        assert_eq!(gamma(4.0), 6.0);
        assert_eq!(gamma(5.0), 24.0);
        assert_eq!(gamma(10.0), 362880.0);
    }

    #[test]
    fn test_half_integers() {
        // Gamma(1/2) = sqrt(pi)
        let sqrt_pi = std::f64::consts::PI.sqrt();
        assert!((gamma(0.5) - sqrt_pi).abs() < EPSILON);

        // Gamma(3/2) = 1/2 * sqrt(pi)
        assert!((gamma(1.5) - 0.5 * sqrt_pi).abs() < EPSILON);

        // Gamma(5/2) = 3/4 * sqrt(pi)
        assert!((gamma(2.5) - 0.75 * sqrt_pi).abs() < EPSILON);
    }

    #[test]
    fn test_recurrence_relation() {
        // Gamma(x + 1) = x * Gamma(x)
        let x = std::f64::consts::PI;
        let lhs = gamma(x + 1.0);
        let rhs = x * gamma(x);
        assert!((lhs - rhs).abs() / lhs < EPSILON);
    }

    #[test]
    fn test_negative_values() {
        // Test a few known negative points using reflection
        // Gamma(-0.5) = -2 * sqrt(pi)
        let expected = -2.0 * std::f64::consts::PI.sqrt();
        assert!((gamma(-0.5) - expected).abs() < EPSILON);
    }

    #[test]
    fn test_special_cases() {
        // Poles (Returns Infinity or NaN based on your implementation)
        assert!(gamma(0.0).is_nan());
        assert!(gamma(-1.0).is_nan());

        // Limits
        assert!(gamma(f64::NAN).is_nan());
        assert!(gamma(f64::INFINITY).is_infinite());
    }
}
