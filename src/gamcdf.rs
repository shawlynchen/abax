use crate::gammainc;
use crate::gammaln::lgammafn;

/// Gamma cumulative distribution function (CDF).
///
/// Returns the probability that a Gamma random variable with shape parameter `a`
/// and scale parameter `b` is less than or equal to `x`.
///
/// # Mathematical Definition
/// For a Gamma distribution with shape <math><mi>a</mi></math> and scale <math><mi>b</mi></math>:
/// - Lower tail (`upper = false`):
/// <math display="block">
///   <mi>F</mi><mo>(</mo><mi>x</mi><mo>;</mo><mi>a</mi><mo>,</mo><mi>b</mi><mo>)</mo>
///   <mo>=</mo>
///   <mi>P</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>/</mo><mi>b</mi><mo>)</mo>
///   <mo>=</mo>
///   <mfrac><mn>1</mn><mrow><mi>Γ</mi><mo>(</mo><mi>a</mi><mo>)</mo></mrow></mfrac>
///   <msubsup><mo>∫</mo><mn>0</mn><mrow><mi>x</mi><mo>/</mo><mi>b</mi></mrow></msubsup>
///   <msup><mi>t</mi><mrow><mi>a</mi><mo>-</mo><mn>1</mn></mrow></msup><msup><mi>e</mi><mrow><mo>-</mo><mi>t</mi></mrow></msup><mi>dt</mi>
/// </math>
/// - Upper tail (`upper = true`):
/// <math display="block">
///   <mn>1</mn><mo>-</mo><mi>F</mi><mo>(</mo><mi>x</mi><mo>;</mo><mi>a</mi><mo>,</mo><mi>b</mi><mo>)</mo>
///   <mo>=</mo>
///   <mi>Q</mi><mo>(</mo><mi>a</mi><mo>,</mo><mi>x</mi><mo>/</mo><mi>b</mi><mo>)</mo>
/// </math>
///
/// # Domain
/// - <math><mi>x</mi><mo>≥</mo><mn>0</mn></math> (Returns 0 for the lower tail if <math><mi>x</mi><mo>&lt;</mo><mn>0</mn></math>).
/// - <math><mi>a</mi><mo>&gt;</mo><mn>0</mn></math> (Shape parameter).
/// - <math><mi>b</mi><mo>&gt;</mo><mn>0</mn></math> (Scale parameter).
/// - Returns `NaN` if `a <= 0` or `b <= 0`.
///
/// # Examples
/// ```
/// use abax::gamcdf;
///
/// // For a=1, Gamma reduces to Exponential distribution: 1 - exp(-x/b)
/// let p = gamcdf(1.0, 1.0, 1.0, false);
/// assert!((p - (1.0 - (-1.0f64).exp())).abs() < 1e-15);
/// ```
pub fn gamcdf(x: f64, a: f64, b: f64, upper: bool) -> f64 {
    let a = if a < 0.0 { f64::NAN } else { a };
    let b = if b <= 0.0 { f64::NAN } else { b };
    let x = if x < 0.0 { 0.0 } else { x };

    let z = x / b;
    gammainc(z, a, !upper, false)
}

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxiliary in log1pmx() and lgamma1p()
 */
fn logcf(x: f64, i: f64, d: f64, eps: f64) -> f64 {
    let mut c1: f64 = 2.0 * d;
    let mut c2: f64 = i + d;
    let mut c4: f64 = c2 + d;
    let mut a1: f64 = c2;
    let mut b1: f64 = i * (c2 - i * x);
    let mut b2: f64 = d * d * x;
    let mut a2: f64 = c4 * c2 - b2;

    b2 = c4 * b1 - i * b2;

    const SQR: fn(f64) -> f64 = |x| x * x;
    let scalefactor: f64 = SQR(SQR(SQR(4294967296.0)));

    while f64::abs(a2 * b1 - a1 * b2) > f64::abs(eps * b1 * b2) {
	    let mut c3 = c2 * c2 * x;
	    c2 += d;
	    c4 += d;
	    a1 = c4 * a2 - c3 * a1;
	    b1 = c4 * b2 - c3 * b1;

	    c3 = c1 * c1 * x;
	    c1 += d;
	    c4 += d;
	    a2 = c4 * a1 - c3 * a2;
	    b2 = c4 * b1 - c3 * b2;

	    if f64::abs(b2) > scalefactor {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        } else if f64::abs (b2) < 1.0 / scalefactor {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
fn log1pmx(x: f64) -> f64 {
    const MINLOG1VALUE: f64 = -0.79149064;

    if x > 1.0 || x < MINLOG1VALUE {
	    return f64::ln_1p(x) - x;
    } else {
	    let r = x / (2.0 + x);
        let y = r * r;
	    if f64::abs(x) < 1e-2 {
	        const TWO: f64 = 2.0;
	        return r * ((((TWO / 9.0 * y + TWO / 7.0) * y + TWO / 5.0) * y + TWO / 3.0) * y - x);
	    } else {
	        const TOL_LOGCF: f64 = 1e-14;
	        return r * (2.0 * y * logcf (y, 3.0, 2.0, TOL_LOGCF) - x);
	    }
    }
}

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
pub(crate) fn lgamma1p (a: f64) -> f64 {
    if f64::abs (a) >= 0.5 {
	    return lgammafn (a + 1.0);
    }

    const EULERS_CONST: f64 = 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const N: i32 = 40;
    const COEFFS: [f64; 40] = [
        0.3224670334241132182362075833230126e-0,
        0.6735230105319809513324605383715000e-1,
        0.2058080842778454787900092413529198e-1,
        0.7385551028673985266273097291406834e-2,
        0.2890510330741523285752988298486755e-2,
        0.1192753911703260977113935692828109e-2,
        0.5096695247430424223356548135815582e-3,
        0.2231547584535793797614188036013401e-3,
        0.9945751278180853371459589003190170e-4,
        0.4492623673813314170020750240635786e-4,
        0.2050721277567069155316650397830591e-4,
        0.9439488275268395903987425104415055e-5,
        0.4374866789907487804181793223952411e-5,
        0.2039215753801366236781900709670839e-5,
        0.9551412130407419832857179772951265e-6,
        0.4492469198764566043294290331193655e-6,
        0.2120718480555466586923135901077628e-6,
        0.1004322482396809960872083050053344e-6,
        0.4769810169363980565760193417246730e-7,
        0.2271109460894316491031998116062124e-7,
        0.1083865921489695409107491757968159e-7,
        0.5183475041970046655121248647057669e-8,
        0.2483674543802478317185008663991718e-8,
        0.1192140140586091207442548202774640e-8,
        0.5731367241678862013330194857961011e-9,
        0.2759522885124233145178149692816341e-9,
        0.1330476437424448948149715720858008e-9,
        0.6422964563838100022082448087644648e-10,
        0.3104424774732227276239215783404066e-10,
        0.1502138408075414217093301048780668e-10,
        0.7275974480239079662504549924814047e-11,
        0.3527742476575915083615072228655483e-11,
        0.1711991790559617908601084114443031e-11,
        0.8315385841420284819798357793954418e-12,
        0.4042200525289440065536008957032895e-12,
        0.1966475631096616490411045679010286e-12,
        0.9573630387838555763782200936508615e-13,
        0.4664076026428374224576492565974577e-13,
        0.2273736960065972320633279596737272e-13,
        0.1109139947083452201658320007192334e-13,
    ];

    const C: f64 = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const TOL_LOGCF: f64 = 1e-14;

    let mut lgam: f64 = C * logcf(-a / 2.0, N as f64 + 2.0, 1.0, TOL_LOGCF);
    for i in (0..N).rev() {
	    lgam = COEFFS[i as usize] - a * lgam;
    }

    return (a * lgam - EULERS_CONST) * a - log1pmx(a);
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamcdf_exponential_identity() {
        // If a=1, Gamma(1, b) is an Exponential distribution with mean b.
        // CDF is 1 - exp(-x/b).
        let a = 1.0;
        let b = 2.5;
        let x = 1.0;
        let p = gamcdf(x, a, b, false);
        let expected = 1.0 - (-x / b).exp();
        assert!((p - expected).abs() < 1e-15);
    }

    #[test]
    fn test_gamcdf_tail_consistency() {
        let x = 2.0;
        let a = 3.0;
        let b = 1.0;
        let p = gamcdf(x, a, b, false);
        let q = gamcdf(x, a, b, true);
        assert!((p + q - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_gamcdf_invalid_params() {
        assert!(gamcdf(1.0, 0.0, 1.0, false).is_nan());
        assert!(gamcdf(1.0, -1.0, 1.0, false).is_nan());
        assert!(gamcdf(1.0, 1.0, 0.0, false).is_nan());
        assert!(gamcdf(1.0, 1.0, -1.0, false).is_nan());
    }
}
