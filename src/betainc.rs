use crate::toms708::bratio;
use crate::rutils::{ISNAN, R_DT_1, R_DT_0, R_FINITE};
use crate::consts::M_LN2;

/// Regularized incomplete beta function <math><msub><mi>I</mi><mi>x</mi></msub><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo></math>.
///
/// See also [`pbeta`](./fn.pbeta.html).
/// 
/// Solves for:
/// - <math><msub><mi>I</mi><mi>x</mi></msub><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo></math> when `lower = true` (regularized lower incomplete beta)
/// - <math><mn>1</mn><mo>-</mo><msub><mi>I</mi><mi>x</mi></msub><mo>(</mo><mi>z</mi><mo>,</mo><mi>w</mi><mo>)</mo></math> when `lower = false` (regularized upper incomplete beta)
/// 
/// The regularized incomplete beta function:
/// 
/// <math display="block">
///   <msub>
///     <mi>I</mi>
///     <mi>x</mi>
///   </msub>
///   <mo stretchy="false">(</mo>
///   <mi>z</mi>
///   <mo>,</mo>
///   <mi>w</mi>
///   <mo stretchy="false">)</mo>
///   <mo>=</mo>
///   <mfrac>
///     <mrow>
///       <mi>B</mi>
///       <mo stretchy="false">(</mo>
///       <mi>x</mi>
///       <mo>;</mo>
///       <mi>z</mi>
///       <mo>,</mo>
///       <mi>w</mi>
///       <mo stretchy="false">)</mo>
///     </mrow>
///     <mrow>
///       <mi>B</mi>
///       <mo stretchy="false">(</mo>
///       <mi>z</mi>
///       <mo>,</mo>
///       <mi>w</mi>
///       <mo stretchy="false">)</mo>
///     </mrow>
///   </mfrac>
///   <mo>=</mo>
///   <mfrac>
///     <mn>1</mn>
///     <mrow>
///       <mi>B</mi>
///       <mo stretchy="false">(</mo>
///       <mi>z</mi>
///       <mo>,</mo>
///       <mi>w</mi>
///       <mo stretchy="false">)</mo>
///     </mrow>
///   </mfrac>
///   <msubsup>
///     <mo>&int;</mo>
///     <mn>0</mn>
///     <mi>x</mi>
///   </msubsup>
///   <msup>
///     <mi>t</mi>
///     <mrow>
///       <mi>z</mi>
///       <mo>&minus;</mo>
///       <mn>1</mn>
///     </mrow>
///   </msup>
///   <msup>
///     <mrow>
///       <mo stretchy="false">(</mo>
///       <mn>1</mn>
///       <mo>&minus;</mo>
///       <mi>t</mi>
///       <mo stretchy="false">)</mo>
///     </mrow>
///     <mrow>
///       <mi>w</mi>
///       <mo>&minus;</mo>
///       <mn>1</mn>
///     </mrow>
///   </msup>
///   <mspace width="0.167em" />
///   <mi>d</mi>
///   <mi>t</mi>
/// </math>
/// 
/// The complete beta function:
/// 
/// <math display="block">
///  <mi>B</mi>
///  <mo stretchy="false">(</mo>
///  <mi>z</mi>
///  <mo>,</mo>
///  <mi>w</mi>
///  <mo stretchy="false">)</mo>
///  <mo>=</mo>
///  <msubsup>
///    <mo>&#x222B;</mo>
///    <mn>0</mn>
///    <mn>1</mn>
///  </msubsup>
///  <msup>
///    <mi>t</mi>
///    <mrow>
///      <mi>z</mi>
///      <mo>&#x2212;</mo>
///      <mn>1</mn>
///    </mrow>
///  </msup>
///  <msup>
///    <mrow>
///      <mo stretchy="false">(</mo>
///      <mn>1</mn>
///      <mo>&#x2212;</mo>
///      <mi>t</mi>
///      <mo stretchy="false">)</mo>
///    </mrow>
///    <mrow>
///      <mi>w</mi>
///      <mo>&#x2212;</mo>
///      <mn>1</mn>
///    </mrow>
///  </msup>
///  <mspace width="0.167em" />
///  <mi>d</mi>
///  <mi>t</mi>
///  <mo>=</mo>
///  <mfrac>
///    <mrow>
///      <mi mathvariant="normal">&#x0393;</mi>
///      <mo stretchy="false">(</mo>
///      <mi>z</mi>
///      <mo stretchy="false">)</mo>
///      <mi mathvariant="normal">&#x0393;</mi>
///      <mo stretchy="false">(</mo>
///      <mi>w</mi>
///      <mo stretchy="false">)</mo>
///    </mrow>
///    <mrow>
///      <mi mathvariant="normal">&#x0393;</mi>
///      <mo stretchy="false">(</mo>
///      <mi>z</mi>
///      <mo>+</mo>
///      <mi>w</mi>
///      <mo stretchy="false">)</mo>
///    </mrow>
///  </mfrac>
///</math>

pub fn betainc(x: f64, z: f64, w: f64, lower: bool) -> f64 {
    pbeta(x, z, w, lower, false)
}

fn pbeta_raw(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    if x >= 1.0 { // may happen when called from qbeta()
        return R_DT_1(lower_tail, log_p);
    }
    // treat limit cases correctly here:
    if a == 0.0 || b == 0.0 || !R_FINITE(a) || !R_FINITE(b) {
	    // NB:  0 <= x < 1 :
	    if a == 0.0 && b == 0.0 { // point mass 1/2 at each of {0,1} :
	        return if log_p { -M_LN2 } else { 0.5 };
        }
	    if a == 0.0 || a/b == 0.0 { // point mass 1 at 0 ==> P(X <= x) = 1, all x >= 0
	        return R_DT_1(lower_tail, log_p);
        }
	    if b == 0.0 || b/a == 0.0 { // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
	        return R_DT_0(lower_tail, log_p);
        }
	    // else, remaining case:  a = b = Inf : point mass 1 at 1/2
	    return if x < 0.5 { R_DT_0(lower_tail, log_p) } else { R_DT_1(lower_tail, log_p) };
    }
    if x <= 0.0 {
        return R_DT_0(lower_tail, log_p);
    }

    // Now:  0 < a < Inf;  0 < b < Inf  and  0 < x < 1
    let x1 = 0.5 - x + 0.5;
    let ans = bratio(a, b, x, x1);
    return match log_p {
        true => if lower_tail { f64::ln(ans.w) } else { f64::ln(ans.w1) },
        false => if lower_tail { ans.w } else { ans.w1 },
    };
}

/// See [`betainc`](./fn.betainc.html) for details.
pub fn pbeta(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    if ISNAN(x) || ISNAN(a) || ISNAN(b) {
        return f64::NAN;
    }

    if a < 0.0 || b < 0.0 {
        return f64::NAN;
    }

    // allowing a==0 and b==0  <==> treat as one- or two-point mass
    return pbeta_raw(x, a, b, lower_tail, log_p);
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_betainc_basic() {
        // I_0.5(1, 1) = 0.5
        assert!((betainc(0.5, 1.0, 1.0, true) - 0.5).abs() < 1e-15);
        // I_0.5(2, 2) = 0.5
        assert!((betainc(0.5, 2.0, 2.0, true) - 0.5).abs() < 1e-15);
        // I_0.2(1, 3) = 1 - (1-0.2)^3 = 0.488
        assert!((betainc(0.2, 1.0, 3.0, true) - 0.488).abs() < 1e-15);
    }

    #[test]
    fn test_betainc_symmetry() {
        let x = 0.3;
        let z = 2.5;
        let w = 1.5;
        let lower = betainc(x, z, w, true);
        let upper = betainc(x, z, w, false);
        assert!((lower + upper - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_betainc_boundaries() {
        assert_eq!(betainc(0.0, 1.0, 1.0, true), 0.0);
        assert_eq!(betainc(1.0, 1.0, 1.0, true), 1.0);
    }

    #[test]
    fn test_betainc_rejects_infinite_shape_parameters() {
        assert!(betainc(0.5, f64::INFINITY, 1.0, true) == 0.0);
        assert!(betainc(0.5, 1.0, f64::INFINITY, true) == 1.0);
        assert!(betainc(0.5, f64::NEG_INFINITY, 1.0, true).is_nan());
        assert!(betainc(0.5, 1.0, f64::NEG_INFINITY, true).is_nan());
    }
}
