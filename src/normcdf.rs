use crate::consts::{M_1_SQRT_2PI, M_SQRT_32};
use crate::erfc;
use crate::rutils::*;
use std::f64::consts::SQRT_2;

/// Normal cumulative distribution function (CDF).
///
/// Given a value `x`, a mean `mu`, and a standard deviation `sigma`,
/// this function returns the probability that a normal random variable
/// is less than or equal to `x`.
///
/// See also [`pnorm`](./fn.pnorm.html)
///
/// # Mathematical Definition
/// For a normal distribution with mean <math><mi>μ</mi></math> and standard deviation <math><mi>σ</mi></math>:
/// - Lower tail (`upper = false`): <math display="block"><mrow><mi>P</mi><mo stretchy="false">(</mo><mi>X</mi><mo>&le;</mo><mi>x</mi><mo stretchy="false">)</mo><mo>=</mo><mi mathvariant="normal">&Phi;</mi><mrow><mo>(</mo><mfrac><mrow><mi>x</mi><mo>&minus;</mo><mi>&mu;</mi></mrow><mi>&sigma;</mi></mfrac><mo>)</mo></mrow><mo>=</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mrow><mo>[</mo><mrow><mn>1</mn><mo>+</mo><mo>erf</mo><mrow><mo>(</mo><mfrac><mrow><mi>x</mi><mo>&minus;</mo><mi>&mu;</mi></mrow><mrow><mi>&sigma;</mi><msqrt><mn>2</mn></msqrt></mrow></mfrac><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></math>
/// - Upper tail (`upper = true`): <math display="block"><mrow><mi>P</mi><mo stretchy="false">(</mo><mi>X</mi><mo>&gt;</mo><mi>x</mi><mo stretchy="false">)</mo><mo>=</mo><mn>1</mn><mo>&minus;</mo><mi>P</mi><mo stretchy="false">(</mo><mi>X</mi><mo>&le;</mo><mi>x</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mo>erfc</mo><mrow><mo>(</mo><mfrac><mrow><mi>x</mi><mo>&minus;</mo><mi>&mu;</mi></mrow><mrow><mi>&sigma;</mi><msqrt><mn>2</mn></msqrt></mrow></mfrac><mo>)</mo></mrow></mrow></math>
///
/// # Examples
/// ```
/// use abax::{normcdf, pnorm};
///
/// // Standard normal median is 0.5
/// assert!((normcdf(0.0, 0.0, 1.0, false) - 0.5).abs() < 1e-15);
/// assert!((pnorm(0.0, 0.0, 1.0, true, false) - 0.5).abs() < 1e-15);
/// // One sigma upper bound
/// assert!((normcdf(1.0, 0.0, 1.0, false) - 0.8413447460685429).abs() < 1e-15);
/// assert!((pnorm(1.0, 0.0, 1.0, true, false) - 0.8413447460685429).abs() < 1e-15);
/// // Upper tail one sigma
/// assert!((normcdf(1.0, 0.0, 1.0, true) - 0.15865525393145705).abs() < 1e-15);
/// assert!((pnorm(1.0, 0.0, 1.0, false, false) - 0.15865525393145705).abs() < 1e-15);
/// ```
pub fn normcdf(x: f64, mu: f64, sigma: f64, upper: bool) -> f64 {
    if x.is_nan() || mu.is_nan() || sigma.is_nan() || sigma < 0.0 {
        return f64::NAN;
    }

    if sigma == 0.0 {
        return if upper {
            if x < mu { 1.0 } else { 0.0 }
        } else {
            if x < mu { 0.0 } else { 1.0 }
        };
    }

    let z = (x - mu) / sigma;
    if upper {
        0.5 * erfc(z / SQRT_2)
    } else {
        0.5 * erfc(-z / SQRT_2)
    }
}

/// See [`normcdf`](./fn.normcdf.html) for details.
pub fn pnorm(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    pnorm5(x, mu, sigma, if lower_tail { 1 } else { 0 }, log_p)
}

fn pnorm5(x: f64, mu: f64, sigma: f64, lower_tail: i32, log_p: bool) -> f64 {
    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return f64::NAN;
    }
    if !x.is_finite() && mu == x {
        return f64::NAN; /* x-mu is NaN */
    }
    if sigma <= 0.0 {
        if sigma < 0.0 {
            return f64::NAN;
        }
        /* sigma = 0 : */
        return if x < mu {
            R_DT_0(if lower_tail != 0 { true } else { false }, log_p)
        } else {
            R_DT_1(if lower_tail != 0 { true } else { false }, log_p)
        };
    }
    let mut p: f64 = (x - mu) / sigma;
    if !R_FINITE(p) {
        return if x < mu {
            R_DT_0(if lower_tail != 0 { true } else { false }, log_p)
        } else {
            R_DT_1(if lower_tail != 0 { true } else { false }, log_p)
        };
    }
    let x = p;

    let mut cp: f64 = 0.0;
    pnorm_both(
        x,
        &mut p,
        &mut cp,
        if lower_tail != 0 { 0 } else { 1 },
        log_p,
    );

    return if lower_tail != 0 { p } else { cp };
}

/// i_tail: 0 = lower, 1 = upper, 2 = both.
/// if(lower) return *cum := P[X ≤ x]
/// if(upper) return **ccum := P[X > x] = 1 - P[x ≤ x]
fn pnorm_both(x: f64, cum: &mut f64, ccum: &mut f64, i_tail: i32, log_p: bool) -> () {
    const A: [f64; 5] = [
        2.2352520354606839287e+0,
        1.6102823106855587881e+2,
        1.0676894854603709582e+3,
        1.8154981253343561249e+4,
        6.5682337918207449113e-2,
    ];
    const B: [f64; 4] = [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956,
    ];
    const C: [f64; 9] = [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8,
    ];
    const D: [f64; 8] = [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727,
    ];
    const P: [f64; 6] = [
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303,
    ];
    const Q: [f64; 5] = [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5,
    ];

    let d_2 = |x: f64| x / 2.0;
    let do_del = |x: f64,
                  lower: bool,
                  upper: bool,
                  log_p: bool,
                  temp: f64,
                  cum: &mut f64,
                  ccum: &mut f64| {
        let xsq = (x * 16.0).trunc() / 16.0;
        let del = (x - xsq) * (x + xsq);
        if log_p {
            *cum = (-xsq * d_2(xsq)) - d_2(del) + f64::ln(temp);
            if (lower && x > 0.0) || (upper && x <= 0.0) {
                *ccum = f64::ln_1p(-f64::exp(-xsq * d_2(xsq)) * f64::exp(-d_2(del)) * temp);
            }
        } else {
            *cum = f64::exp(-xsq * d_2(xsq)) * f64::exp(-d_2(del)) * temp;
            *ccum = 1.0 - *cum;
        }
    };
    let swap_tail = |x: f64, lower: bool, cum: &mut f64, ccum: &mut f64| {
        if x > 0.0 {
            /* swap ccum <--> cum */
            let temp = *cum;
            if lower {
                *cum = *ccum;
            }
            *ccum = temp;
        }
    };

    if x.is_nan() {
        *cum = f64::NAN;
        *ccum = f64::NAN;
        return;
    }

    /* Consider changing these : */
    const EPS: f64 = f64::EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    let lower = i_tail != 1;
    let upper = i_tail != 0;

    let y = f64::abs(x);
    if y <= 0.67448975 {
        let mut xnum: f64;
        let mut xden: f64;
        if y > EPS {
            let xsq = x * x;
            xnum = A[4] * xsq;
            xden = xsq;
            for i in 0..3 {
                xnum = (xnum + A[i]) * xsq;
                xden = (xden + B[i]) * xsq;
            }
        } else {
            xnum = 0.0;
            xden = 0.0;
        }

        let temp = x * (xnum + A[3]) / (xden + B[3]);
        if lower {
            *cum = 0.5 + temp;
        }
        if upper {
            *ccum = 0.5 - temp;
        }
        if log_p {
            if lower {
                *cum = f64::ln(*cum);
            }
            if upper {
                *ccum = f64::ln(*ccum);
            }
        }
    } else if y <= M_SQRT_32 {
        /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
        let mut xnum: f64 = C[8] * y;
        let mut xden: f64 = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);

        do_del(y, lower, upper, log_p, temp, cum, ccum);
        swap_tail(x, lower, cum, ccum);
    } else if (log_p && y < 1e170)
        || (lower && -38.4674 < x && x < 8.2924)
        || (upper && -8.2924 < x && x < 38.4674)
    {
        /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
        let xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }
        let temp = xsq * (xnum + P[4]) / (xden + Q[4]);
        let temp = (M_1_SQRT_2PI - temp) / y;

        do_del(x, lower, upper, log_p, temp, cum, ccum);
        swap_tail(x, lower, cum, ccum);
    } else {
        /* large |x| such that probs are 0 or 1 */
        if x > 0.0 {
            *cum = R_D__1(log_p);
            *ccum = R_D__0(log_p);
        } else {
            *cum = R_D__0(log_p);
            *ccum = R_D__1(log_p);
        }
    }

    return;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normcdf_standard() {
        let tol = 1e-14;
        assert!((normcdf(0.0, 0.0, 1.0, false) - 0.5).abs() < tol);
        assert!((pnorm(0.0, 0.0, 1.0, true, false) - 0.5).abs() < tol);
        assert!((normcdf(1.959963984540054, 0.0, 1.0, false) - 0.975).abs() < tol);
        assert!((pnorm(1.959963984540054, 0.0, 1.0, true, false) - 0.975).abs() < tol);
        assert!((normcdf(-1.959963984540054, 0.0, 1.0, false) - 0.025).abs() < tol);
        assert!((pnorm(-1.959963984540054, 0.0, 1.0, true, false) - 0.025).abs() < tol);
    }

    #[test]
    fn test_normcdf_upper() {
        let tol = 1e-14;
        assert!((normcdf(0.0, 0.0, 1.0, true) - 0.5).abs() < tol);
        assert!((pnorm(0.0, 0.0, 1.0, false, false) - 0.5).abs() < tol);
        assert!((normcdf(1.959963984540054, 0.0, 1.0, true) - 0.025).abs() < tol);
        assert!((pnorm(1.959963984540054, 0.0, 1.0, false, false) - 0.025).abs() < tol);
        assert!((normcdf(-1.959963984540054, 0.0, 1.0, true) - 0.975).abs() < tol);
        assert!((pnorm(-1.959963984540054, 0.0, 1.0, false, false) - 0.975).abs() < tol);
    }

    #[test]
    fn test_normcdf_zero_sigma() {
        assert_eq!(normcdf(5.0, 5.0, 0.0, false), 1.0);
        assert_eq!(pnorm(5.0, 5.0, 0.0, true, false), 1.0);
        assert_eq!(normcdf(4.99, 5.0, 0.0, false), 0.0);
        assert_eq!(pnorm(4.99, 5.0, 0.0, true, false), 0.0);
        assert_eq!(normcdf(4.99, 5.0, 0.0, true), 1.0);
        assert_eq!(pnorm(4.99, 5.0, 0.0, false, false), 1.0);
    }

    #[test]
    fn test_normcdf_invalid() {
        assert!(normcdf(0.5, 0.0, -1.0, false).is_nan());
        assert!(pnorm(0.5, 0.0, -1.0, true, false).is_nan());
    }
}
