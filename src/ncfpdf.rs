use crate::fpdf::fpdf;
use crate::stirlerr::stirlerr;
use crate::binodeviance::binodeviance;
use std::f64::consts::PI;

const LNSR2PI: f64 = 0.918_938_533_204_672_8; // 0.5*log(2*pi)

/// Noncentral F probability density function.
///
/// `nu1`: numerator degrees of freedom
/// `nu2`: denominator degrees of freedom
/// `delta`: noncentrality parameter
pub fn ncfpdf(x: f64, nu1: f64, nu2: f64, delta: f64) -> f64 {
    // NaN handling
    if x.is_nan() || nu1.is_nan() || nu2.is_nan() || delta.is_nan() {
        return f64::NAN;
    }

    // Invalid parameters -> NaN
    let k = nu1 <= 0.0 || nu2 <= 0.0 || delta < 0.0;
    if k {
        return f64::NAN;
    }

    let mut y = 0.0;

    // Edge case x == 0
    if x == 0.0 {
        if nu1 < 2.0 {
            return f64::INFINITY;
        } else if nu1 == 2.0 {
            return (-delta / 2.0).exp();
        }
        // nu1 > 2: density is 0
        return 0.0;
    }

    // Central distribution
    if delta == 0.0 {
        return fpdf(x, nu1, nu2);
    }

    // Pre-divide by 2 for convenience
    let nu1 = nu1 / 2.0;
    let nu2 = nu2 / 2.0;
    let delta = delta / 2.0;

    // z, z1, xs
    let z = nu1 * x / (nu1 * x + nu2);
    let z1 = nu2 / (nu1 * x + nu2);
    let xs = delta * z;

    let intnu2 = nu2 == nu2.round();

    let mut big_k: f64 = 0.0;
    let mut termk: f64 = 0.0;

    if intnu2 {
        let smallx = xs <= nu1 / nu2;
        let largex = xs >= nu2 * (nu1 + nu2 - 1.0) && !smallx;
        if largex {
            big_k = nu2;
        } else if !smallx {
            let d = 0.5 * (1.0 - xs - nu1);
            big_k = (d + (d * d + xs * (nu2 + 1.0)).sqrt()).floor();
        }
        // K==0 case
        let keqnu2 = big_k == nu2;
        let zerok = (nu1 + big_k) <= 1.0;

        if keqnu2 {
            if z < 0.9 {
                termk = (nu1 + nu2 - 1.0) * z.ln();
            } else {
                termk = (nu1 + nu2 - 1.0) * (-z1).ln_1p();
            }
        } else if zerok {
            termk = stirlerr(nu1 + nu2)
                - stirlerr(nu1)
                - stirlerr(nu2)
                - binodeviance(nu1, (nu1 + nu2) * z)
                - binodeviance(nu2, (nu1 + nu2) * z1);
        } else {
            termk = stirlerr(nu1 + nu2 - 1.0)
                - stirlerr(nu1 + big_k - 1.0)
                - stirlerr(nu2 - big_k)
                - binodeviance(nu1 + big_k - 1.0, (nu1 + nu2 - 1.0) * z)
                - binodeviance(nu2 - big_k, (nu1 + nu2 - 1.0) * z1);
        }

        // Poisson density for leading term
        let x1 = delta * z1;
        let mut rsum = 0.0_f64;
        if big_k <= x1 * f64::MIN_POSITIVE {
            y = termk - x1;
        } else {
            y = termk
                - LNSR2PI
                - 0.5 * big_k.ln()
                - stirlerr(big_k)
                - binodeviance(big_k, x1);
        }

        // Sum recursively downwards
        let mut term = 1.0;
        let mut kk = big_k;
        while kk > 0.0 {
            kk -= 1.0;
            term = term * (kk + 1.0) * (kk + nu1) / (nu2 - kk) / xs;
            if term < f64::EPSILON * rsum.max(1.0) {
                break;
            }
            rsum += term;
        }

        // Sum recursively upwards
        term = 1.0;
        let mut kk = big_k;
        while kk < nu2 {
            term = term * xs * (nu2 - kk) / (kk + nu1) / (kk + 1.0);
            if term < f64::EPSILON * rsum.max(1.0) {
                break;
            }
            rsum += term;
            kk += 1.0;
        }

        // Combine
        let zerok = (nu1 + big_k) <= 1.0;
        let keqnu2 = !zerok && intnu2 && big_k == nu2;
        if zerok {
            y = y.exp()
                * (1.0 + rsum)
                * (nu1 * nu2 / (nu1 + nu2) / (2.0 * PI)).sqrt()
                / x;
        } else if keqnu2 {
            y = y.exp() * (1.0 + rsum) * nu1 * z1;
        } else {
            y = y.exp()
                * (1.0 + rsum)
                * nu1
                * z1
                * ((nu1 + nu2 - 1.0) / (nu2 - big_k) / (nu1 + big_k - 1.0) / (2.0 * PI))
                    .sqrt();
        }
    } else {
        // Non-integer nu2
        let largex = xs > nu1 / (nu1 + nu2);
        if largex {
            let d = 0.5 * (1.0 + xs - nu1);
            big_k = (d + (d * d + xs * (nu1 + nu2 - 1.0)).sqrt()).floor();
        }

        let zerok = (nu1 + big_k) <= 1.0;

        if zerok {
            termk = stirlerr(nu1 + nu2)
                - stirlerr(nu1)
                - stirlerr(nu2)
                - binodeviance(nu1, (nu1 + nu2) * z)
                - binodeviance(nu2, (nu1 + nu2) * z1);
        } else {
            termk = stirlerr(nu1 + nu2 + big_k - 1.0)
                - stirlerr(nu1 + big_k - 1.0)
                - stirlerr(nu2)
                - binodeviance(
                    nu1 + big_k - 1.0,
                    (nu1 + nu2 + big_k - 1.0) * z,
                )
                - binodeviance(nu2, (nu1 + nu2 + big_k - 1.0) * z1);
        }

        // Poisson density at leading term
        if big_k <= delta * f64::MIN_POSITIVE {
            y = termk - delta;
        } else {
            y = termk
                - LNSR2PI
                - 0.5 * big_k.ln()
                - stirlerr(big_k)
                - binodeviance(big_k, delta);
        }

        // Sum recursively downwards
        let mut rsum = 0.0_f64;
        let mut term = 1.0_f64;
        let mut kk = big_k;
        while kk > 0.0 {
            kk -= 1.0;
            term = term * (kk + 1.0) * (kk + nu1) / (kk + nu1 + nu2) / xs;
            if term < f64::EPSILON * rsum.max(1.0) {
                break;
            }
            rsum += term;
        }

        // Sum recursively upwards
        term = 1.0;
        let mut kk = big_k;
        loop {
            term = term * xs * (kk + nu1 + nu2) / (kk + nu1) / (kk + 1.0);
            if term < f64::EPSILON * rsum.max(1.0) {
                break;
            }
            rsum += term;
            kk += 1.0;
        }

        let zerok = (nu1 + big_k) <= 1.0;
        if zerok {
            y = y.exp()
                * (1.0 + rsum)
                * (nu1 * nu2 / (nu1 + nu2) / (2.0 * PI)).sqrt()
                / x;
        } else {
            y = y.exp()
                * (1.0 + rsum)
                * nu1
                * z1
                * ((nu1 + nu2 + big_k - 1.0) / nu2 / (nu1 + big_k - 1.0) / (2.0 * PI))
                    .sqrt();
        }
    }

    y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ncfpdf() {
        let x = 2.0;
        let nu1 = 5.0;
        let nu2 = 10.0;
        let delta = 3.0;
        let result = ncfpdf(x, nu1, nu2, delta); // 0.2492817442645997
        println!("ncfpdf({}, {}, {}, {}) = {}", x, nu1, nu2, delta, result); 
    }
}