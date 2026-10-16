use crate::amos::{d1mach, zabs, zdiv, zexp, zmlt, zsqrt};

/// Computes the I Bessel function for `Re(z) >= 0` by means of the
/// asymptotic expansion for large `|z|` in the region
/// `|z| > max(rl, fnu*fnu/2)`.
///
/// Rust port of the SLATEC `ZASYI` subroutine.
///
/// # Arguments
/// * `zr`, `zi` - real and imaginary parts of the argument, `Re(z) >= 0`
/// * `fnu`      - real order
/// * `kode`     - scaling flag: `1` unscaled, `2` scaled by `exp(-z)`
/// * `n`        - number of components
/// * `rl`       - order threshold (see prologue)
/// * `tol`      - tolerance
/// * `elim`     - `-log` of the largest representable number
/// * `alim`     - `-log` of the smallest representable number
///
/// # Returns
/// `(yr, yi, nz)` where `yr`/`yi` are the real/imaginary parts of the `n`
/// I Bessel values, and `nz` is the error flag:
///   * `0`  - normal return
///   * `-1` - overflow on `kode = 1` (`|Re(z)| > elim`)
///   * `-2` - series did not converge within `jl` terms
///
/// # Panics
/// Panics if `n < 1`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zasyi(zr: f64, zi: f64, fnu: f64, kode: i32, n: i32, yr: &mut [f64], yi: &mut [f64], nz: &mut i32, rl: f64, tol: f64, elim: f64, alim: f64) {
    debug_assert!(n >= 1, "ZASYI: N must be >= 1");
    debug_assert!(n as usize <= yr.len());
    debug_assert!(n as usize <= yi.len());
    debug_assert_eq!(yr.len(), yi.len());

    const PI: f64 = 3.14159265358979324;
    const RTPI: f64 = 0.159154943091895336;

    let zeror = 0.0_f64;
    let zeroi = 0.0_f64;
    let coner = 1.0_f64;
    let conei = 0.0_f64;

    *nz = 0;

    let az = zabs(zr, zi);

    let arm  = 1.0e3 * d1mach(1);
    let rtr1 = arm.sqrt();

    let il   = 2.min(n);
    let dfnu = fnu + (n - il) as f64;

    // OVERFLOW TEST
    let raz = 1.0 / az;
    let (mut str, mut sti) = (zr * raz, -zi * raz);
    let (mut ak1r, mut ak1i) = (RTPI * str * raz, RTPI * sti * raz);
    zsqrt(ak1r, ak1i, &mut ak1r, &mut ak1i);

    let (czr, czi) = if kode == 2 { (zeror, zi) } else { (zr, zi) };

    if czr.abs() > elim {
        *nz = -1;
        return;
    }

    let dnu2 = dfnu + dfnu;
    let mut koded: i32 = 1;

    match czr.abs() > alim && n > 2 {
        true => {},
        false => {
            koded = 0;
            zexp(czr, czi, &mut str, &mut sti);
            zmlt(ak1r, ak1i, str, sti, &mut ak1r, &mut ak1i);
        },
    }

    let mut fdn = 0.0_f64;
    if dnu2 > rtr1 {
        fdn = dnu2 * dnu2;
    }
    let (ezr, ezi) = (zr * 8.0, zi * 8.0);

    let aez = 8.0 * az;
    let s = tol / aez;
    let jl = (rl + rl).trunc() as i32 + 2;

    let (mut p1r, mut p1i) = (zeror, zeroi);
    if zi != 0.0 {
        // CALCULATE EXP(PI*(0.5+FNU+N-IL)*I)
        let mut inu = fnu.trunc() as i32;
        let arg = (fnu - inu as f64) * PI;
        inu = inu + n - il;
        let ak = -arg.sin();
        let bk = if zi < 0.0 { -arg.cos() } else { arg.cos() };
        (p1r, p1i) = (ak, bk);
        if inu % 2 != 0 {
            (p1r, p1i) = (-p1r, -p1i);
        }
    }

    for k in 1..=il {
        let mut sqk = fdn - 1.0;
        let atol = s * sqk.abs();
        let mut sgn = 1.0_f64;
        let (mut cs1r, mut cs1i) = (coner, conei);
        let (mut cs2r, mut cs2i) = (coner, conei);
        let (mut ckr, mut cki)   = (coner, conei);
        let mut ak = 0.0_f64;
        let mut aa = 1.0_f64;
        let mut bb = aez;
        let (mut dkr, mut dki) = (ezr, ezi);

        let mut converged = false;
        for _ in 0..jl {
            zdiv(ckr, cki, dkr, dki, &mut str, &mut sti);
            (ckr, cki) = (str * sqk, sti * sqk);
            (cs2r, cs2i) = (cs2r + ckr, cs2i + cki);
            sgn = -sgn;
            (cs1r, cs1i) = (cs1r + ckr * sgn, cs1i + cki * sgn);
            (dkr, dki) = (dkr + ezr, dki + ezi);
            aa = aa * sqk.abs() / bb;
            bb += aez;
            ak += 8.0;
            sqk -= ak;
            if aa <= atol {
                converged = true;
                break;
            }
        }

        if !converged {
            *nz = -2;
            return;
        }

        let (mut s2r, mut s2i) = (cs1r, cs1i);

        if zr + zr < elim {
            let (tzr, tzi) = (zr + zr, zi + zi);
            zexp(-tzr, -tzi, &mut str, &mut sti);
            zmlt(str, sti, p1r, p1i, &mut str, &mut sti);
            zmlt(str, sti, cs2r, cs2i, &mut str, &mut sti);
            s2r += str;
            s2i += sti;
        }

        fdn = fdn + 8.0 * dfnu + 4.0;
        p1r = -p1r;
        p1i = -p1i;

        let m = (n - il + k) as usize;
        yr[m - 1] = s2r * ak1r - s2i * ak1i;
        yi[m - 1] = s2r * ak1i + s2i * ak1r;
    }

    if n <= 2 {
        return;
    }

    let nn = n;
    let mut k = nn - 2;
    let mut ak = k as f64;
    (str, sti) = (zr * raz, -zi * raz);
    let (rzr, rzi) = ((str + str) * raz, (sti + sti) * raz);

    for _ in 3..=nn {
        let kk = k as usize; // 1-based index into YR/YI
        yr[kk - 1] = (ak + fnu) * (rzr * yr[kk] - rzi * yi[kk]) + yr[kk + 1];
        yi[kk - 1] = (ak + fnu) * (rzr * yi[kk] + rzi * yr[kk]) + yi[kk + 1];
        ak -= 1.0;
        k -= 1;
    }

    if koded == 0 {
        return;
    }

    let (mut ckr, mut cki) = (0.0_f64, 0.0_f64);
    zexp(czr, czi, &mut ckr, &mut cki);

    for i in 0..nn as usize {
        str = yr[i] * ckr - yi[i] * cki;
        yi[i] = yr[i] * cki + yi[i] * ckr;
        yr[i] = str;
    }

    return;
}