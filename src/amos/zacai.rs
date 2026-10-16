use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zasyi::zasyi;
use crate::amos::zbknu::zbknu;
use crate::amos::zmlri::zmlri;
use crate::amos::zseri::zseri;
use crate::amos::zs1s2::zs1s2;

/// ZACAI applies the analytic continuation formula
///
///     K(FNU,ZN*EXP(MP)) = K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
///             MP = PI*MR*CMPLX(0.0,1.0)
///
/// to continue the K function from the right half to the left
/// half Z plane for use with ZAIRY where FNU=1/3 or 2/3 and N=1.
/// ZACAI is the same as ZACON with the parts for larger orders and
/// recurrence removed. A recursive call to ZACON can result if ZACON
/// is called from ZAIRY.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zacai(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    mr: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    rl: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    const PI: f64 = 3.14159265358979324;

    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;
    let znr = -zr;
    let zni = -zi;
    let az = zabs(zr, zi);
    let nn = n;
    let dfnu = fnu + ((n - 1) as f64);

    // -----------------------------------------------------------------
    // Dispatch: power series vs. large-|Z| branches
    // -----------------------------------------------------------------
    if az > 2.0 && az * az * 0.25 > dfnu + 1.0 {
        // Large-|Z| branch.
        let mut nw: i32 = 0;
        if az >= rl {
            // Asymptotic expansion for large Z for the I function.
            zasyi(znr, zni, fnu, kode, nn, yr, yi, &mut nw, rl, tol, elim, alim);
        } else {
            // Miller algorithm normalized by the series for the I function.
            zmlri(znr, zni, fnu, kode, nn, yr, yi, &mut nw, tol);
        }
        if nw < 0 {
            *nz = if nw == -2 { -2 } else { -1 };
            return;
        }
    } else {
        // Power series for the I function.
        let mut nw: i32 = 0;
        zseri(znr, zni, fnu, kode, nn, yr, yi, &mut nw, tol, elim, alim);
        // Note: NW is not checked after ZSERI, matching the Fortran.
    }

    // -----------------------------------------------------------------
    // Analytic continuation to the left half plane for the K function
    // -----------------------------------------------------------------
    let mut cyr = [0.0f64; 2];
    let mut cyi = [0.0f64; 2];
    let mut nw: i32 = 0;
    zbknu(znr, zni, fnu, kode, 1, &mut cyr, &mut cyi, &mut nw, tol, elim, alim);
    if nw != 0 {
        *nz = if nw == -2 { -2 } else { -1 };
        return;
    }

    let fmr = mr as f64;
    let sgn = -PI.copysign(fmr);
    let mut csgnr = 0.0f64;
    let mut csgni = sgn;

    if kode != 1 {
        let yy = -zni;
        csgnr = -csgni * yy.sin();
        csgni *= yy.cos();
    }

    // -----------------------------------------------------------------
    // Calculate CSPN = EXP(FNU*PI*I) to minimize loss of significance
    // when FNU is large.
    // -----------------------------------------------------------------
    let inu = fnu as i32;
    let arg = (fnu - (inu as f64)) * sgn;
    let mut cspnr = arg.cos();
    let mut cspni = arg.sin();

    if inu % 2 != 0 {
        cspnr = -cspnr;
        cspni = -cspni;
    }

    let mut c1r = cyr[0];
    let mut c1i = cyi[0];
    let mut c2r = yr[0];
    let mut c2i = yi[0];

    if kode != 1 {
        let mut iuf: i32 = 0;
        let ascle = 1.0e+3 * d1mach(1) / tol;
        let mut nw_s1s2: i32 = 0;
        zs1s2(znr, zni, &mut c1r, &mut c1i, &mut c2r, &mut c2i, &mut nw_s1s2, ascle, alim, &mut iuf);
        *nz += nw_s1s2;
    }

    yr[0] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
    yi[0] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
}