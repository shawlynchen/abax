use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zbinu::zbinu;
use crate::amos::zdiv::zdiv;
use crate::amos::zsqrt::zsqrt;

// ---------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------

const CONER: f64 = 1.0;
const CONEI: f64 = 0.0;

/// ZBIRY computes Airy functions BI(Z) and DBI(Z) for complex Z.
///
/// On KODE=1, ZBIRY computes BI(Z) or its derivative DBI(Z)/DZ on
/// ID=0 or ID=1 respectively. On KODE=2, a scaling option
/// CEXP(-AXZTA)*BI(Z) or CEXP(-AXZTA)*DBI(Z)/DZ removes the
/// exponential behavior in both half planes where
/// ZTA = (2/3)*Z*CSQRT(Z) and AXZTA = ABS(RE(ZTA)).
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbiry(
    zr: f64,
    zi: f64,
    id: i32,
    kode: i32,
    bir: &mut f64,
    bii: &mut f64,
    ierr: &mut i32,
) {
    const TTH: f64 = 6.66666666666666667e-01;
    const C1: f64 = 6.14926627446000736e-01;
    const C2: f64 = 4.48288357353826359e-01;
    const COEF: f64 = 5.77350269189625765e-01;
    const PI: f64 = 3.14159265358979324;

    *ierr = 0;
    let mut nz: i32 = 0;

    if id < 0 || id > 1 {
        *ierr = 1;
    }
    if kode < 1 || kode > 2 {
        *ierr = 1;
    }
    if *ierr != 0 {
        return;
    }

    let az = zabs(zr, zi);
    let tol = d1mach(4).max(1.0e-18);
    let fid = id as f64;

    if az <= 1.0 {
        // ----------------------------------------------------------------
        // Power series for |Z| <= 1.
        // ----------------------------------------------------------------
        let mut s1r = CONER;
        let mut s1i = CONEI;
        let mut s2r = CONER;
        let mut s2i = CONEI;

        if az < tol {
            // Extremely small |Z|.
            let aa = C1 * (1.0 - fid) + fid * C2;
            *bir = aa;
            *bii = 0.0;
            return;
        }

        let aa = az * az;
        if aa >= tol / az {
            // Series accumulation.
            let mut trm1r = CONER;
            let mut trm1i = CONEI;
            let mut trm2r = CONER;
            let mut trm2i = CONEI;
            let mut atrm = 1.0;
            let mut str = zr * zr - zi * zi;
            let mut sti = zr * zi + zi * zr;
            let z3r = str * zr - sti * zi;
            let z3i = str * zi + sti * zr;
            let az3 = az * aa;
            let mut ak = 2.0 + fid;
            let mut bk = 3.0 - fid - fid;
            let ck = 4.0 - fid;
            let _ = ck;
            let dk = 3.0 + fid + fid;
            let mut d1 = ak * dk;
            let mut d2 = bk * ck;
            let mut ad = d1.min(d2);
            ak = 24.0 + 9.0 * fid;
            bk = 30.0 - 9.0 * fid;

            for _k in 1..=25 {
                str = (trm1r * z3r - trm1i * z3i) / d1;
                trm1i = (trm1r * z3i + trm1i * z3r) / d1;
                trm1r = str;
                s1r += trm1r;
                s1i += trm1i;
                str = (trm2r * z3r - trm2i * z3i) / d2;
                trm2i = (trm2r * z3i + trm2i * z3r) / d2;
                trm2r = str;
                s2r += trm2r;
                s2i += trm2i;
                atrm = atrm * az3 / ad;
                d1 += ak;
                d2 += bk;
                ad = d1.min(d2);
                if atrm < tol * ad {
                    break;
                }
                ak += 18.0;
                bk += 18.0;
            }
        }

        // Series end.
        if id != 1 {
            *bir = C1 * s1r + C2 * (zr * s2r - zi * s2i);
            *bii = C1 * s1i + C2 * (zr * s2i + zi * s2r);
            if kode == 1 {
                return;
            }
            let (mut str, mut sti) = (0.0f64, 0.0f64);
            zsqrt(zr, zi, &mut str, &mut sti);
            let ztar = TTH * (zr * str - zi * sti);
            let _ztai = TTH * (zr * sti + zi * str);
            let aa = -ztar.abs();
            let eaa = aa.exp();
            *bir *= eaa;
            *bii *= eaa;
            return;
        }

        // ID=1 branch.
        *bir = s2r * C2;
        *bii = s2i * C2;
        if az > tol {
            let cc = C1 / (1.0 + fid);
            let str = s1r * zr - s1i * zi;
            let sti = s1r * zi + s1i * zr;
            *bir += cc * (str * zr - sti * zi);
            *bii += cc * (str * zi + sti * zr);
        }

        if kode == 1 {
            return;
        }
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        zsqrt(zr, zi, &mut str, &mut sti);
        let ztar = TTH * (zr * str - zi * sti);
        let _ztai = TTH * (zr * sti + zi * str);
        let aa = -ztar.abs();
        let eaa = aa.exp();
        *bir *= eaa;
        *bii *= eaa;
        return;
    }

    // -------------------------------------------------------------------
    // Case for |Z| > 1.0
    // -------------------------------------------------------------------
    let mut fnu = (1.0 + fid) / 3.0;

    // -------------------------------------------------------------------
    // Set parameters related to machine constants.
    // -------------------------------------------------------------------
    let k1 = i1mach(15);
    let k2 = i1mach(16);
    let r1m5 = d1mach(5);
    let k = k1.abs().min(k2.abs());
    let elim = 2.303 * ((k as f64) * r1m5 - 3.0);
    let k1b = i1mach(14) - 1;
    let mut aa = r1m5 * (k1b as f64);
    let dig = aa.min(18.0);
    aa *= 2.303;
    let alim = elim + (-aa).max(-41.45);
    let rl = 1.2 * dig + 3.0;
    let fnul = 10.0 + 6.0 * (dig - 3.0);

    // -------------------------------------------------------------------
    // Test for range.
    // -------------------------------------------------------------------
    aa = 0.5 / tol;
    let bb = (i1mach(9) as f64) * 0.5;
    aa = aa.min(bb);
    aa = aa.powf(TTH);
    if az > aa {
        *ierr = 4;
        nz = 0;
        return;
    }
    aa = aa.sqrt();
    if az > aa {
        *ierr = 3;
    }

    let (mut csqr, mut csqi) = (0.0f64, 0.0f64);
    zsqrt(zr, zi, &mut csqr, &mut csqi);
    let mut ztar = TTH * (zr * csqr - zi * csqi);
    let mut ztai = TTH * (zr * csqi + zi * csqr);

    // -------------------------------------------------------------------
    // Re(ZTA) <= 0 when Re(Z) < 0, especially when Im(Z) is small.
    // -------------------------------------------------------------------
    let mut sfac = 1.0;
    let ak = ztai;
    if zr < 0.0 {
        let bk = ztar;
        let ck = -bk.abs();
        ztar = ck;
        ztai = ak;
    }

    if zi == 0.0 && zr <= 0.0 {
        ztar = 0.0;
        ztai = ak;
    }

    aa = ztar;

    if kode != 2 {
        // Overflow test.
        let bb = aa.abs();
        if bb >= alim {
            let bb2 = bb + 0.25 * az.ln();
            sfac = tol;
            if bb2 > elim {
                *ierr = 2;
                nz = 0;
                return;
            }
        }
    }

    let mut fmr = 0.0;
    if !(aa >= 0.0 && zr > 0.0) {
        fmr = PI;
        if zi < 0.0 {
            fmr = -PI;
        }
        ztar = -ztar;
        ztai = -ztai;
    }

    // -------------------------------------------------------------------
    // AA = factor for analytic continuation of I(FNU,ZTA).
    // KODE=2 returns EXP(-ABS(XZTA))*I(FNU,ZTA) from ZBINU.
    // -------------------------------------------------------------------
    let mut cyr = [0.0f64; 2];
    let mut cyi = [0.0f64; 2];
    zbinu(ztar, ztai, fnu, kode, 1, &mut cyr, &mut cyi, &mut nz, rl, fnul, tol, elim, alim);
    if nz < 0 {
        if nz == -1 {
            *ierr = 2;
            nz = 0;
            return;
        }
        nz = 0;
        *ierr = 5;
        return;
    }

    aa = fmr * fnu;
    let z3r = sfac;
    let str = aa.cos();
    let sti = aa.sin();
    let mut s1r = (str * cyr[0] - sti * cyi[0]) * z3r;
    let mut s1i = (str * cyi[0] + sti * cyr[0]) * z3r;

    fnu = (2.0 - fid) / 3.0;
    zbinu(ztar, ztai, fnu, kode, 2, &mut cyr, &mut cyi, &mut nz, rl, fnul, tol, elim, alim);

    cyr[0] *= z3r;
    cyi[0] *= z3r;
    cyr[1] *= z3r;
    cyi[1] *= z3r;

    // -------------------------------------------------------------------
    // Backward recur one step for orders -1/3 or -2/3.
    // -------------------------------------------------------------------
    let (mut str, mut sti) = (0.0f64, 0.0f64);
    zdiv(cyr[0], cyi[0], ztar, ztai, &mut str, &mut sti);
    let s2r = (fnu + fnu) * str + cyr[1];
    let s2i = (fnu + fnu) * sti + cyi[1];
    aa = fmr * (fnu - 1.0);
    let str = aa.cos();
    let sti = aa.sin();
    s1r = COEF * (s1r + s2r * str - s2i * sti);
    s1i = COEF * (s1i + s2r * sti + s2i * str);

    if id != 1 {
        let str = csqr * s1r - csqi * s1i;
        s1i = csqr * s1i + csqi * s1r;
        s1r = str;
        *bir = s1r / sfac;
        *bii = s1i / sfac;
        return;
    }

    let str = zr * s1r - zi * s1i;
    s1i = zr * s1i + zi * s1r;
    s1r = str;
    *bir = s1r / sfac;
    *bii = s1i / sfac;
}