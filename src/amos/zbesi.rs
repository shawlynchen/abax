use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zbinu::zbinu;

/// ZBESI computes I-Bessel functions of complex argument.
///
/// On KODE=1, ZBESI computes an N-member sequence of complex Bessel
/// functions CY(J)=I(FNU+J-1,Z) for real nonnegative orders and
/// complex Z in the cut plane -PI < ARG(Z) <= PI. On KODE=2, ZBESI
/// returns the scaled functions
///
///     CY(J) = EXP(-ABS(X))*I(FNU+J-1,Z),  J=1,...,N,  X=REAL(Z)
///
/// with the exponential growth removed in both half planes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbesi(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    nz: &mut i32,
    ierr: &mut i32,
) {
    const PI: f64 = 3.14159265358979324;
    const CONER: f64 = 1.0;
    const CONEI: f64 = 0.0;

    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");

    *ierr = 0;
    *nz = 0;

    if fnu < 0.0 {
        *ierr = 1;
    }
    if kode < 1 || kode > 2 {
        *ierr = 1;
    }
    if n < 1 {
        *ierr = 1;
    }
    if *ierr != 0 {
        return;
    }

    // -------------------------------------------------------------------
    // Set parameters related to machine constants.
    // -------------------------------------------------------------------
    let tol = d1mach(4).max(1.0e-18);
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
    // Test for proper range.
    // -------------------------------------------------------------------
    let az = zabs(zr, zi);
    let fn_val = fnu + ((n - 1) as f64);
    let mut aa = 0.5 / tol;
    let bb = (i1mach(9) as f64) * 0.5;
    aa = aa.min(bb);
    if az > aa {
        *nz = 0;
        *ierr = 4;
        return;
    }
    if fn_val > aa {
        *nz = 0;
        *ierr = 4;
        return;
    }
    aa = aa.sqrt();
    if az > aa {
        *ierr = 3;
    }
    if fn_val > aa {
        *ierr = 3;
    }

    let mut znr = zr;
    let mut zni = zi;
    let mut csgnr = CONER;
    let mut csgni = CONEI;

    if zr < 0.0 {
        znr = -zr;
        zni = -zi;

        // -------------------------------------------------------------
        // Calculate CSGN = EXP(FNU*PI*I) to minimize losses of
        // significance when FNU is large.
        // -------------------------------------------------------------
        let inu = fnu as i32;
        let mut arg = (fnu - (inu as f64)) * PI;
        if zi < 0.0 {
            arg = -arg;
        }
        csgnr = arg.cos();
        csgni = arg.sin();
        if inu % 2 != 0 {
            csgnr = -csgnr;
            csgni = -csgni;
        }
    }

    zbinu(znr, zni, fnu, kode, n, cyr, cyi, nz, rl, fnul, tol, elim, alim);

    if *nz < 0 {
        if *nz == -2 {
            *nz = 0;
            *ierr = 5;
            return;
        }
        *nz = 0;
        *ierr = 2;
        return;
    }

    if zr >= 0.0 {
        return;
    }

    // -------------------------------------------------------------------
    // Analytic continuation to the left half plane.
    // -------------------------------------------------------------------
    let nn = n - *nz;
    if nn == 0 {
        return;
    }

    let rtol = 1.0 / tol;
    let ascle = d1mach(1) * rtol * 1.0e+3;

    for i in 0..(nn as usize) {
        let mut a = cyr[i];
        let mut b = cyi[i];
        let mut atol = 1.0;
        if a.abs().max(b.abs()) <= ascle {
            a *= rtol;
            b *= rtol;
            atol = tol;
        }
        let str = a * csgnr - b * csgni;
        let sti = a * csgni + b * csgnr;
        cyr[i] = str * atol;
        cyi[i] = sti * atol;
        csgnr = -csgnr;
        csgni = -csgni;
    }
}