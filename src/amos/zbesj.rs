use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zbinu::zbinu;

/// ZBESJ computes the J-Bessel function of a complex argument.
///
/// On KODE=1, ZBESJ computes an N-member sequence of complex Bessel
/// functions CY(I)=J(FNU+I-1,Z) for real nonnegative orders and
/// complex Z in the cut plane -PI < ARG(Z) <= PI. On KODE=2, ZBESJ
/// returns the scaled functions
///
/// CY(I) = EXP(-ABS(Y))*J(FNU+I-1,Z),  I=1,...,N,  Y=AIMAG(Z)
///
/// which remove the exponential growth in both half planes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbesj(
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
    const HPI: f64 = 1.57079632679489662;

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

    // -------------------------------------------------------------------
    // Calculate CSGN = EXP(FNU*HPI*I) to minimize losses of significance
    // when FNU is large.
    // -------------------------------------------------------------------
    let mut cii = 1.0;
    let inu = fnu as i32;
    let inuh = inu / 2;
    let ir = inu - 2 * inuh;
    let arg = (fnu - ((inu - ir) as f64)) * HPI;
    let mut csgnr = arg.cos();
    let mut csgni = arg.sin();
    if inuh % 2 != 0 {
        csgnr = -csgnr;
        csgni = -csgni;
    }

    // -------------------------------------------------------------------
    // ZN is in the right half plane.
    // -------------------------------------------------------------------
    let mut znr = zi;
    let mut zni = -zr;
    if zi < 0.0 {
        znr = -znr;
        zni = -zni;
        csgni = -csgni;
        cii = -cii;
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

    let nl = n - *nz;
    if nl == 0 {
        return;
    }

    let rtol = 1.0 / tol;
    let ascle = d1mach(1) * rtol * 1.0e+3;

    for i in 0..(nl as usize) {
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
        let str = -csgni * cii;
        csgni = csgnr * cii;
        csgnr = str;
    }
}