use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zbesh::zbesh;

/// ZBESY computes the Y-Bessel function of a complex argument.
///
/// On KODE=1, ZBESY computes an N-member sequence of complex Bessel
/// functions CY(I)=Y(FNU+I-1,Z) for real nonnegative orders and
/// complex Z in the cut plane -PI < ARG(Z) <= PI.
///
/// On KODE=2, ZBESY returns the scaled functions
///
/// CY(I) = EXP(-ABS(Y))*Y(FNU+I-1,Z),  I=1,...,N,  Y=AIMAG(Z)
///
/// which remove the exponential growth in both half planes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbesy(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    nz: &mut i32,
    cwrkr: &mut [f64],
    cwrki: &mut [f64],
    ierr: &mut i32,
) {
    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");
    assert!(cwrkr.len() >= n_usize, "cwrkr must have at least n elements");
    assert!(cwrki.len() >= n_usize, "cwrki must have at least n elements");

    *ierr = 0;
    *nz = 0;

    if zr == 0.0 && zi == 0.0 {
        *ierr = 1;
    }
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

    let hcii = 0.5;

    let mut nz1: i32 = 0;
    zbesh(zr, zi, fnu, kode, 1, n, cyr, cyi, &mut nz1, ierr);
    if *ierr != 0 && *ierr != 3 {
        *nz = 0;
        return;
    }

    let mut nz2: i32 = 0;
    zbesh(zr, zi, fnu, kode, 2, n, cwrkr, cwrki, &mut nz2, ierr);
    if *ierr != 0 && *ierr != 3 {
        *nz = 0;
        return;
    }

    *nz = nz1.min(nz2);

    if kode != 2 {
        // Unscaled path.
        for i in 0..n_usize {
            let str = cwrkr[i] - cyr[i];
            let sti = cwrki[i] - cyi[i];
            cyr[i] = -sti * hcii;
            cyi[i] = str * hcii;
        }
        return;
    }

    // -------------------------------------------------------------------
    // KODE=2 scaled path.
    // -------------------------------------------------------------------
    let tol = d1mach(4).max(1.0e-18);
    let k1 = i1mach(15);
    let k2 = i1mach(16);
    let k = k1.abs().min(k2.abs());
    let r1m5 = d1mach(5);

    // ELIM is the approximate exponential under- and overflow limit.
    let elim = 2.303 * ((k as f64) * r1m5 - 3.0);
    let exr = zr.cos();
    let exi = zr.sin();
    let mut ey = 0.0;
    let tay = (zi + zi).abs();
    if tay < elim {
        ey = (-tay).exp();
    }

    let (c1r, c1i, c2r, c2i);
    if zi < 0.0 {
        // Lower half plane.
        c1r = exr;
        c1i = exi;
        c2r = exr * ey;
        c2i = -exi * ey;
    } else {
        // Upper half plane.
        c1r = exr * ey;
        c1i = exi * ey;
        c2r = exr;
        c2i = -exi;
    }

    // Scaling loop.
    *nz = 0;
    let rtol = 1.0 / tol;
    let ascle = d1mach(1) * rtol * 1.0e+3;

    for i in 0..n_usize {
        let mut aa = cwrkr[i];
        let mut bb = cwrki[i];
        let mut atol = 1.0;
        if aa.abs().max(bb.abs()) <= ascle {
            aa *= rtol;
            bb *= rtol;
            atol = tol;
        }
        let mut str = (aa * c2r - bb * c2i) * atol;
        let mut sti = (aa * c2i + bb * c2r) * atol;

        aa = cyr[i];
        bb = cyi[i];
        atol = 1.0;
        if aa.abs().max(bb.abs()) <= ascle {
            aa *= rtol;
            bb *= rtol;
            atol = tol;
        }
        str -= (aa * c1r - bb * c1i) * atol;
        sti -= (aa * c1i + bb * c1r) * atol;

        cyr[i] = -sti * hcii;
        cyi[i] = str * hcii;

        if str == 0.0 && sti == 0.0 && ey == 0.0 {
            *nz += 1;
        }
    }
}