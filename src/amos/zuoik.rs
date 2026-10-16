use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zlog::zlog;
use crate::amos::zuchk::zuchk;
use crate::amos::zunhj::zunhj;
use crate::amos::zunik::zunik;

/// ZUOIK computes the leading terms of the uniform asymptotic
/// expansions for the I and K functions and compares them
/// (in logarithmic form) to ALIM and ELIM for over and underflow
/// where ALIM < ELIM. If the magnitude, based on the leading
/// exponential, is less than ALIM or greater than -ALIM, then
/// the result is on scale. If not, then a refined test using other
/// multipliers (in logarithmic form) is made based on ELIM. Here
/// exp(-ELIM) = smallest machine number*1.0E+3 and
/// exp(-ALIM) = exp(-ELIM)/TOL.
///
/// IKFLG=1 means the I sequence is tested
///      =2 means the K sequence is tested
/// NUF = 0 means the last member of the sequence is on scale
///     =-1 means an overflow would occur
/// IKFLG=1 and NUF>0 means the last NUF Y values were set to zero;
///         the first N-NUF values must be set by another routine
/// IKFLG=2 and NUF=N means all Y values were set to zero
/// IKFLG=2 and 0<NUF<N not considered. Y must be set by
///         another routine
#[allow(clippy::too_many_arguments)]
pub(crate) fn zuoik(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    ikflg: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nuf: &mut i32,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const AIC: f64 = 1.265512123484645396e+00;

    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nuf = 0;
    let mut nn = n;

    let (zrr, zri) = if zr < 0.0 {
        (-zr, -zi)
    } else {
        (zr, zi)
    };

    let zbr = zrr;
    let zbi = zri;
    let ax = zr.abs() * 1.7321;
    let ay = zi.abs();
    let iform = if ay > ax { 2 } else { 1 };
    let mut gnu = fnu.max(1.0);
    if ikflg != 1 {
        let fnn = nn as f64;
        let gnn = fnu + fnn - 1.0;
        gnu = gnn.max(fnn);
    }

    // -------------------------------------------------------------------
    // Only the magnitude of ARG and PHI are needed along with the
    // real parts of ZETA1, ZETA2 and ZB. No attempt is made to get
    // the sign of the imaginary part correct.
    // -------------------------------------------------------------------
    let mut phir = 0.0f64;
    let mut phii = 0.0f64;
    let mut argr = 0.0f64;
    let mut argi = 0.0f64;
    let mut zeta1r = 0.0f64;
    let mut zeta1i = 0.0f64;
    let mut zeta2r = 0.0f64;
    let mut zeta2i = 0.0f64;
    let mut asumr = 0.0f64;
    let mut asumi = 0.0f64;
    let mut bsumr = 0.0f64;
    let mut bsumi = 0.0f64;
    let mut sumr = 0.0f64;
    let mut sumi = 0.0f64;
    let mut cwrkr = [0.0f64; 16];
    let mut cwrki = [0.0f64; 16];

    // ZNR/ZNI are only used on the ZUNHJ path.
    let mut znr = zri;
    let zni = -zrr;

    let mut czr;
    let mut czi;
    let mut aarg = 0.0f64;

    if iform != 2 {
        let mut init: i32 = 0;
        zunik(
            zrr, zri, gnu, ikflg, 1, tol, &mut init, &mut phir, &mut phii, &mut zeta1r,
            &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut sumr, &mut sumi, &mut cwrkr,
            &mut cwrki,
        );
        czr = -zeta1r + zeta2r;
        czi = -zeta1i + zeta2i;
    } else {
        if zi > 0.0 {
            znr = -znr;
        }
        zunhj(
            znr, zni, gnu, 1, tol, &mut phir, &mut phii, &mut argr, &mut argi, &mut zeta1r,
            &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut asumr, &mut asumi, &mut bsumr,
            &mut bsumi,
        );
        czr = -zeta1r + zeta2r;
        czi = -zeta1i + zeta2i;
        aarg = zabs(argr, argi);
    }

    if kode != 1 {
        czr -= zbr;
        czi -= zbi;
    }

    if ikflg != 1 {
        czr = -czr;
        czi = -czi;
    }

    let aphi = zabs(phir, phii);
    let mut rcz = czr;

    // -------------------------------------------------------------------
    // Overflow test
    // -------------------------------------------------------------------
    if rcz > elim {
        *nuf = -1;
        return;
    }

    if rcz >= alim {
        rcz += aphi.ln();
        if iform == 2 {
            rcz = rcz - 0.25 * aarg.ln() - AIC;
        }
        if rcz > elim {
            *nuf = -1;
            return;
        }
        zuoik_continue(
            ikflg, n, &mut nn, nuf, yr, yi, fnu, kode, iform, tol, elim, alim, zbr, zbi,
            zrr, zri, znr, zni,
        );
        return;
    }

    // -------------------------------------------------------------------
    // Underflow test
    // -------------------------------------------------------------------
    if rcz <= -elim {
        for i in 0..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        *nuf = nn;
        return;
    }

    if rcz > -alim {
        // On scale.
        zuoik_continue(
            ikflg, n, &mut nn, nuf, yr, yi, fnu, kode, iform, tol, elim, alim, zbr, zbi,
            zrr, zri, znr, zni,
        );
        return;
    }

    rcz += aphi.ln();
    if iform == 2 {
        rcz = rcz - 0.25 * aarg.ln() - AIC;
    }
    if rcz <= -elim {
        for i in 0..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        *nuf = nn;
        return;
    }

    // Refined underflow test (label 110).
    let ascle = 1.0e+3 * d1mach(1) / tol;
    let (mut str, mut sti) = (0.0f64, 0.0f64);
    let _ = zlog(phir, phii, &mut str, &mut sti);
    czr += str;
    czi += sti;
    if iform != 1 {
        let _ = zlog(argr, argi, &mut str, &mut sti);
        czr = czr - 0.25 * str - AIC;
        czi -= 0.25 * sti;
    }
    let ax = rcz.exp() / tol;
    let ay = czi;
    czr = ax * ay.cos();
    czi = ax * ay.sin();
    let mut nw: i32 = 0;
    zuchk(czr, czi, &mut nw, ascle, tol);
    if nw != 0 {
        for i in 0..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        *nuf = nn;
        return;
    }

    // On scale.
    zuoik_continue(
        ikflg, n, &mut nn, nuf, yr, yi, fnu, kode, iform, tol, elim, alim, zbr, zbi, zrr,
        zri, znr, zni,
    );
}

// ---------------------------------------------------------------------
// The "on scale" continuation (Fortran label 130 onwards).
//
// In Fortran this was a shared tail reachable from several GOTOs, using
// host association for all state. In Rust we pass the mutable state
// explicitly.
// ---------------------------------------------------------------------
#[allow(clippy::too_many_arguments)]
fn zuoik_continue(
    ikflg: i32,
    n: i32,
    nn: &mut i32,
    nuf: &mut i32,
    yr: &mut [f64],
    yi: &mut [f64],
    fnu: f64,
    kode: i32,
    iform: i32,
    tol: f64,
    elim: f64,
    alim: f64,
    zbr: f64,
    zbi: f64,
    zrr: f64,
    zri: f64,
    znr: f64,
    zni: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const AIC: f64 = 1.265512123484645396e+00;

    if ikflg == 2 {
        return;
    }
    if n == 1 {
        return;
    }

    // -------------------------------------------------------------------
    // Set underflows on I sequence
    // -------------------------------------------------------------------
    let mut phir = 0.0f64;
    let mut phii = 0.0f64;
    let mut argr = 0.0f64;
    let mut argi = 0.0f64;
    let mut zeta1r = 0.0f64;
    let mut zeta1i = 0.0f64;
    let mut zeta2r = 0.0f64;
    let mut zeta2i = 0.0f64;
    let mut asumr = 0.0f64;
    let mut asumi = 0.0f64;
    let mut bsumr = 0.0f64;
    let mut bsumi = 0.0f64;
    let mut sumr = 0.0f64;
    let mut sumi = 0.0f64;
    let mut cwrkr = [0.0f64; 16];
    let mut cwrki = [0.0f64; 16];

    loop {
        let gnu = fnu + ((*nn - 1) as f64);

        let (czr, czi, aarg);
        if iform != 2 {
            let mut init: i32 = 0;
            zunik(
                zrr, zri, gnu, ikflg, 1, tol, &mut init, &mut phir, &mut phii, &mut zeta1r,
                &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut sumr, &mut sumi, &mut cwrkr,
                &mut cwrki,
            );
            czr = -zeta1r + zeta2r;
            czi = -zeta1i + zeta2i;
            aarg = 0.0;
        } else {
            if znr == 0.0 && zni == 0.0 {
                // keep compiler happy; znr/zni set by caller already
            }
            zunhj(
                znr, zni, gnu, 1, tol, &mut phir, &mut phii, &mut argr, &mut argi,
                &mut zeta1r, &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut asumr, &mut asumi,
                &mut bsumr, &mut bsumi,
            );
            czr = -zeta1r + zeta2r;
            czi = -zeta1i + zeta2i;
            aarg = zabs(argr, argi);
        }

        let mut czr = czr;
        let mut czi = czi;
        if kode != 1 {
            czr -= zbr;
            czi -= zbi;
        }

        let aphi = zabs(phir, phii);
        let mut rcz = czr;

        if rcz < -elim {
            // Underflow — zero and continue.
            yr[(*nn - 1) as usize] = ZEROR;
            yi[(*nn - 1) as usize] = ZEROI;
            *nn -= 1;
            *nuf += 1;
            if *nn == 0 {
                return;
            }
            continue;
        }

        if rcz > -alim {
            return;
        }

        rcz += aphi.ln();
        if iform == 2 {
            rcz = rcz - 0.25 * aarg.ln() - AIC;
        }
        if rcz <= -elim {
            yr[(*nn - 1) as usize] = ZEROR;
            yi[(*nn - 1) as usize] = ZEROI;
            *nn -= 1;
            *nuf += 1;
            if *nn == 0 {
                return;
            }
            continue;
        }

        // Refined underflow test (label 190).
        let ascle = 1.0e+3 * d1mach(1) / tol;
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        let _ = zlog(phir, phii, &mut str, &mut sti);
        czr += str;
        czi += sti;
        if iform != 1 {
            let _ = zlog(argr, argi, &mut str, &mut sti);
            czr = czr - 0.25 * str - AIC;
            czi -= 0.25 * sti;
        }
        let ax = rcz.exp() / tol;
        let ay = czi;
        czr = ax * ay.cos();
        czi = ax * ay.sin();
        let mut nw: i32 = 0;
        zuchk(czr, czi, &mut nw, ascle, tol);
        if nw != 0 {
            yr[(*nn - 1) as usize] = ZEROR;
            yi[(*nn - 1) as usize] = ZEROI;
            *nn -= 1;
            *nuf += 1;
            if *nn == 0 {
                return;
            }
            continue;
        }

        // On scale — done.
        return;
    }
}