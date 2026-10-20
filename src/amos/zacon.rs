use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zbinu::zbinu;
use crate::amos::zbknu::zbknu;
use crate::amos::zmlt::zmlt;
use crate::amos::zs1s2::zs1s2;

/// ZACON applies the analytic continuation formula
///
/// K(FNU,ZN*EXP(MP)) = K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
/// MP = PI*MR*CMPLX(0.0,1.0)
///
/// to continue the K function from the right half to the left
/// half Z plane.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zacon(
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
    fnul: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    const PI: f64 = 3.14159265358979324;
    const ZEROR: f64 = 0.0;
    const CONER: f64 = 1.0;

    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;
    let znr = -zr;
    let zni = -zi;
    let mut nn = n;

    let mut nw: i32 = 0;
    zbinu(znr, zni, fnu, kode, nn, yr, yi, &mut nw, rl, fnul, tol, elim, alim);
    if nw < 0 {
        *nz = if nw == -2 { -2 } else { -1 };
        return;
    }

    // -------------------------------------------------------------------
    // Analytic continuation to the left half plane for the K function.
    // -------------------------------------------------------------------
    nn = 2.min(n);
    let mut cyr = [0.0f64; 2];
    let mut cyi = [0.0f64; 2];
    zbknu(znr, zni, fnu, kode, nn, &mut cyr, &mut cyi, &mut nw, tol, elim, alim);
    if nw != 0 {
        *nz = if nw == -2 { -2 } else { -1 };
        return;
    }

    let s1r = cyr[0];
    let s1i = cyi[0];
    let fmr = mr as f64;
    let sgn = -PI.copysign(fmr);
    let mut csgnr = ZEROR;
    let mut csgni = sgn;

    if kode != 1 {
        let yy = -zni;
        let cpn = yy.cos();
        let spn = yy.sin();
        let (mut outr, mut outi) = (0.0f64, 0.0f64);
        zmlt(csgnr, csgni, cpn, spn, &mut outr, &mut outi);
        csgnr = outr;
        csgni = outi;
    }

    // -------------------------------------------------------------------
    // Calculate CSPN = EXP(FNU*PI*I) to minimize losses of significance
    // when FNU is large.
    // -------------------------------------------------------------------
    let inu = fnu as i32;
    let arg = (fnu - (inu as f64)) * sgn;
    let cpn = arg.cos();
    let spn = arg.sin();
    let mut cspnr = cpn;
    let mut cspni = spn;
    if inu % 2 != 0 {
        cspnr = -cspnr;
        cspni = -cspni;
    }

    let mut iuf: i32 = 0;
    let mut c1r = s1r;
    let mut c1i = s1i;
    let mut c2r = yr[0];
    let mut c2i = yi[0];
    let ascle = 1.0e+3 * d1mach(1) / tol;

    let mut sc1r = 0.0f64;
    let mut sc1i = 0.0f64;
    let mut sc2r = 0.0f64;
    let mut sc2i = 0.0f64;

    if kode != 1 {
        zs1s2(znr, zni, &mut c1r, &mut c1i, &mut c2r, &mut c2i, &mut nw, ascle, alim, &mut iuf);
        *nz += nw;
        sc1r = c1r;
        sc1i = c1i;
    }

    let (mut str, mut sti) = (0.0f64, 0.0f64);
    zmlt(cspnr, cspni, c1r, c1i, &mut str, &mut sti);
    let (mut ptr, mut pti) = (0.0f64, 0.0f64);
    zmlt(csgnr, csgni, c2r, c2i, &mut ptr, &mut pti);
    yr[0] = str + ptr;
    yi[0] = sti + pti;

    if n == 1 {
        return;
    }

    cspnr = -cspnr;
    cspni = -cspni;
    let mut s2r = cyr[1];
    let mut s2i = cyi[1];
    c1r = s2r;
    c1i = s2i;
    c2r = yr[1];
    c2i = yi[1];

    if kode != 1 {
        zs1s2(znr, zni, &mut c1r, &mut c1i, &mut c2r, &mut c2i, &mut nw, ascle, alim, &mut iuf);
        *nz += nw;
        sc2r = c1r;
        sc2i = c1i;
    }

    zmlt(cspnr, cspni, c1r, c1i, &mut str, &mut sti);
    zmlt(csgnr, csgni, c2r, c2i, &mut ptr, &mut pti);
    yr[1] = str + ptr;
    yi[1] = sti + pti;

    if n == 2 {
        return;
    }

    cspnr = -cspnr;
    cspni = -cspni;
    let azn = zabs(znr, zni);
    let razn = 1.0 / azn;
    let str = znr * razn;
    let sti = -zni * razn;
    let rzr = (str + str) * razn;
    let rzi = (sti + sti) * razn;
    let fn_val = fnu + 1.0;
    let mut ckr = fn_val * rzr;
    let mut cki = fn_val * rzi;

    // -------------------------------------------------------------------
    // Scale near exponent extremes during recurrence on K functions.
    // -------------------------------------------------------------------
    let cscl = 1.0 / tol;
    let cscr = tol;
    let cssr = [cscl, CONER, cscr];
    let csrr = [cscr, CONER, cscl];
    let bry = [ascle, 1.0 / ascle, d1mach(2)];

    let as2 = zabs(s2r, s2i);
    let mut kflag: i32 = 2;
    if as2 <= bry[0] {
        kflag = 1;
    } else if as2 >= bry[1] {
        kflag = 3;
    }

    let mut bscle = bry[(kflag - 1) as usize];
    let mut s1r = s1r * cssr[(kflag - 1) as usize];
    let mut s1i = s1i * cssr[(kflag - 1) as usize];
    s2r *= cssr[(kflag - 1) as usize];
    s2i *= cssr[(kflag - 1) as usize];
    let mut csr = csrr[(kflag - 1) as usize];

    let mut i = 3;
    while i <= n {
        let str = s2r;
        let sti = s2i;
        s2r = ckr * str - cki * sti + s1r;
        s2i = ckr * sti + cki * str + s1i;
        s1r = str;
        s1i = sti;
        c1r = s2r * csr;
        c1i = s2i * csr;
        let str = c1r;
        let sti = c1i;
        c2r = yr[(i - 1) as usize];
        c2i = yi[(i - 1) as usize];

        if kode != 1 && iuf >= 0 {
            zs1s2(znr, zni, &mut c1r, &mut c1i, &mut c2r, &mut c2i, &mut nw, ascle, alim, &mut iuf);
            *nz += nw;
            sc1r = sc2r;
            sc1i = sc2i;
            sc2r = c1r;
            sc2i = c1i;
            if iuf == 3 {
                iuf = -4;
                s1r = sc1r * cssr[(kflag - 1) as usize];
                s1i = sc1i * cssr[(kflag - 1) as usize];
                s2r = sc2r * cssr[(kflag - 1) as usize];
                s2i = sc2i * cssr[(kflag - 1) as usize];
                // Note: `str`/`sti` were already destructured above; the
                // Fortran reassigns them here but they're only read in the
                // KFLAG bump block below, which we handle inline.
            }
        }

        let ptr = cspnr * c1r - cspni * c1i;
        let pti = cspnr * c1i + cspni * c1r;
        yr[(i - 1) as usize] = ptr + csgnr * c2r - csgni * c2i;
        yi[(i - 1) as usize] = pti + csgnr * c2i + csgni * c2r;
        ckr += rzr;
        cki += rzi;
        cspnr = -cspnr;
        cspni = -cspni;

        if kflag < 3 {
            let c1m = c1r.abs().max(c1i.abs());
            if c1m > bscle {
                kflag += 1;
                bscle = bry[(kflag - 1) as usize];
                s1r *= csr;
                s1i *= csr;
                s2r = str;
                s2i = sti;
                s1r *= cssr[(kflag - 1) as usize];
                s1i *= cssr[(kflag - 1) as usize];
                s2r *= cssr[(kflag - 1) as usize];
                s2i *= cssr[(kflag - 1) as usize];
                csr = csrr[(kflag - 1) as usize];
            }
        }
        i += 1;
    }
}