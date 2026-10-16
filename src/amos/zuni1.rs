use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zuchk::zuchk;
use crate::amos::zunik::zunik;
use crate::amos::zuoik::zuoik;

/// ZUNI1 computes I(FNU,Z) by means of the uniform asymptotic
/// expansion for I(FNU,Z) in -PI/3 <= ARG Z <= PI/3.
///
/// FNUL is the smallest order permitted for the asymptotic
/// expansion. NLAST=0 means all of the Y values were set.
/// NLAST/=0 is the number left to be computed by another
/// formula for orders FNU to FNU+NLAST-1 because FNU+NLAST-1 < FNUL.
/// Y(I)=CZERO for I=NLAST+1,N
#[allow(clippy::too_many_arguments)]
pub(crate) fn zuni1(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    nlast: &mut i32,
    fnul: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const CONER: f64 = 1.0;

    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;
    let mut nd = n;
    *nlast = 0;

    // -------------------------------------------------------------------
    // Computed values with exponents between ALIM and ELIM in magnitude
    // are scaled to keep intermediate arithmetic on scale,
    // exp(ALIM) = exp(ELIM)*TOL
    // -------------------------------------------------------------------
    let cscl = 1.0 / tol;
    let crsc = tol;
    let cssr = [cscl, CONER, crsc];
    let csrr = [crsc, CONER, cscl];
    let mut bry = [0.0f64; 3];
    bry[0] = 1.0e+3 * d1mach(1) / tol;

    // -------------------------------------------------------------------
    // Check for underflow and overflow on first member
    // -------------------------------------------------------------------
    let mut fn_val = fnu.max(1.0);
    let mut init: i32 = 0;

    let mut phir = 0.0f64;
    let mut phii = 0.0f64;
    let mut zeta1r = 0.0f64;
    let mut zeta1i = 0.0f64;
    let mut zeta2r = 0.0f64;
    let mut zeta2i = 0.0f64;
    let mut sumr = 0.0f64;
    let mut sumi = 0.0f64;
    let mut cwrkr = [0.0f64; 16];
    let mut cwrki = [0.0f64; 16];

    zunik(
        zr, zi, fn_val, 1, 1, tol, &mut init, &mut phir, &mut phii, &mut zeta1r,
        &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut sumr, &mut sumi, &mut cwrkr, &mut cwrki,
    );

    let (mut s1r, mut s1i);
    if kode != 1 {
        let mut str = zr + zeta2r;
        let mut sti = zi + zeta2i;
        let rast = fn_val / zabs(str, sti);
        str = str * rast * rast;
        sti = -sti * rast * rast;
        s1r = -zeta1r + str;
        s1i = -zeta1i + sti;
    } else {
        s1r = -zeta1r + zeta2r;
        s1i = -zeta1i + zeta2i;
    }

    let mut rs1 = s1r;

    if rs1.abs() > elim {
        if rs1 > 0.0 {
            *nz = -1;
            return;
        }
        *nz = n;
        for i in 0..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        return;
    }

    let mut iflag: i32 = 0;
    let mut cyr = [0.0f64; 2];
    let mut cyi = [0.0f64; 2];

    // -------------------------------------------------------------------
    // Main computation loop (Fortran label 30).
    // -------------------------------------------------------------------
    'main_loop: loop {
        let nn = 2.min(nd);

        let mut i = 1;
        while i <= nn {
            fn_val = fnu + ((nd - i) as f64);
            init = 0;
            zunik(
                zr, zi, fn_val, 1, 0, tol, &mut init, &mut phir, &mut phii, &mut zeta1r,
                &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut sumr, &mut sumi, &mut cwrkr,
                &mut cwrki,
            );

            if kode != 1 {
                let mut str = zr + zeta2r;
                let mut sti = zi + zeta2i;
                let rast = fn_val / zabs(str, sti);
                str = str * rast * rast;
                sti = -sti * rast * rast;
                s1r = -zeta1r + str;
                s1i = -zeta1i + sti + zi;
            } else {
                s1r = -zeta1r + zeta2r;
                s1i = -zeta1i + zeta2i;
            }

            // -------------------------------------------------------------
            // Test for underflow and overflow
            // -------------------------------------------------------------
            rs1 = s1r;
            if rs1.abs() > elim {
                // Underflow/overflow handling (Fortran label 110).
                if let Some(control) = underflow_110(rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu, kode, tol, elim, alim) {
                    match control {
                        Flow::Return => return,
                        Flow::Continue => continue 'main_loop,
                        Flow::SetNlastAndReturn => return,
                    }
                }
                // (underflow_110 always returns Some; the None branch is
                // unreachable, but kept for clarity.)
                unreachable!();
            }

            if i == 1 {
                iflag = 2;
            }

            if rs1.abs() < alim {
                // On scale — skip the refine test.
            } else {
                // ---------------------------------------------------------
                // Refine test and scale
                // ---------------------------------------------------------
                let aphi = zabs(phir, phii);
                rs1 += aphi.ln();
                if rs1.abs() > elim {
                    if let Some(control) = underflow_110(rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu, kode, tol, elim, alim) {
                        match control {
                            Flow::Return => return,
                            Flow::Continue => continue 'main_loop,
                            Flow::SetNlastAndReturn => return,
                        }
                    }
                    unreachable!();
                }
                if i == 1 {
                    iflag = 1;
                }
                if rs1 >= 0.0 {
                    if i == 1 {
                        iflag = 3;
                    }
                }
            }

            // -------------------------------------------------------------
            // Scale S1 if |S1| < ASCLE
            // -------------------------------------------------------------
            let mut s2r = phir * sumr - phii * sumi;
            let mut s2i = phir * sumi + phii * sumr;
            let str = rs1.exp() * cssr[(iflag - 1) as usize];
            s1r = str * s1i.cos();
            s1i = str * s1i.sin();
            let str2 = s2r * s1r - s2i * s1i;
            s2i = s2r * s1i + s2i * s1r;
            s2r = str2;

            if iflag == 1 {
                let mut nw: i32 = 0;
                zuchk(s2r, s2i, &mut nw, bry[0], tol);
                if nw != 0 {
                    if let Some(control) = underflow_110(rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu, kode, tol, elim, alim) {
                        match control {
                            Flow::Return => return,
                            Flow::Continue => continue 'main_loop,
                            Flow::SetNlastAndReturn => return,
                        }
                    }
                    unreachable!();
                }
            }

            cyr[(i - 1) as usize] = s2r;
            cyi[(i - 1) as usize] = s2i;
            let m = nd - i + 1;
            yr[(m - 1) as usize] = s2r * csrr[(iflag - 1) as usize];
            yi[(m - 1) as usize] = s2i * csrr[(iflag - 1) as usize];

            i += 1;
        }

        if nd <= 2 {
            return;
        }

        // -----------------------------------------------------------------
        // Forward recurrence with scaling.
        // -----------------------------------------------------------------
        let rast = 1.0 / zabs(zr, zi);
        let str = zr * rast;
        let sti = -zi * rast;
        let rzr = (str + str) * rast;
        let rzi = (sti + sti) * rast;
        bry[1] = 1.0 / bry[0];
        bry[2] = d1mach(2);

        let mut s1r = cyr[0];
        let mut s1i = cyi[0];
        let mut s2r = cyr[1];
        let mut s2i = cyi[1];
        let mut c1r = csrr[(iflag - 1) as usize];
        let mut ascle = bry[(iflag - 1) as usize];
        let mut k = nd - 2;
        fn_val = k as f64;

        let mut i = 3;
        while i <= nd {
            let c2r_in = s2r;
            let c2i_in = s2i;
            s2r = s1r + (fnu + fn_val) * (rzr * c2r_in - rzi * c2i_in);
            s2i = s1i + (fnu + fn_val) * (rzr * c2i_in + rzi * c2r_in);
            s1r = c2r_in;
            s1i = c2i_in;
            let c2r = s2r * c1r;
            let c2i = s2i * c1r;
            yr[(k - 1) as usize] = c2r;
            yi[(k - 1) as usize] = c2i;
            k -= 1;
            fn_val -= 1.0;

            if iflag < 3 {
                let c2m = c2r.abs().max(c2i.abs());
                if c2m > ascle {
                    iflag += 1;
                    ascle = bry[(iflag - 1) as usize];
                    s1r *= c1r;
                    s1i *= c1r;
                    s2r = c2r;
                    s2i = c2i;
                    s1r *= cssr[(iflag - 1) as usize];
                    s1i *= cssr[(iflag - 1) as usize];
                    s2r *= cssr[(iflag - 1) as usize];
                    s2i *= cssr[(iflag - 1) as usize];
                    c1r = csrr[(iflag - 1) as usize];
                }
            }
            i += 1;
        }

        return;
    }
}

// ---------------------------------------------------------------------
// Enum for the underflow/overflow continuation.
// ---------------------------------------------------------------------
enum Flow {
    Return,
    Continue,
    SetNlastAndReturn,
}

/// Fortran label-110 underflow/overflow handling.
///
/// Returns `Some(Flow::...)` indicating what the caller should do:
/// - `Flow::Return`: return from `zuni1` immediately.
/// - `Flow::Continue`: `continue 'main_loop`.
/// - `Flow::SetNlastAndReturn`: set `nlast = nd` and return.
#[allow(clippy::too_many_arguments)]
fn underflow_110(
    rs1: f64,
    nd: &mut i32,
    nz: &mut i32,
    fn_val: &mut f64,
    fnul: f64,
    nlast: &mut i32,
    yr: &mut [f64],
    yi: &mut [f64],
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    tol: f64,
    elim: f64,
    alim: f64,
) -> Option<Flow> {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;

    if rs1 > 0.0 {
        *nz = -1;
        return Some(Flow::Return);
    }
    yr[(*nd - 1) as usize] = ZEROR;
    yi[(*nd - 1) as usize] = ZEROI;
    *nz += 1;
    *nd -= 1;
    if *nd == 0 {
        return Some(Flow::Return);
    }

    let mut nuf: i32 = 0;
    zuoik(zr, zi, fnu, kode, 1, *nd, yr, yi, &mut nuf, tol, elim, alim);
    if nuf < 0 {
        *nz = -1;
        return Some(Flow::Return);
    }
    *nd -= nuf;
    *nz += nuf;
    if *nd == 0 {
        return Some(Flow::Return);
    }

    *fn_val = fnu + ((*nd - 1) as f64);
    if *fn_val >= fnul {
        return Some(Flow::Continue);
    }
    *nlast = *nd;
    Some(Flow::SetNlastAndReturn)
}