use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zuni1::zuni1;
use crate::amos::zuni2::zuni2;

/// ZBUNI computes the I Bessel function for large |Z| > FNUL and
/// FNU+N-1 < FNUL. The order is increased from FNU+N-1 greater
/// than FNUL by adding NUI and computing according to the uniform
/// asymptotic expansion for I(FNU,Z) on IFORM=1 and the expansion
/// for J(FNU,Z) on IFORM=2.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbuni(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    nui: i32,
    nlast: &mut i32,
    fnul: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;

    let ax = zr.abs() * 1.7321;
    let ay = zi.abs();
    let iform = if ay > ax { 2 } else { 1 };

    if nui != 0 {
        // -------------------------------------------------------------
        // Order-increased path.
        // -------------------------------------------------------------
        let mut fnui = nui as f64;
        let dfnu = fnu + ((n - 1) as f64);
        let gnu = dfnu + fnui;

        let mut cyr = [0.0f64; 2];
        let mut cyi = [0.0f64; 2];
        let mut nw: i32 = 0;

        if iform != 2 {
            zuni1(
                zr, zi, gnu, kode, 2, &mut cyr, &mut cyi, &mut nw, nlast, fnul, tol,
                elim, alim,
            );
        } else {
            zuni2(
                zr, zi, gnu, kode, 2, &mut cyr, &mut cyi, &mut nw, nlast, fnul, tol,
                elim, alim,
            );
        }

        if nw < 0 {
            *nz = if nw == -2 { -2 } else { -1 };
            return;
        }
        if nw != 0 {
            *nlast = n;
            return;
        }

        let str = zabs(cyr[0], cyi[0]);

        // -------------------------------------------------------------
        // Scale backward recurrence. BRY(3) is defined but never used.
        // -------------------------------------------------------------
        let mut bry = [0.0f64; 3];
        bry[0] = 1.0e+3 * d1mach(1) / tol;
        bry[1] = 1.0 / bry[0];
        bry[2] = bry[1];

        let mut iflag: i32 = 2;
        let mut ascle = bry[1];
        let mut csclr = 1.0;

        if str <= bry[0] {
            iflag = 1;
            ascle = bry[0];
            csclr = 1.0 / tol;
        } else if str >= bry[1] {
            iflag = 3;
            ascle = bry[2];
            csclr = tol;
        }

        let mut cscrr = 1.0 / csclr;
        let mut s1r = cyr[1] * csclr;
        let mut s1i = cyi[1] * csclr;
        let mut s2r = cyr[0] * csclr;
        let mut s2i = cyi[0] * csclr;
        let raz = 1.0 / zabs(zr, zi);
        let str = zr * raz;
        let sti = -zi * raz;
        let rzr = (str + str) * raz;
        let rzi = (sti + sti) * raz;

        for _i in 1..=nui {
            let str = s2r;
            let sti = s2i;
            s2r = (dfnu + fnui) * (rzr * str - rzi * sti) + s1r;
            s2i = (dfnu + fnui) * (rzr * sti + rzi * str) + s1i;
            s1r = str;
            s1i = sti;
            fnui -= 1.0;

            if iflag < 3 {
                let str = s2r * cscrr;
                let sti = s2i * cscrr;
                let c1r = str.abs();
                let c1i = sti.abs();
                let c1m = c1r.max(c1i);
                if c1m > ascle {
                    iflag += 1;
                    ascle = bry[(iflag - 1) as usize];
                    s1r *= cscrr;
                    s1i *= cscrr;
                    s2r = str;
                    s2i = sti;
                    csclr *= tol;
                    cscrr = 1.0 / csclr;
                    s1r *= csclr;
                    s1i *= csclr;
                    s2r *= csclr;
                    s2i *= csclr;
                }
            }
        }

        yr[(n - 1) as usize] = s2r * cscrr;
        yi[(n - 1) as usize] = s2i * cscrr;
        if n == 1 {
            return;
        }

        let nl = n - 1;
        fnui = nl as f64;
        let mut k = nl;

        for _i in 1..=nl {
            let str = s2r;
            let sti = s2i;
            s2r = (fnu + fnui) * (rzr * str - rzi * sti) + s1r;
            s2i = (fnu + fnui) * (rzr * sti + rzi * str) + s1i;
            s1r = str;
            s1i = sti;
            let str = s2r * cscrr;
            let sti = s2i * cscrr;
            yr[(k - 1) as usize] = str;
            yi[(k - 1) as usize] = sti;
            fnui -= 1.0;
            k -= 1;

            if iflag < 3 {
                let c1r = str.abs();
                let c1i = sti.abs();
                let c1m = c1r.max(c1i);
                if c1m > ascle {
                    iflag += 1;
                    ascle = bry[(iflag - 1) as usize];
                    s1r *= cscrr;
                    s1i *= cscrr;
                    s2r = str;
                    s2i = sti;
                    csclr *= tol;
                    cscrr = 1.0 / csclr;
                    s1r *= csclr;
                    s1i *= csclr;
                    s2r *= csclr;
                    s2i *= csclr;
                }
            }
        }
        return;
    }

    // -------------------------------------------------------------------
    // NUI = 0 path.
    // -------------------------------------------------------------------
    let mut nw: i32 = 0;
    if iform != 2 {
        zuni1(zr, zi, fnu, kode, n, yr, yi, &mut nw, nlast, fnul, tol, elim, alim);
    } else {
        zuni2(zr, zi, fnu, kode, n, yr, yi, &mut nw, nlast, fnul, tol, elim, alim);
    }

    if nw < 0 {
        *nz = if nw == -2 { -2 } else { -1 };
        return;
    }
    *nz = nw;
}