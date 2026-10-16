use crate::amos::dgamln::dgamln;
use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zdiv::zdiv;
use crate::amos::zlog::zlog;
use crate::amos::zmlt::zmlt;
use crate::amos::zuchk::zuchk;

/// ZSERI computes the I Bessel function for Re(Z) >= 0 by
/// means of the power series for large |Z| in the region
/// |Z| <= 2*sqrt(FNU+1). NZ=0 is a normal return.
/// NZ>0 means that the last NZ components were set to zero
/// due to underflow. NZ<0 means underflow occurred, but the
/// condition |Z| <= 2*sqrt(FNU+1) was violated and the
/// computation must be completed in another routine with
/// N = N - abs(NZ).
#[allow(clippy::too_many_arguments)]
pub(crate) fn zseri(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const CONER: f64 = 1.0;
    const CONEI: f64 = 0.0;

    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;

    let az = zabs(zr, zi);

    // AZ == 0 branch (Fortran label 160).
    if az == 0.0 {
        yr[0] = ZEROR;
        yi[0] = ZEROI;
        if fnu == 0.0 {
            yr[0] = CONER;
            yi[0] = CONEI;
        }
        if n == 1 {
            return;
        }
        for i in 1..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        return;
    }

    let arm = 1.0e+3 * d1mach(1);
    let rtr1 = arm.sqrt();
    let mut crscr = 1.0;
    let mut iflag = 0;

    // AZ < ARM branch (Fortran label 150).
    if az < arm {
        *nz = n;
        if fnu == 0.0 {
            *nz -= 1;
        }
        yr[0] = ZEROR;
        yi[0] = ZEROI;
        if fnu == 0.0 {
            yr[0] = CONER;
            yi[0] = CONEI;
        }
        if n == 1 {
            return;
        }
        for i in 1..n_usize {
            yr[i] = ZEROR;
            yi[i] = ZEROI;
        }
        return;
    }

    let hzr = 0.5 * zr;
    let hzi = 0.5 * zi;
    let (mut czr, mut czi) = (ZEROR, ZEROI);
    if az > rtr1 {
        zmlt(hzr, hzi, hzr, hzi, &mut czr, &mut czi);
    }

    let acz = zabs(czr, czi);
    let mut nn = n;

    let (mut ckr, mut cki) = (0.0f64, 0.0f64);
    let mut idum: i32 = 0;
    let _ = zlog(hzr, hzi, &mut ckr, &mut cki);

    // Work arrays for the scaled recurrence.
    let mut wr = [0.0f64; 2];
    let mut wi = [0.0f64; 2];

    // -------------------------------------------------------------------
    // Main loop over decreasing NN.
    // -------------------------------------------------------------------
    'main_loop: loop {
        let mut dfnu = fnu + ((nn - 1) as f64);
        let mut fnup = dfnu + 1.0;

        // Underflow test
        let mut ak1r = ckr * dfnu;
        let ak1i = cki * dfnu;
        let ak = dgamln(fnup, &mut idum);
        ak1r -= ak;
        if kode == 2 {
            ak1r -= zr;
        }

        if ak1r <= -elim {
            // Underflow branch (Fortran label 30).
            *nz += 1;
            yr[(nn - 1) as usize] = ZEROR;
            yi[(nn - 1) as usize] = ZEROI;

            if acz > dfnu {
                *nz = -*nz;
                return;
            }
            nn -= 1;
            if nn == 0 {
                return;
            }
            continue 'main_loop;
        }

        // Fortran label 40: continuation.
        let mut ss = 0.0;
        let mut ascle = 0.0;
        if ak1r <= -alim {
            iflag = 1;
            ss = 1.0 / tol;
            crscr = tol;
            ascle = arm * ss;
        }

        let mut aa = ak1r.exp();
        if iflag == 1 {
            aa *= ss;
        }
        let mut coefr = aa * ak1i.cos();
        let mut coefi = aa * ak1i.sin();

        let atol = tol * acz / fnup;
        let il = nn.min(2);

        for i in 1..=il {
            dfnu = fnu + ((nn - i) as f64);
            fnup = dfnu + 1.0;
            let (mut s1r, mut s1i) = (CONER, CONEI);

            if acz >= tol * fnup {
                // Series accumulation (Fortran label 60).
                let mut ak1r = CONER;
                let mut ak1i = CONEI;
                let mut ak = fnup + 2.0;
                let mut s = fnup;
                let mut aa = 2.0;

                loop {
                    let rs = 1.0 / s;
                    let str = ak1r * czr - ak1i * czi;
                    let sti = ak1r * czi + ak1i * czr;
                    ak1r = str * rs;
                    ak1i = sti * rs;
                    s1r += ak1r;
                    s1i += ak1i;
                    s += ak;
                    ak += 2.0;
                    aa = aa * acz * rs;
                    if aa <= atol {
                        break;
                    }
                }
            }

            let s2r = s1r * coefr - s1i * coefi;
            let s2i = s1r * coefi + s1i * coefr;
            wr[(i - 1) as usize] = s2r;
            wi[(i - 1) as usize] = s2i;

            if iflag != 0 {
                let mut nw: i32 = 0;
                zuchk(s2r, s2i, &mut nw, ascle, tol);
                if nw != 0 {
                    // Underflow branch (Fortran label 30).
                    *nz += 1;
                    yr[(nn - 1) as usize] = ZEROR;
                    yi[(nn - 1) as usize] = ZEROI;

                    if acz > dfnu {
                        *nz = -*nz;
                        return;
                    }
                    nn -= 1;
                    if nn == 0 {
                        return;
                    }
                    continue 'main_loop;
                }
            }

            let m = nn - i + 1;
            yr[(m - 1) as usize] = s2r * crscr;
            yi[(m - 1) as usize] = s2i * crscr;

            if i != il {
                let (mut str, mut sti) = (0.0f64, 0.0f64);
                zdiv(coefr, coefi, hzr, hzi, &mut str, &mut sti);
                coefr = str * dfnu;
                coefi = sti * dfnu;
            }
        }

        if nn <= 2 {
            return;
        }

        let mut k = nn - 2;
        let mut ak = k as f64;
        let raz = 1.0 / az;
        let str = zr * raz;
        let sti = -zi * raz;
        let rzr = (str + str) * raz;
        let rzi = (sti + sti) * raz;

        if iflag != 1 {
            // Unscaled backward recurrence (Fortran label 100).
            let ib = 3;
            for _i in ib..=nn {
                let kk = k as usize;
                yr[kk - 1] = (ak + fnu) * (rzr * yr[kk] - rzi * yi[kk]) + yr[kk + 1];
                yi[kk - 1] = (ak + fnu) * (rzr * yi[kk] + rzi * yr[kk]) + yi[kk + 1];
                ak -= 1.0;
                k -= 1;
            }
            return;
        }

        // -----------------------------------------------------------------
        // Recur backward with scaled values (Fortran label 120).
        //
        // exp(-ALIM) = exp(-ELIM)/TOL = approx. one precision above the
        // underflow limit = ASCLE = D1MACH(1)*SS*1.0D+3
        // -----------------------------------------------------------------
        let mut s1r = wr[0];
        let mut s1i = wi[0];
        let mut s2r = wr[1];
        let mut s2i = wi[1];

        let mut ib = nn + 1; // default: no continuation
        for l in 3..=nn {
            let ckr = s2r;
            let cki = s2i;
            s2r = s1r + (ak + fnu) * (rzr * ckr - rzi * cki);
            s2i = s1i + (ak + fnu) * (rzr * cki + rzi * ckr);
            s1r = ckr;
            s1i = cki;
            let ckr = s2r * crscr;
            let cki = s2i * crscr;
            yr[(k - 1) as usize] = ckr;
            yi[(k - 1) as usize] = cki;
            ak -= 1.0;
            k -= 1;

            if zabs(ckr, cki) > ascle {
                ib = l + 1;
                break;
            }
        }

        if ib > nn {
            return;
        }

        // Continue the unscaled recurrence from IB..NN (Fortran label 100).
        for _i in ib..=nn {
            let kk = k as usize;
            yr[kk - 1] = (ak + fnu) * (rzr * yr[kk] - rzi * yi[kk]) + yr[kk + 1];
            yi[kk - 1] = (ak + fnu) * (rzr * yi[kk] + rzi * yr[kk]) + yi[kk + 1];
            ak -= 1.0;
            k -= 1;
        }
        return;
    }
}