use crate::amos::{zabs, zlog, zuchk};

/// Set K functions to zero on underflow, continue recurrence
/// on scaled functions until two members come on scale, then
/// return with min(nz+2, n) values scaled by 1/tol.
///
/// * `zrr, zri`  – real and imaginary part of the scale factor ZR
/// * `fnu`       – order
/// * `n`         – number of terms
/// * `yr, yi`    – real/imag arrays; on input they hold the starting
///                 values, on output the scaled results
/// * `nz`        – number of underflows (output)
/// * `rzr, rzi`  – real/imag part of the recurrence multiplier
/// * `ascle`     – underflow / overflow threshold
/// * `tol`       – scaling tolerance
/// * `elim`      – exponent limit
pub(crate) fn zkscl(
    zrr: f64,
    zri: f64,
    fnu: f64,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    rzr: f64,
    rzi: f64,
    ascle: f64,
    tol: f64,
    elim: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;

    let n_usize = n as usize;

    debug_assert!(n >= 0);
    debug_assert!(yr.len() >= n_usize);
    debug_assert!(yi.len() >= n_usize);

    *nz = 0;
    let mut ic: i32 = 0;
    let nn = i32::min(2, n);

    // cyr / cyi hold copies of the first min(2, n) values.
    let mut cyr = [0.0f64; 2];
    let mut cyi = [0.0f64; 2];

    // ---- First loop: handle the first min(2, n) entries. -------------
    // Fortran: DO 10 I = 1, NN  (1-based index, array index = I)
    for i in 0..(nn as usize) {
        let s1r = yr[i];
        let s1i = yi[i];
        cyr[i] = s1r;
        cyi[i] = s1i;

        let a = zabs(s1r, s1i);
        let acs = -zrr + a.ln();
        *nz += 1;
        yr[i] = ZEROR;
        yi[i] = ZEROI;

        if acs < -elim {
            continue;
        }

        let mut csr = 0.0f64;
        let mut csi = 0.0f64;
        let _idum = zlog(s1r, s1i, &mut csr, &mut csi);
        csr -= zrr;
        csi -= zri;

        let str = csr.exp() / tol;
        csr = str * csi.cos();
        csi = str * csi.sin();

        let mut nw: i32 = 0;
        zuchk(csr, csi, &mut nw, ascle, tol);
        if nw != 0 {
            continue;
        }

        yr[i] = csr;
        yi[i] = csi;
        ic = (i + 1) as i32; // Fortran IC = I (1-based)
        *nz -= 1;
    }

    if n == 1 {
        return;
    }

    // ---- Handle the "IC > 1" branch. ---------------------------------
    // Fortran label 20: IF (IC.GT.1) GO TO 20 ; then zero out Y(1).
    if ic <= 1 {
        yr[0] = ZEROR;
        yi[0] = ZEROI;
        *nz = 2;
    }

    if n == 2 {
        return;
    }
    if *nz == 0 {
        return;
    }

    // ---- Recurrence setup. -------------------------------------------
    let fn_val = fnu + 1.0;
    let mut ckr = fn_val * rzr;
    let mut cki = fn_val * rzi;

    let mut s1r = cyr[0];
    let mut s1i = cyi[0];
    let mut s2r = cyr[1];
    let mut s2i = cyi[1];

    let helim = 0.5 * elim;
    let elm = (-elim).exp();
    let celmr = elm;
    let mut zdr = zrr;
    let zdi = zri;

    let mut kk: i32;
    let mut broke_early = false;

    // ---- Main recurrence loop: find two consecutive values on scale. --
    // Fortran: DO 30 I = 3, N
    for i in 2..(n as usize) {
        kk = (i + 1) as i32; // Fortran KK = I

        let csr = s2r;
        let csi = s2i;

        s2r = ckr * csr - cki * csi + s1r;
        s2i = cki * csr + ckr * csi + s1i;

        s1r = csr;
        s1i = csi;

        ckr += rzr;
        cki += rzi;

        let a = zabs(s2r, s2i);
        let alas = a.ln();
        let acs = -zdr + alas;
        *nz += 1;

        yr[i] = ZEROR;
        yi[i] = ZEROI;

        // Fortran label 25 check: IF (ACS.LT.(-ELIM)) GO TO 25
        // Otherwise do the scaling/log/zuchk block.
        if acs >= -elim {
            let mut csr2 = 0.0f64;
            let mut csi2 = 0.0f64;
            let _ = zlog(s2r, s2i, &mut csr2, &mut csi2);
            csr2 -= zdr;
            csi2 -= zdi;

            let str = csr2.exp() / tol;
            csr2 = str * csi2.cos();
            csi2 = str * csi2.sin();

            let mut nw: i32 = 0;
            zuchk(csr2, csi2, &mut nw, ascle, tol);
            if nw == 0 {
                yr[i] = csr2;
                yi[i] = csi2;
                *nz -= 1;

                // Fortran: IF (IC.EQ.KK-1) GO TO 40  =>  NZ = KK - 2, done.
                if ic == kk - 1 {
                    *nz = kk - 2;
                    broke_early = true;
                    break;
                }
                ic = kk;
                continue;
            }
        }

        // Fortran label 25 body.
        if alas < helim {
            continue;
        }

        zdr -= elim;
        s1r *= celmr;
        s1i *= celmr;
        s2r *= celmr;
        s2i *= celmr;
    }

    // ---- Final NZ adjustment (Fortran label 45 logic). ---------------
    // Fortran: NZ = N; IF (IC.EQ.N) NZ = N-1
    // This is only executed when we did NOT jump to label 40.
    if !broke_early {
        *nz = n;
        if ic == n {
            *nz = n - 1;
        }
    }

    // ---- Zero out the first NZ entries. ------------------------------
    // Fortran: DO 50 I = 1, NZ
    let nz_final = (*nz).max(0) as usize;
    for i in 0..nz_final {
        yr[i] = ZEROR;
        yi[i] = ZEROI;
    }
}