use crate::amos::zabs::zabs;
use crate::amos::zdiv::zdiv;

/// ZRATI computes ratios of I Bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, Vol 77B, p111-114, September, 1973,
/// Bessel Functions I and J of Complex Argument and Integer Order,
/// by D. J. Sookne.
///
/// * `zr, zi`    – real and imaginary part of Z
/// * `fnu`       – order
/// * `n`         – number of terms
/// * `cyr, cyi`  – real/imag arrays (output)
/// * `tol`       – tolerance
pub(crate) fn zrati(
    zr: f64,
    zi: f64,
    fnu: f64,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    tol: f64,
) {
    const CZEROR: f64 = 0.0;
    const CZEROI: f64 = 0.0;
    const CONER: f64 = 1.0;
    const CONEI: f64 = 0.0;
    const RT2: f64 = 1.41421356237309505;

    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");

    let az = zabs(zr, zi);
    let inu = fnu as i32;
    let idnu = inu + n - 1;
    let magz = az as i32;
    let amagz = (magz + 1) as f64;
    let fdnu = idnu as f64;
    let fnup = amagz.max(fdnu);
    let mut id = idnu - magz - 1;

    let mut itime = 1;
    let mut k: i32 = 1;

    let ptr = 1.0 / az;
    let rzr = ptr * (zr + zr) * ptr;
    let rzi = -ptr * (zi + zi) * ptr;

    let mut t1r = rzr * fnup;
    let mut t1i = rzi * fnup;
    let mut p2r = -t1r;
    let mut p2i = -t1i;
    let mut p1r = CONER;
    let mut p1i = CONEI;
    t1r += rzr;
    t1i += rzi;

    if id > 0 {
        id = 0;
    }

    let mut ap2 = zabs(p2r, p2i);
    let mut ap1 = zabs(p1r, p1i);

    // -------------------------------------------------------------------
    // The overflow test on K(FNU+I-1,Z) before the call to CBKNU
    // guarantees that P2 is on scale. Scale TEST1 and all subsequent
    // P2 values by AP1 to ensure that an overflow does not occur
    // prematurely.
    // -------------------------------------------------------------------
    let arg = (ap2 + ap2) / (ap1 * tol);
    let test1 = arg.sqrt();
    let mut test = test1;
    let rap1 = 1.0 / ap1;
    p1r *= rap1;
    p1i *= rap1;
    p2r *= rap1;
    p2i *= rap1;
    ap2 *= rap1;

    // Forward recurrence to determine the starting index.
    // Fortran: DO 10 (infinite loop with GO TO 20 to break).
    loop {
        k += 1;
        ap1 = ap2;
        let (ptr, pti) = (p2r, p2i);
        p2r = p1r - (t1r * ptr - t1i * pti);
        p2i = p1i - (t1r * pti + t1i * ptr);
        p1r = ptr;
        p1i = pti;
        t1r += rzr;
        t1i += rzi;
        ap2 = zabs(p2r, p2i);

        if ap1 <= test {
            continue;
        }
        if itime == 2 {
            break;
        }

        let ak = zabs(t1r, t1i) * 0.5;
        let flam = ak + (ak * ak - 1.0).sqrt();
        let rho = (ap2 / ap1).min(flam);
        test = test1 * (rho / (rho * rho - 1.0)).sqrt();
        itime = 2;
    }

    // Backward recurrence.
    let kk = k + 1 - id;
    let ak = kk as f64;
    let mut t1r = ak;
    let _t1i = CZEROI;
    let dfnu = fnu + ((n - 1) as f64);

    p1r = 1.0 / ap2;
    p1i = CZEROI;
    p2r = CZEROR;
    p2i = CZEROI;

    for _ in 1..=kk {
        let (ptr, pti) = (p1r, p1i);
        let rap1 = dfnu + t1r;
        let ttr = rzr * rap1;
        let tti = rzi * rap1;
        p1r = (ptr * ttr - pti * tti) + p2r;
        p1i = (ptr * tti + pti * ttr) + p2i;
        p2r = ptr;
        p2i = pti;
        t1r -= CONER;
    }

    if p1r == CZEROR && p1i == CZEROI {
        p1r = tol;
        p1i = tol;
    }

    let (mut cyr_n, mut cyi_n) = (0.0f64, 0.0f64);
    zdiv(p2r, p2i, p1r, p1i, &mut cyr_n, &mut cyi_n);
    cyr[(n - 1) as usize] = cyr_n;
    cyi[(n - 1) as usize] = cyi_n;

    if n == 1 {
        return;
    }

    let mut k = n - 1;
    let ak = k as f64;
    let mut t1r = ak;
    let t1i = CZEROI;
    let cdfnur = fnu * rzr;
    let cdfnui = fnu * rzi;

    // Fortran: DO 60 I=2,N, with K decrementing.
    for _ in 2..=n {
        let mut ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[(k + 1) as usize];
        let mut pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[(k + 1) as usize];
        let mut ak = zabs(ptr, pti);

        if ak == CZEROR {
            ptr = tol;
            pti = tol;
            ak = tol * RT2;
        }

        let rak = CONER / ak;
        cyr[k as usize] = rak * ptr * rak;
        cyi[k as usize] = -rak * pti * rak;
        t1r -= CONER;
        k -= 1;
    }
}