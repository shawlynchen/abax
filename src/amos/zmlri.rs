use crate::amos::dgamln::dgamln;
use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zexp::zexp;
use crate::amos::zlog::zlog;
use crate::amos::zmlt::zmlt;

pub(crate) fn zmlri(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    tol: f64,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const CONER: f64 = 1.0;
    const CONEI: f64 = 0.0;

    let n_usize = n as usize;
    debug_assert!(n >= 0);
    debug_assert!(yr.len() >= n_usize, "yr must have at least n elements");
    debug_assert!(yi.len() >= n_usize, "yi must have at least n elements");

    let scle = d1mach(1) / tol;
    *nz = 0;

    let az = zabs(zr, zi);
    let iaz = az as i32;
    let ifnu = fnu as i32;
    let inu = ifnu + n - 1;

    let mut at = (iaz as f64) + 1.0;
    let raz = 1.0 / az;
    let (str, sti) = (zr * raz, -zi * raz);
    let mut ckr = str * at * raz;
    let mut cki = sti * at * raz;
    let rzr = (str + str) * raz;
    let rzi = (sti + sti) * raz;

    let (mut p1r, mut p1i) = (ZEROR, ZEROI);
    let (mut p2r, mut p2i) = (CONER, CONEI);

    let ack = (at + 1.0) * raz;
    let rho = ack + (ack * ack - 1.0).sqrt();
    let rho2 = rho * rho;
    let mut tst = (rho2 + rho2) / ((rho2 - 1.0) * (rho - 1.0));
    tst /= tol;

    // Compute relative truncation error index for series
    let mut ak = at;
    let mut i: i32 = 81; // sentinel: 81 == loop exhausted, matching Fortran
    for idx in 1..=80 {
        let (ptr, pti) = (p2r, p2i);
        p2r = p1r - (ckr * ptr - cki * pti);
        p2i = p1i - (cki * ptr + ckr * pti);
        p1r = ptr;
        p1i = pti;
        ckr += rzr;
        cki += rzi;
        let ap = zabs(p2r, p2i);
        if ap > tst * ak * ak {
            i = idx;
            break;
        }
        ak += 1.0;
    }

    if i == 81 {
        *nz = -2;
        return;
    }

    i += 1;
    let mut k: i32 = 0;

    // Compute relative truncation error for ratios (only when INU >= IAZ)
    if inu >= iaz {
        p1r = ZEROR;
        p1i = ZEROI;
        p2r = CONER;
        p2i = CONEI;
        at = (inu as f64) + 1.0;
        let (str2, sti2) = (zr * raz, -zi * raz);
        ckr = str2 * at * raz;
        cki = sti2 * at * raz;
        let ack2 = at * raz;
        tst = (ack2 / tol).sqrt();

        let mut itime = 1;
        k = 81; // sentinel
        for idx in 1..=80 {
            let (ptr, pti) = (p2r, p2i);
            p2r = p1r - (ckr * ptr - cki * pti);
            p2i = p1i - (ckr * pti + cki * ptr);
            p1r = ptr;
            p1i = pti;
            ckr += rzr;
            cki += rzi;
            let ap = zabs(p2r, p2i);
            if ap < tst {
                continue;
            }
            if itime == 2 {
                k = idx;
                break;
            }
            let ack = zabs(ckr, cki);
            let flam = ack + (ack * ack - 1.0).sqrt();
            let fkap = ap / zabs(p1r, p1i);
            let rho = flam.min(fkap);
            tst = tst * (rho / (rho * rho - 1.0)).sqrt();
            itime = 2;
        }

        if k == 81 {
            *nz = -2;
            return;
        }
    }

    // Backward recurrence and sum normalizing relation
    k += 1;
    let kk = (i + iaz).max(k + inu);
    let mut fkk = kk as f64;
    p1r = ZEROR;
    p1i = ZEROI;

    // Scale P2 and sum by SCLE
    p2r = scle;
    p2i = ZEROI;
    let fnf = fnu - (ifnu as f64);
    let tfnf = fnf + fnf;

    let mut idum: i32 = 0;
    let mut bk = dgamln(fkk + tfnf + 1.0, &mut idum)
        - dgamln(fkk + 1.0, &mut idum)
        - dgamln(tfnf + 1.0, &mut idum);
    bk = bk.exp();

    let (mut sumr, mut sumi) = (ZEROR, ZEROI);
    let km = kk - inu;

    for _ in 1..=km {
        let (ptr, pti) = (p2r, p2i);
        p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
        p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
        p1r = ptr;
        p1i = pti;
        let ak = 1.0 - tfnf / (fkk + tfnf);
        let ack = bk * ak;
        sumr += (ack + bk) * p1r;
        sumi += (ack + bk) * p1i;
        bk = ack;
        fkk -= 1.0;
    }

    yr[(n - 1) as usize] = p2r;
    yi[(n - 1) as usize] = p2i;

    if n != 1 {
        for idx in 2..=n {
            let (ptr, pti) = (p2r, p2i);
            p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
            p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
            p1r = ptr;
            p1i = pti;
            let ak = 1.0 - tfnf / (fkk + tfnf);
            let ack = bk * ak;
            sumr += (ack + bk) * p1r;
            sumi += (ack + bk) * p1i;
            bk = ack;
            fkk -= 1.0;
            let m = n - idx + 1;
            yr[(m - 1) as usize] = p2r;
            yi[(m - 1) as usize] = p2i;
        }
    }

    if ifnu > 0 {
        for _ in 1..=ifnu {
            let (ptr, pti) = (p2r, p2i);
            p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
            p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
            p1r = ptr;
            p1i = pti;
            let ak = 1.0 - tfnf / (fkk + tfnf);
            let ack = bk * ak;
            sumr += (ack + bk) * p1r;
            sumi += (ack + bk) * p1i;
            bk = ack;
            fkk -= 1.0;
        }
    }

    // Normalize the result
    let (mut ptr, pti) = (zr, zi);
    if kode == 2 {
        ptr = ZEROR;
    }

    let (mut str_out, mut sti_out) = (0.0f64, 0.0f64);
    let _ = zlog(rzr, rzi, &mut str_out, &mut sti_out);
    let p1r = -fnf * str_out + ptr;
    let p1i = -fnf * sti_out + pti;

    let ap = dgamln(1.0 + fnf, &mut idum);
    let ptr = p1r - ap;
    let pti = p1i;

    // The division cexp(PT)/(SUM+P2) is altered to avoid overflow
    // in the denominator by squaring large quantities.
    p2r += sumr;
    p2i += sumi;
    let ap = zabs(p2r, p2i);
    let p1r_inv = 1.0 / ap;

    let (mut str_e, mut sti_e) = (0.0f64, 0.0f64);
    zexp(ptr, pti, &mut str_e, &mut sti_e);

    let ckr = str_e * p1r_inv;
    let cki = sti_e * p1r_inv;
    let ptr_inv = p2r * p1r_inv;
    let pti_inv = -p2i * p1r_inv;

    let (mut cnormr, mut cnormi) = (0.0f64, 0.0f64);
    zmlt(ckr, cki, ptr_inv, pti_inv, &mut cnormr, &mut cnormi);

    for i in 0..n_usize {
        let str_final = yr[i] * cnormr - yi[i] * cnormi;
        yi[i] = yr[i] * cnormi + yi[i] * cnormr;
        yr[i] = str_final;
    }
}