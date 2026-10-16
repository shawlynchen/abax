use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zairy::zairy;
use crate::amos::zuchk::zuchk;
use crate::amos::zunhj::zunhj;
use crate::amos::zuoik::zuoik;

const ZEROR: f64 = 0.0;
const ZEROI: f64 = 0.0;
const CONER: f64 = 1.0;
const HPI: f64 = 1.57079632679489662;
const AIC: f64 = 1.265512123484645396;

const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
const CIPI: [f64; 4] = [0.0, 1.0, 0.0, -1.0];

enum Flow {
    Return,
    Continue,
    SetNlastAndReturn,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn zuni2(
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
    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    *nz = 0;
    let mut nd = n;
    *nlast = 0;

    let cscl = 1.0 / tol;
    let crsc = tol;
    let cssr = [cscl, CONER, crsc];
    let csrr = [crsc, CONER, cscl];
    let mut bry = [0.0f64; 3];
    bry[0] = 1.0e+3 * d1mach(1) / tol;

    let mut znr = zi;
    let zni = -zr;
    let zbr = zr;
    let mut zbi = zi;
    let mut cidi = -CONER;
    let inu = fnu as i32;
    let ang = HPI * (fnu - (inu as f64));
    let mut c2r = ang.cos();
    let mut c2i = ang.sin();
    let car = c2r;
    let sar = c2i;

    let mut idx = inu + n - 1;
    idx = idx.rem_euclid(4) + 1;
    let str = c2r * CIPR[(idx - 1) as usize] - c2i * CIPI[(idx - 1) as usize];
    c2i = c2r * CIPI[(idx - 1) as usize] + c2i * CIPR[(idx - 1) as usize];
    c2r = str;

    if zi <= 0.0 {
        znr = -znr;
        zbi = -zbi;
        cidi = -cidi;
        c2i = -c2i;
    }

    let mut fn_val = fnu.max(1.0);

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

    zunhj(
        znr, zni, fn_val, 1, tol, &mut phir, &mut phii, &mut argr, &mut argi,
        &mut zeta1r, &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut asumr, &mut asumi,
        &mut bsumr, &mut bsumi,
    );

    let (s1r, s1i);
    if kode != 1 {
        let mut str = zbr + zeta2r;
        let mut sti = zbi + zeta2i;
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

    'main_loop: loop {
        let nn = 2.min(nd);

        let mut i = 1;
        while i <= nn {
            fn_val = fnu + ((nd - i) as f64);
            zunhj(
                znr, zni, fn_val, 0, tol, &mut phir, &mut phii, &mut argr, &mut argi,
                &mut zeta1r, &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut asumr, &mut asumi,
                &mut bsumr, &mut bsumi,
            );

            let (mut s1r, mut s1i);
            if kode != 1 {
                let mut str = zbr + zeta2r;
                let mut sti = zbi + zeta2i;
                let rast = fn_val / zabs(str, sti);
                str = str * rast * rast;
                sti = -sti * rast * rast;
                s1r = -zeta1r + str;
                s1i = -zeta1i + sti + zi.abs();
            } else {
                s1r = -zeta1r + zeta2r;
                s1i = -zeta1i + zeta2i;
            }

            rs1 = s1r;
            if rs1.abs() > elim {
                match underflow_120(
                    rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu, kode,
                    tol, elim, alim, inu, car, sar, &mut c2r, &mut c2i,
                ) {
                    Flow::Return => return,
                    Flow::Continue => continue 'main_loop,
                    Flow::SetNlastAndReturn => return,
                }
            }

            if i == 1 {
                iflag = 2;
            }

            if rs1.abs() < alim {
                // On scale.
            } else {
                let aphi = zabs(phir, phii);
                let aarg = zabs(argr, argi);
                rs1 = rs1 + aphi.ln() - 0.25 * aarg.ln() - AIC;
                if rs1.abs() > elim {
                    match underflow_120(
                        rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu,
                        kode, tol, elim, alim, inu, car, sar, &mut c2r, &mut c2i,
                    ) {
                        Flow::Return => return,
                        Flow::Continue => continue 'main_loop,
                        Flow::SetNlastAndReturn => return,
                    }
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

            let mut air = 0.0f64;
            let mut aii = 0.0f64;
            let mut dair = 0.0f64;
            let mut daii = 0.0f64;
            let mut nai: i32 = 0;
            let mut ndai: i32 = 0;
            let mut idum: i32 = 0;
            zairy(argr, argi, 0, 2, &mut air, &mut aii, &mut nai, &mut idum);
            zairy(argr, argi, 1, 2, &mut dair, &mut daii, &mut ndai, &mut idum);

            let mut str = dair * bsumr - daii * bsumi;
            let mut sti = dair * bsumi + daii * bsumr;
            str += air * asumr - aii * asumi;
            sti += air * asumi + aii * asumr;

            let mut s2r = phir * str - phii * sti;
            let mut s2i = phir * sti + phii * str;

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
                    match underflow_120(
                        rs1, &mut nd, nz, &mut fn_val, fnul, nlast, yr, yi, zr, zi, fnu,
                        kode, tol, elim, alim, inu, car, sar, &mut c2r, &mut c2i,
                    ) {
                        Flow::Return => return,
                        Flow::Continue => continue 'main_loop,
                        Flow::SetNlastAndReturn => return,
                    }
                }
            }

            if zi <= 0.0 {
                s2i = -s2i;
            }
            let str = s2r * c2r - s2i * c2i;
            s2i = s2r * c2i + s2i * c2r;
            s2r = str;

            cyr[(i - 1) as usize] = s2r;
            cyi[(i - 1) as usize] = s2i;
            let j = nd - i + 1;
            yr[(j - 1) as usize] = s2r * csrr[(iflag - 1) as usize];
            yi[(j - 1) as usize] = s2i * csrr[(iflag - 1) as usize];

            let str = -c2i * cidi;
            c2i = c2r * cidi;
            c2r = str;

            i += 1;
        }

        if nd <= 2 {
            return;
        }

        let raz = 1.0 / zabs(zr, zi);
        let str = zr * raz;
        let sti = -zi * raz;
        let rzr = (str + str) * raz;
        let rzi = (sti + sti) * raz;
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

#[allow(clippy::too_many_arguments)]
fn underflow_120(
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
    inu: i32,
    car: f64,
    sar: f64,
    c2r: &mut f64,
    c2i: &mut f64,
) -> Flow {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
    const CIPI: [f64; 4] = [0.0, 1.0, 0.0, -1.0];

    if rs1 > 0.0 {
        *nz = -1;
        return Flow::Return;
    }

    yr[(*nd - 1) as usize] = ZEROR;
    yi[(*nd - 1) as usize] = ZEROI;
    *nz += 1;
    *nd -= 1;
    if *nd == 0 {
        return Flow::Return;
    }

    let mut nuf: i32 = 0;
    zuoik(zr, zi, fnu, kode, 1, *nd, yr, yi, &mut nuf, tol, elim, alim);
    if nuf < 0 {
        *nz = -1;
        return Flow::Return;
    }
    *nd -= nuf;
    *nz += nuf;
    if *nd == 0 {
        return Flow::Return;
    }

    *fn_val = fnu + ((*nd - 1) as f64);
    if *fn_val < fnul {
        *nlast = *nd;
        return Flow::SetNlastAndReturn;
    }

    let mut idx = inu + *nd - 1;
    idx = idx.rem_euclid(4);
    *c2r = car * CIPR[idx as usize] - sar * CIPI[idx as usize];
    *c2i = car * CIPI[idx as usize] + sar * CIPR[idx as usize];
    if zi <= 0.0 {
        *c2i = -*c2i;
    }

    Flow::Continue
}