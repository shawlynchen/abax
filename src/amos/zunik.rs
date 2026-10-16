use crate::amos::d1mach::d1mach;
use crate::amos::zdiv::zdiv;
use crate::amos::zlog::zlog;
use crate::amos::zsqrt::zsqrt;

/// ZUNIK computes parameters for the uniform asymptotic
/// expansions of the I and K functions on IKFLG = 1 or 2
/// respectively by
///
///     W(FNU,ZR) = PHI*EXP(ZETA)*SUM
///
/// where  ZETA = -ZETA1 + ZETA2  or  ZETA1 - ZETA2.
///
/// The first call must have INIT=0. Subsequent calls with the
/// same ZR and FNU will return the I or K function on IKFLG =
/// 1 or 2 with no change in INIT. CWRK is a complex work
/// array. IPMTR=0 computes all parameters. IPMTR=1 computes PHI,
/// ZETA1, ZETA2.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zunik(
    zrr: f64,
    zri: f64,
    fnu: f64,
    ikflg: i32,
    ipmtr: i32,
    tol: f64,
    init: &mut i32,
    phir: &mut f64,
    phii: &mut f64,
    zeta1r: &mut f64,
    zeta1i: &mut f64,
    zeta2r: &mut f64,
    zeta2i: &mut f64,
    sumr: &mut f64,
    sumi: &mut f64,
    cwrkr: &mut [f64; 16],
    cwrki: &mut [f64; 16],
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;
    const CONER: f64 = 1.0;
    const CONEI: f64 = 0.0;

    // CON(1), CON(2)
    const CON: [f64; 2] = [3.98942280401432678e-01, 1.25331413731550025e+00];

    // Coefficients of the asymptotic expansion (1-based in Fortran;
    // shifted to 0-based indexing here).
    const C: [f64; 120] = [
        1.00000000000000000e+00,  -2.08333333333333333e-01,
        1.25000000000000000e-01,   3.34201388888888889e-01,
       -4.01041666666666667e-01,   7.03125000000000000e-02,
       -1.02581259645061728e+00,   1.84646267361111111e+00,
       -8.91210937500000000e-01,   7.32421875000000000e-02,
        4.66958442342624743e+00,  -1.12070026162229938e+01,
        8.78912353515625000e+00,  -2.36408691406250000e+00,
        1.12152099609375000e-01,  -2.82120725582002449e+01,
        8.46362176746007346e+01,  -9.18182415432400174e+01,
        4.25349987453884549e+01,  -7.36879435947963170e+00,
        2.27108001708984375e-01,   2.12570130039217123e+02,
       -7.65252468141181642e+02,   1.05999045252799988e+03,
       -6.99579627376132541e+02,   2.18190511744211590e+02,
       -2.64914304869515555e+01,   5.72501420974731445e-01,
       -1.91945766231840700e+03,   8.06172218173730938e+03,
       -1.35865500064341374e+04,   1.16553933368645332e+04,
       -5.30564697861340311e+03,   1.20090291321635246e+03,
       -1.08090919788394656e+02,   1.72772750258445740e+00,
        2.02042913309661486e+04,  -9.69805983886375135e+04,
        1.92547001232531532e+05,  -2.03400177280415534e+05,
        1.22200464983017460e+05,  -4.11926549688975513e+04,
        7.10951430248936372e+03,  -4.93915304773088012e+02,
        6.07404200127348304e+00,  -2.42919187900551333e+05,
        1.31176361466297720e+06,  -2.99801591853810675e+06,
        3.76327129765640400e+06,  -2.81356322658653411e+06,
        1.26836527332162478e+06,  -3.31645172484563578e+05,
        4.52187689813627263e+04,  -2.49983048181120962e+03,
        2.43805296995560639e+01,   3.28446985307203782e+06,
       -1.97068191184322269e+07,   5.09526024926646422e+07,
       -7.41051482115326577e+07,   6.63445122747290267e+07,
       -3.75671766607633513e+07,   1.32887671664218183e+07,
       -2.78561812808645469e+06,   3.08186404612662398e+05,
       -1.38860897537170405e+04,   1.10017140269246738e+02,
       -4.93292536645099620e+07,   3.25573074185765749e+08,
       -9.39462359681578403e+08,   1.55359689957058006e+09,
       -1.62108055210833708e+09,   1.10684281682301447e+09,
       -4.95889784275030309e+08,   1.42062907797533095e+08,
       -2.44740627257387285e+07,   2.24376817792244943e+06,
       -8.40054336030240853e+04,   5.51335896122020586e+02,
        8.14789096118312115e+08,  -5.86648149205184723e+09,
        1.86882075092958249e+10,  -3.46320433881587779e+10,
        4.12801855797539740e+10,  -3.30265997498007231e+10,
        1.79542137311556001e+10,  -6.56329379261928433e+09,
        1.55927986487925751e+09,  -2.25105661889415278e+08,
        1.73951075539781645e+07,  -5.49842327572288687e+05,
        3.03809051092238427e+03,  -1.46792612476956167e+10,
        1.14498237732025810e+11,  -3.99096175224466498e+11,
        8.19218669548577329e+11,  -1.09837515608122331e+12,
        1.00815810686538209e+12,  -6.45364869245376503e+11,
        2.87900649906150589e+11,  -8.78670721780232657e+10,
        1.76347306068349694e+10,  -2.16716498322379509e+09,
        1.43157876718888981e+08,  -3.87183344257261262e+06,
        1.82577554742931747e+04,   2.86464035717679043e+11,
       -2.40629790002850396e+12,   9.10934118523989896e+12,
       -2.05168994109344374e+13,   3.05651255199353206e+13,
       -3.16670885847851584e+13,   2.33483640445818409e+13,
       -1.23204913055982872e+13,   4.61272578084913197e+12,
       -1.19655288019618160e+12,   2.05914503232410016e+11,
       -2.18229277575292237e+10,   1.24700929351271032e+09,
       -2.91883881222208134e+07,   1.18838426256783253e+05,
    ];

    if *init == 0 {
        // -----------------------------------------------------------------
        // Initialize all variables
        // -----------------------------------------------------------------
        let rfn = 1.0 / fnu;

        // -----------------------------------------------------------------
        // Overflow test (ZR/FNU too small)
        // -----------------------------------------------------------------
        let test = d1mach(1) * 1.0e+3;
        let ac = fnu * test;
        if zrr.abs() <= ac && zri.abs() <= ac {
            *zeta1r = 2.0 * test.abs().ln() + fnu;
            *zeta1i = 0.0;
            *zeta2r = fnu;
            *zeta2i = 0.0;
            *phir = 1.0;
            *phii = 0.0;
            return;
        }

        let tr = zrr * rfn;
        let ti = zri * rfn;
        let sr = CONER + (tr * tr - ti * ti);
        let si = CONEI + (tr * ti + ti * tr);

        let (mut srr, mut sri) = (0.0f64, 0.0f64);
        zsqrt(sr, si, &mut srr, &mut sri);

        let str = CONER + srr;
        let sti = CONEI + sri;

        let (mut znr, mut zni) = (0.0f64, 0.0f64);
        zdiv(str, sti, tr, ti, &mut znr, &mut zni);

        let (mut str, mut sti) = (0.0f64, 0.0f64);
        let _ = zlog(znr, zni, &mut str, &mut sti);
        *zeta1r = fnu * str;
        *zeta1i = fnu * sti;
        *zeta2r = fnu * srr;
        *zeta2i = fnu * sri;

        let (mut tr, mut ti) = (0.0f64, 0.0f64);
        zdiv(CONER, CONEI, srr, sri, &mut tr, &mut ti);
        let srr = tr * rfn;
        let sri = ti * rfn;

        zsqrt(srr, sri, &mut cwrkr[15], &mut cwrki[15]);
        *phir = cwrkr[15] * CON[(ikflg - 1) as usize];
        *phii = cwrki[15] * CON[(ikflg - 1) as usize];

        if ipmtr != 0 {
            return;
        }

        let (mut t2r, mut t2i) = (0.0f64, 0.0f64);
        zdiv(CONER, CONEI, sr, si, &mut t2r, &mut t2i);

        cwrkr[0] = CONER;
        cwrki[0] = CONEI;
        let mut crfnr = CONER;
        let mut crfni = CONEI;
        let mut ac = 1.0;
        let mut l: usize = 0; // 0-based index into C; Fortran L starts at 1

        *init = 15;
        for k in 2..=15 {
            let mut sr = ZEROR;
            let mut si = ZEROI;
            for _j in 1..=k {
                l += 1;
                let str = sr * t2r - si * t2i + C[l - 1]; // Fortran C(L), 1-based
                si = sr * t2i + si * t2r;
                sr = str;
            }
            let str = crfnr * srr - crfni * sri;
            crfni = crfnr * sri + crfni * srr;
            crfnr = str;
            cwrkr[k - 1] = crfnr * sr - crfni * si;
            cwrki[k - 1] = crfnr * si + crfni * sr;
            ac *= rfn;
            let test = cwrkr[k - 1].abs() + cwrki[k - 1].abs();
            if ac < tol && test < tol {
                *init = k as i32;
                break;
            }
        }
    }

    // -------------------------------------------------------------------
    // Compute sum for the I function (IKFLG == 1)
    // -------------------------------------------------------------------
    if ikflg != 2 {
        let mut sr = ZEROR;
        let mut si = ZEROI;
        for i in 0..(*init as usize) {
            sr += cwrkr[i];
            si += cwrki[i];
        }
        *sumr = sr;
        *sumi = si;
        *phir = cwrkr[15] * CON[0];
        *phii = cwrki[15] * CON[0];
        return;
    }

    // -------------------------------------------------------------------
    // Compute sum for the K function (IKFLG == 2)
    // -------------------------------------------------------------------
    let mut sr = ZEROR;
    let mut si = ZEROI;
    let mut tr = CONER;
    for i in 0..(*init as usize) {
        sr += tr * cwrkr[i];
        si += tr * cwrki[i];
        tr = -tr;
    }
    *sumr = sr;
    *sumi = si;
    *phir = cwrkr[15] * CON[1];
    *phii = cwrki[15] * CON[1];
}