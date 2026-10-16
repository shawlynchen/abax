use crate::amos::zabs::zabs;
use crate::amos::zexp::zexp;
use crate::amos::zlog::zlog;

/// ZS1S2 tests for a possible underflow resulting from the
/// addition of the I and K functions in the analytic continuation
/// formula where S1 = K function and S2 = I function.
/// On KODE=1 the I and K functions are different orders of
/// magnitude, but for KODE=2 they can be of the same order
/// of magnitude and the maximum must be at least one
/// precision above the underflow limit.
///
/// * `zrr, zri`  – real and imaginary part of ZR
/// * `s1r, s1i`  – real/imag part of S1 (inout)
/// * `s2r, s2i`  – real/imag part of S2 (inout)
/// * `nz`        – underflow count (out)
/// * `ascle`     – underflow / overflow threshold
/// * `alim`      – argument limit
/// * `iuf`       – underflow flag (inout)
pub(crate) fn zs1s2(
    zrr: f64,
    zri: f64,
    s1r: &mut f64,
    s1i: &mut f64,
    s2r: &mut f64,
    s2i: &mut f64,
    nz: &mut i32,
    ascle: f64,
    alim: f64,
    iuf: &mut i32,
) {
    const ZEROR: f64 = 0.0;
    const ZEROI: f64 = 0.0;

    *nz = 0;

    let mut as1 = zabs(*s1r, *s1i);
    let as2 = zabs(*s2r, *s2i);

    // First branch: handle the S1 underflow test.
    // Fortran:
    //   IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
    //   IF (AS1.EQ.0.0D0) GO TO 10
    //   ... body ...
    //   IF (ALN.LT.(-ALIM)) GO TO 10
    //   ... inner body ...
    //  10 CONTINUE
    if !(*s1r == ZEROR && *s1i == ZEROI) && as1 != ZEROR {
        let aln = -zrr - zrr + as1.ln();
        let s1dr = *s1r;
        let s1di = *s1i;
        *s1r = ZEROR;
        *s1i = ZEROI;
        as1 = ZEROR;

        if aln >= -alim {
            let (mut c1r, mut c1i) = (0.0f64, 0.0f64);
            let _ = zlog(s1dr, s1di, &mut c1r, &mut c1i);
            c1r = c1r - zrr - zrr;
            c1i = c1i - zri - zri;

            zexp(c1r, c1i, s1r, s1i);
            as1 = zabs(*s1r, *s1i);
            *iuf += 1;
        }
    }

    // Final underflow check on the max of AS1, AS2.
    let aa = as1.max(as2);
    if aa > ascle {
        return;
    }

    *s1r = ZEROR;
    *s1i = ZEROI;
    *s2r = ZEROR;
    *s2i = ZEROI;
    *nz = 1;
    *iuf = 0;
}