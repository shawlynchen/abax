use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zbknu::zbknu;
use crate::amos::zrati::zrati;

/// ZWRSK computes the I Bessel function for Re(Z) >= 0 by
/// normalizing the I function ratios from ZRATI by the Wronskian.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zwrsk(
    zrr: f64,
    zri: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    cwr: &mut [f64; 2],
    cwi: &mut [f64; 2],
    tol: f64,
    elim: f64,
    alim: f64,
) {
    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    // -------------------------------------------------------------------
    // I(FNU+I-1,Z) by backward recurrence for ratios
    // Y(I) = I(FNU+I,Z)/I(FNU+I-1,Z) from CRATI normalized by the
    // Wronskian with K(FNU,Z) and K(FNU+1,Z) from CBKNU.
    // -------------------------------------------------------------------
    *nz = 0;
    let mut nw: i32 = 0;
    zbknu(zrr, zri, fnu, kode, 2, cwr, cwi, &mut nw, tol, elim, alim);
    if nw != 0 {
        *nz = -1;
        if nw == -2 {
            *nz = -2;
        }
        return;
    }

    zrati(zrr, zri, fnu, n, yr, yi, tol);

    // -------------------------------------------------------------------
    // Recur forward on I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    // R(FNU+J-1,Z) = Y(J), J = 1,...,N
    // -------------------------------------------------------------------
    let (mut cinur, mut cinui) = if kode != 1 {
        (zri.cos(), zri.sin())
    } else {
        (1.0, 0.0)
    };

    // -------------------------------------------------------------------
    // On low exponent machines the K functions can be close to both
    // the under and overflow limits and the normalization must be
    // scaled to prevent over or underflow. ZUOIK has determined that
    // the result is on scale.
    // -------------------------------------------------------------------
    let acw = zabs(cwr[1], cwi[1]);
    let mut ascle = 1.0e+3 * d1mach(1) / tol;
    let csclr;

    if acw <= ascle {
        csclr = 1.0 / tol;
    } else {
        ascle = 1.0 / ascle;
        if acw >= ascle {
            csclr = tol;
        } else {
            csclr = 1.0;
        }
    }

    let c1r = cwr[0] * csclr;
    let c1i = cwi[0] * csclr;
    let c2r = cwr[1] * csclr;
    let c2i = cwi[1] * csclr;
    let mut str = yr[0];
    let mut sti = yi[0];

    // -------------------------------------------------------------------
    // CINU = CINU*(CONJG(CT)/|CT|)*(1/|CT|) prevents under- or overflow
    // prematurely by squaring |CT|.
    // -------------------------------------------------------------------
    let mut ptr = str * c1r - sti * c1i;
    let mut pti = str * c1i + sti * c1r;
    ptr += c2r;
    pti += c2i;
    let mut ctr = zrr * ptr - zri * pti;
    let mut cti = zrr * pti + zri * ptr;
    let act = zabs(ctr, cti);
    let ract = 1.0 / act;
    ctr *= ract;
    cti = -cti * ract;
    ptr = cinur * ract;
    pti = cinui * ract;
    cinur = ptr * ctr - pti * cti;
    cinui = ptr * cti + pti * ctr;
    yr[0] = cinur * csclr;
    yi[0] = cinui * csclr;

    if n == 1 {
        return;
    }

    for i in 1..n_usize {
        ptr = str * cinur - sti * cinui;
        cinui = str * cinui + sti * cinur;
        cinur = ptr;
        str = yr[i];
        sti = yi[i];
        yr[i] = cinur * csclr;
        yi[i] = cinui * csclr;
    }
}