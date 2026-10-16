use crate::amos::zunk1::zunk1;
use crate::amos::zunk2::zunk2;

/// ZBUNK computes the K Bessel function for FNU > FNUL.
/// According to the uniform asymptotic expansion for K(FNU,Z)
/// in ZUNK1 and the expansion for H(2,FNU,Z) in ZUNK2.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbunk(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    mr: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
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

    if ay <= ax {
        // Asymptotic expansion for K(FNU,Z) for large FNU applied in
        // -PI/3 <= ARG(Z) <= PI/3.
        zunk1(zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim, alim);
    } else {
        // Asymptotic expansion for H(2,FNU,Z*EXP(M*HPI)) for large FNU
        // applied in PI/3 < |ARG(Z)| <= PI/2 where M=+I or -I and HPI=PI/2.
        zunk2(zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim, alim);
    }
}