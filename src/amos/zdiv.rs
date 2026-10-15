use crate::amos::zabs;

pub(crate) fn zdiv(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64) {
    // BM = 1 / |B|
    let bm = 1.0_f64 / zabs(br, bi);

    // Normalize B by |B|
    let cc = br * bm; // real part of B / |B|
    let cd = bi * bm; // imag part of B / |B|

    // C = (A * conj(B)) / |B|^2, computed with the normalized B
    *cr = (ar * cc + ai * cd) * bm;
    *ci = (ai * cc - ar * cd) * bm;
}