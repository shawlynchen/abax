/// Tests whether a scaled complex quantity `Y = (yr, yi)` would
/// underflow when scaled (by `tol`) to its proper value.
///
/// `Y` enters as a scaled quantity whose magnitude is greater than
/// `exp(-alim) = ascle = 1.0e+3 * D1MACH(1) / tol`. The test checks
/// whether the magnitude of the real or imaginary part would underflow
/// when `Y` is scaled (by `tol`) to its proper value. `Y` is accepted
/// if the underflow is at least one precision below the magnitude of
/// the largest component; otherwise the phase angle does not have
/// absolute accuracy and an underflow is assumed.
///
/// Returns `nz`:
/// - `nz == 0` if `Y` is acceptable (no underflow),
/// - `nz == 1` if an underflow is assumed.
///
/// This is a Rust port of the SLATEC `ZUCHK` routine.
pub fn zuchk(yr: f64, yi: f64, nz: &mut i32, ascle: f64, tol: f64) {
    *nz = 0;

    let wr = yr.abs();
    let wi = yi.abs();

    // ST = min(|yr|, |yi|)
    let mut st = wr.min(wi);

    // If the smaller component is already above ASCLE, no underflow.
    if st > ascle {
        return;
    }

    // SS = max(|yr|, |yi|)
    let ss = wr.max(wi);

    // ST = min(|yr|, |yi|) / TOL  (the component after scaling back)
    st /= tol;

    // If the largest component (in scaled form) is smaller than the
    // scaled-down smaller component, the phase angle is unreliable.
    if ss < st {
        *nz = 1;
    }
}
