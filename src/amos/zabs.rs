// f64::hypot(zr, zi)
pub(crate) fn zabs(zr: f64, zi: f64) -> f64 {
    let u = zr.abs();
    let v = zi.abs();
    let s = u + v;

    if s == 0.0 {
        return 0.0;
    }

    if u > v {
        let q = v / u;
        u * (1.0 + q * q).sqrt()
    } else {
        let q = u / v;
        v * (1.0 + q * q).sqrt()
    }
}