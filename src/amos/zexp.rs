pub(crate) fn zexp(ar: f64, ai: f64, br: &mut f64, bi: &mut f64) {
    let zm = ar.exp();
    *br = zm * ai.cos();
    *bi = zm * ai.sin();
}