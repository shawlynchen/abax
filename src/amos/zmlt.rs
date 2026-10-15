pub(crate) fn zmlt(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64) {
    *cr = ar * br - ai * bi;
    *ci = ar * bi + ai * br;
}