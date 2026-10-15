pub(crate) fn zshch(zr: f64, zi: f64, cshr: &mut f64, cshi: &mut f64, cchr: &mut f64, cchi: &mut f64) {
    let sh = zr.sinh();
    let ch = zr.cosh();
    let sn = zi.sin();
    let cn = zi.cos();

    *cshr = sh * cn;
    *cshi = ch * sn;
    *cchr = ch * cn;
    *cchi = sh * sn;
}