use crate::amos::zabs;

pub(crate) fn zsqrt(ar: f64, ai: f64, br: &mut f64, bi: &mut f64) {
    const DRT: f64 = 7.071067811865475244008443621e-1_f64;
    const DPI: f64 = 3.141592653589793238462643383_f64;

    let zm = zabs(ar, ai).sqrt();

    // AR == 0
    if ar == 0.0 {
        if ai > 0.0 {
            *br = zm * DRT;
            *bi = zm * DRT;
            return;
        }

        if ai < 0.0 {
            *br = zm * DRT;
            *bi = -zm * DRT;
            return;
        }

        *br = 0.0;
        *bi = 0.0;
        return;
    }

    // AI == 0
    if ai == 0.0 {
        if ar > 0.0 {
            *br = ar.sqrt();
            *bi = 0.0;
        } else {
            *br = 0.0;
            *bi = ar.abs().sqrt();
        }
        return;
    }

    // General case
    let mut dtheta = (ai / ar).atan();

    if dtheta <= 0.0 {
        if ar < 0.0 {
            dtheta += DPI;
        }
    } else if ar < 0.0 {
        dtheta -= DPI;
    }

    dtheta *= 0.5;

    *br = zm * dtheta.cos();
    *bi = zm * dtheta.sin();
}