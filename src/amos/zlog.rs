use crate::amos::zabs;

pub(crate) fn zlog(ar: f64, ai: f64, br: &mut f64, bi: &mut f64) -> i32 {
    const DPI: f64 = 3.141592653589793238462643383_f64;
    const DHPI: f64 = 1.570796326794896619231321696_f64;

    if ar == 0.0 {
        if ai == 0.0 {
            return 1;
        }

        *bi = if ai < 0.0 { -DHPI } else { DHPI };
        *br = ai.abs().ln();
        return 0;
    }

    if ai == 0.0 {
        if ar > 0.0 {
            *br = ar.ln();
            *bi = 0.0;
        } else {
            *br = ar.abs().ln();
            *bi = DPI;
        }
        return 0;
    }

    let mut dtheta = (ai / ar).atan();

    if dtheta > 0.0 {
        if ar < 0.0 {
            dtheta -= DPI;
        }
    } else if ar < 0.0 {
        dtheta += DPI;
    }

    let zm = zabs(ar, ai);
    *br = zm.ln();
    *bi = dtheta;

    0
}