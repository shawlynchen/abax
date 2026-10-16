use crate::amos::zabs::zabs;
use crate::amos::zasyi::zasyi;
use crate::amos::zbuni::zbuni;
use crate::amos::zmlri::zmlri;
use crate::amos::zseri::zseri;
use crate::amos::zuoik::zuoik;
use crate::amos::zwrsk::zwrsk;

/// ZBINU computes the I function in the right half Z plane.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbinu(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    nz: &mut i32,
    rl: f64,
    fnul: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");

    *nz = 0;
    let az = zabs(zr, zi);
    let mut nn = n;
    let mut dfnu = fnu + ((n - 1) as f64);

    // -------------------------------------------------------------------
    // Power series dispatch.
    // -------------------------------------------------------------------
    let use_power_series = !(az > 2.0 && az * az * 0.25 > dfnu + 1.0);

    if use_power_series {
        let mut nw: i32 = 0;
        zseri(zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, tol, elim, alim);
        let inw = nw.abs();
        *nz += inw;
        nn -= inw;
        if nn == 0 {
            return;
        }
        if nw >= 0 {
            return;
        }
        dfnu = fnu + ((nn - 1) as f64);
    }

    // -------------------------------------------------------------------
    // Large-|Z| state machine starting at label 20.
    // -------------------------------------------------------------------
    let mut cwr = [0.0f64; 2];
    let mut cwi = [0.0f64; 2];
    let mut nw: i32 = 0;
    let mut nlast: i32 = 0;

    let mut d: i32 = 20;
    loop {
        match d {
            20 => {
                // Fortran label 20.
                if az < rl {
                    if dfnu <= 1.0 {
                        d = 70;
                    } else {
                        d = 50;
                    }
                } else if dfnu <= 1.0 {
                    d = 30;
                } else if az + az < dfnu * dfnu {
                    d = 50;
                } else {
                    d = 30;
                }
            }

            30 => {
                // Fortran label 30: asymptotic expansion for large Z.
                zasyi(zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, rl, tol, elim, alim);
                if nw < 0 {
                    d = 130;
                } else {
                    d = 120;
                }
            }

            50 => {
                // Fortran label 50: overflow/underflow test on I sequence.
                zuoik(zr, zi, fnu, kode, 1, nn, cyr, cyi, &mut nw, tol, elim, alim);
                if nw < 0 {
                    d = 130;
                } else {
                    *nz += nw;
                    nn -= nw;
                    if nn == 0 {
                        d = 120;
                    } else {
                        dfnu = fnu + ((nn - 1) as f64);
                        if dfnu > fnul || az > fnul {
                            d = 110;
                        } else {
                            d = 60;
                        }
                    }
                }
            }

            60 => {
                // Fortran label 60.
                if az > rl {
                    d = 80;
                } else {
                    d = 70;
                }
            }

            70 => {
                // Fortran label 70: Miller algorithm normalized by series.
                zmlri(zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, tol);
                if nw < 0 {
                    d = 130;
                } else {
                    d = 120;
                }
            }

            80 => {
                // Fortran label 80: Miller algorithm normalized by Wronskian.
                // First, overflow test on K functions.
                zuoik(zr, zi, fnu, kode, 2, 2, &mut cwr, &mut cwi, &mut nw, tol, elim, alim);
                if nw < 0 {
                    // Zero out and return (Fortran label 90 region).
                    *nz = nn;
                    for i in 0..(nn as usize) {
                        cyr[i] = 0.0;
                        cyi[i] = 0.0;
                    }
                    d = 120;
                } else if nw > 0 {
                    d = 130;
                } else {
                    zwrsk(
                        zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, &mut cwr, &mut cwi, tol,
                        elim, alim,
                    );
                    if nw < 0 {
                        d = 130;
                    } else {
                        d = 120;
                    }
                }
            }

            110 => {
                // Fortran label 110: increment FNU+NN-1 up to FNUL.
                let mut nui = (fnul - dfnu) as i32 + 1;
                nui = nui.max(0);
                zbuni(
                    zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, nui, &mut nlast, fnul, tol,
                    elim, alim,
                );
                if nw < 0 {
                    d = 130;
                } else {
                    *nz += nw;
                    if nlast == 0 {
                        d = 120;
                    } else {
                        nn = nlast;
                        d = 60;
                    }
                }
            }

            120 => {
                // Fortran label 120: normal return.
                return;
            }

            130 => {
                // Fortran label 130: error return.
                *nz = if nw == -2 { -2 } else { -1 };
                return;
            }

            _ => unreachable!("zbinu: invalid dispatch state {}", d),
        }
    }
}