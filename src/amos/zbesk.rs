use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zacon::zacon;
use crate::amos::zbknu::zbknu;
use crate::amos::zbunk::zbunk;
use crate::amos::zuoik::zuoik;

/// ZBESK computes K-Bessel functions of complex argument.
///
/// On KODE=1, ZBESK computes an N-member sequence of complex Bessel
/// functions CY(J)=K(FNU+J-1,Z) for real nonnegative orders and
/// complex Z /= 0 in the cut plane -PI < ARG(Z) <= PI.
///
/// On KODE=2, ZBESK returns the scaled K functions
///
/// CY(J) = EXP(Z)*K(FNU+J-1,Z),  J=1,...,N
///
/// which remove the exponential behavior in both the left and right
/// half planes for Z to infinity.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbesk(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    nz: &mut i32,
    ierr: &mut i32,
) {
    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");

    *ierr = 0;
    *nz = 0;

    if zi == 0.0 && zr == 0.0 {
        *ierr = 1;
    }
    if fnu < 0.0 {
        *ierr = 1;
    }
    if kode < 1 || kode > 2 {
        *ierr = 1;
    }
    if n < 1 {
        *ierr = 1;
    }
    if *ierr != 0 {
        return;
    }

    let mut nn = n;

    // -------------------------------------------------------------------
    // Set parameters related to machine constants.
    // -------------------------------------------------------------------
    let tol = d1mach(4).max(1.0e-18);
    let k1 = i1mach(15);
    let k2 = i1mach(16);
    let r1m5 = d1mach(5);
    let k = k1.abs().min(k2.abs());
    let elim = 2.303 * ((k as f64) * r1m5 - 3.0);
    let k1b = i1mach(14) - 1;
    let mut aa = r1m5 * (k1b as f64);
    let dig = aa.min(18.0);
    aa *= 2.303;
    let alim = elim + (-aa).max(-41.45);
    let fnul = 10.0 + 6.0 * (dig - 3.0);
    let rl = 1.2 * dig + 3.0;

    // -------------------------------------------------------------------
    // Test for proper range.
    // -------------------------------------------------------------------
    let az = zabs(zr, zi);
    let fn_val = fnu + ((nn - 1) as f64);
    let mut aa = 0.5 / tol;
    let bb = (i1mach(9) as f64) * 0.5;
    aa = aa.min(bb);
    if az > aa {
        *nz = 0;
        *ierr = 4;
        return;
    }
    if fn_val > aa {
        *nz = 0;
        *ierr = 4;
        return;
    }
    aa = aa.sqrt();
    if az > aa {
        *ierr = 3;
    }
    if fn_val > aa {
        *ierr = 3;
    }

    // -------------------------------------------------------------------
    // Overflow test on the last member of the sequence.
    // -------------------------------------------------------------------
    let ufl = d1mach(1) * 1.0e+3;

    if az < ufl {
        // Too small.
        *nz = 0;
        *ierr = 2;
        return;
    }

    if fnu <= fnul {
        // Regular path.
        if fn_val > 1.0 {
            if fn_val > 2.0 {
                // ZUOIK test.
                let mut nuf: i32 = 0;
                zuoik(zr, zi, fnu, kode, 2, nn, cyr, cyi, &mut nuf, tol, elim, alim);
                if nuf < 0 {
                    *nz = 0;
                    *ierr = 2;
                    return;
                }
                *nz += nuf;
                nn -= nuf;
                if nn == 0 {
                    if zr < 0.0 {
                        *nz = 0;
                        *ierr = 2;
                        return;
                    }
                    return;
                }
            } else {
                // fn in (1, 2]: test for small |Z|.
                if az <= tol {
                    let arg = 0.5 * az;
                    let aln = -fn_val * arg.ln();
                    if aln > elim {
                        *nz = 0;
                        *ierr = 2;
                        return;
                    }
                }
            }
        }

        // Right/left half plane dispatch.
        if zr >= 0.0 {
            // Right half plane computation.
            let mut nw: i32 = 0;
            zbknu(zr, zi, fnu, kode, nn, cyr, cyi, &mut nw, tol, elim, alim);
            if nw < 0 {
                if nw == -1 {
                    *nz = 0;
                    *ierr = 2;
                    return;
                }
                *nz = 0;
                *ierr = 5;
                return;
            }
            *nz = nw;
            return;
        }

        // Left half plane computation.
        if *nz != 0 {
            *nz = 0;
            *ierr = 2;
            return;
        }

        let mr = if zi < 0.0 { -1 } else { 1 };
        let mut nw: i32 = 0;
        zacon(zr, zi, fnu, kode, mr, nn, cyr, cyi, &mut nw, rl, fnul, tol, elim, alim);
        if nw < 0 {
            if nw == -1 {
                *nz = 0;
                *ierr = 2;
                return;
            }
            *nz = 0;
            *ierr = 5;
            return;
        }
        *nz = nw;
        return;
    }

    // -------------------------------------------------------------------
    // FNU > FNUL: uniform asymptotic expansions.
    // -------------------------------------------------------------------
    let mr = if zr < 0.0 {
        if zi < 0.0 { -1 } else { 1 }
    } else {
        0
    };

    let mut nw: i32 = 0;
    zbunk(zr, zi, fnu, kode, mr, nn, cyr, cyi, &mut nw, tol, elim, alim);
    if nw < 0 {
        if nw == -1 {
            *nz = 0;
            *ierr = 2;
            return;
        }
        *nz = 0;
        *ierr = 5;
        return;
    }
    *nz += nw;
}