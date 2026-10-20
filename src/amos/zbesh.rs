use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zacon::zacon;
use crate::amos::zbknu::zbknu;
use crate::amos::zbunk::zbunk;
use crate::amos::zuoik::zuoik;

/// ZBESH computes the H-Bessel (Hankel) functions of a complex argument.
///
/// On KODE=1, ZBESH computes an N-member sequence of complex Hankel
/// functions CY(J)=H(M,FNU+J-1,Z) for kinds M=1 or 2, real nonnegative
/// orders, and complex Z /= 0 in the cut plane -PI < ARG(Z) <= PI.
///
/// On KODE=2, ZBESH returns the scaled Hankel functions
///
/// CY(I) = EXP(-MM*Z*I)*H(M,FNU+J-1,Z),  MM = 3-2*M,  I**2 = -1
///
/// which removes the exponential behavior in both the upper and lower
/// half planes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbesh(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    m: i32,
    n: i32,
    cyr: &mut [f64],
    cyi: &mut [f64],
    nz: &mut i32,
    ierr: &mut i32,
) {
    const HPI: f64 = 1.57079632679489662;

    let n_usize = n as usize;
    assert!(cyr.len() >= n_usize, "cyr must have at least n elements");
    assert!(cyi.len() >= n_usize, "cyi must have at least n elements");

    *ierr = 0;
    *nz = 0;

    if zr == 0.0 && zi == 0.0 {
        *ierr = 1;
    }
    if fnu < 0.0 {
        *ierr = 1;
    }
    if m < 1 || m > 2 {
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

    let mut fn_val = fnu + ((nn - 1) as f64);
    let mm = 3 - m - m;
    let fmm = mm as f64;
    let mut znr = fmm * zi;
    let mut zni = -fmm * zr;

    // -------------------------------------------------------------------
    // Test for proper range.
    // -------------------------------------------------------------------
    let az = zabs(zr, zi);
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
                zuoik(znr, zni, fnu, kode, 2, nn, cyr, cyi, &mut nuf, tol, elim, alim);
                if nuf < 0 {
                    *nz = 0;
                    *ierr = 2;
                    return;
                }
                *nz += nuf;
                nn -= nuf;
                if nn == 0 {
                    if znr < 0.0 {
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
        if znr < 0.0 || (znr == 0.0 && zni < 0.0 && m == 2) {
            // Left half plane computation.
            let mr = -mm;
            let mut nw: i32 = 0;
            zacon(znr, zni, fnu, kode, mr, nn, cyr, cyi, &mut nw, rl, fnul, tol, elim, alim);
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
        } else {
            // Right half plane computation.
            zbknu(znr, zni, fnu, kode, nn, cyr, cyi, nz, tol, elim, alim);
        }
    } else {
        // -------------------------------------------------------------------
        // FNU > FNUL: uniform asymptotic expansions.
        // -------------------------------------------------------------------
        let mut mr = 0;
        if !((znr >= 0.0) && (znr != 0.0 || zni >= 0.0 || m != 2)) {
            mr = -mm;
            if znr == 0.0 && zni < 0.0 {
                znr = -znr;
                zni = -zni;
            }
        }

        let mut nw: i32 = 0;
        zbunk(znr, zni, fnu, kode, mr, nn, cyr, cyi, &mut nw, tol, elim, alim);
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

    // -------------------------------------------------------------------
    // H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
    // ZT = EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
    // -------------------------------------------------------------------
    let sgn = HPI.copysign(-fmm);

    // -------------------------------------------------------------------
    // Calculate EXP(FNU*HPI*I) to minimize losses of significance
    // when FNU is large.
    // -------------------------------------------------------------------
    let inu = fnu as i32;
    let inuh = inu / 2;
    let ir = inu - 2 * inuh;
    let arg = (fnu - ((inu - ir) as f64)) * sgn;
    let rhpi = 1.0 / sgn;
    let mut csgni = rhpi * arg.cos();
    let mut csgnr = -rhpi * arg.sin();
    if inuh % 2 != 0 {
        csgnr = -csgnr;
        csgni = -csgni;
    }

    let zti = -fmm;
    let rtol = 1.0 / tol;
    let ascle = ufl * rtol;

    for i in 0..(nn as usize) {
        let mut aa = cyr[i];
        let mut bb = cyi[i];
        let mut atol = 1.0;
        if aa.abs().max(bb.abs()) <= ascle {
            aa *= rtol;
            bb *= rtol;
            atol = tol;
        }
        let str = aa * csgnr - bb * csgni;
        let sti = aa * csgni + bb * csgnr;
        cyr[i] = str * atol;
        cyi[i] = sti * atol;
        let str = -csgni * zti;
        csgni = csgnr * zti;
        csgnr = str;
    }
}