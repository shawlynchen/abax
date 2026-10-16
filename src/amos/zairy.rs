use crate::amos::d1mach::d1mach;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zacai::zacai;
use crate::amos::zbknu::zbknu;
use crate::amos::zexp::zexp;
use crate::amos::zsqrt::zsqrt;

// ---------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------

const TTH: f64 = 6.66666666666666667e-01;
const C1: f64 = 3.55028053887817240e-01;
const C2: f64 = 2.58819403792806799e-01;
const COEF: f64 = 1.83776298473930683e-01;

const ZEROR: f64 = 0.0;
const ZEROI: f64 = 0.0;
const CONER: f64 = 1.0;
const CONEI: f64 = 0.0;

// ---------------------------------------------------------------------
// State shared across helpers.
// ---------------------------------------------------------------------

struct ZairyState {
    zr: f64,
    zi: f64,
    id: i32,
    kode: i32,

    // Outputs
    air: f64,
    aii: f64,
    nz: i32,
    ierr: i32,

    // Shared working values
    az: f64,
    tol: f64,
    fid: f64,
    csqr: f64,
    csqi: f64,
    iflag: i32,
    sfac: f64,
}

/// ZAIRY computes the complex Airy function AI(Z) or its derivative
/// DAI(Z)/DZ for complex Z.
///
/// On KODE=1, ZAIRY computes AI(Z) or DAI(Z)/DZ on ID=0 or ID=1
/// respectively. On KODE=2, a scaling option CEXP(ZTA)*AI(Z) or
/// CEXP(ZTA)*DAI(Z)/DZ is provided to remove the exponential decay
/// in -PI/3 < ARG(Z) < PI/3 and the exponential growth in
/// PI/3 < |ARG(Z)| < PI where ZTA = (2/3)*Z*CSQRT(Z).
pub(crate) fn zairy(
    zr: f64,
    zi: f64,
    id: i32,
    kode: i32,
    air: &mut f64,
    aii: &mut f64,
    nz: &mut i32,
    ierr: &mut i32,
) {
    let mut s = ZairyState {
        zr,
        zi,
        id,
        kode,
        air: 0.0,
        aii: 0.0,
        nz: 0,
        ierr: 0,
        az: 0.0,
        tol: 0.0,
        fid: 0.0,
        csqr: 0.0,
        csqi: 0.0,
        iflag: 0,
        sfac: 1.0,
    };

    s.run();

    *air = s.air;
    *aii = s.aii;
    *nz = s.nz;
    *ierr = s.ierr;
}

impl ZairyState {
    fn run(&mut self) {
        self.ierr = 0;
        self.nz = 0;

        if self.id < 0 || self.id > 1 {
            self.ierr = 1;
        }
        if self.kode < 1 || self.kode > 2 {
            self.ierr = 1;
        }
        if self.ierr != 0 {
            return;
        }

        self.az = zabs(self.zr, self.zi);
        self.tol = d1mach(4).max(1.0e-18);
        self.fid = self.id as f64;

        if self.az > 1.0 {
            self.large_z();
        } else {
            self.power_series();
        }
    }

    // -----------------------------------------------------------------
    // Power series for |Z| <= 1.
    // -----------------------------------------------------------------
    fn power_series(&mut self) {
        let mut s1r = CONER;
        let mut s1i = CONEI;
        let mut s2r = CONER;
        let mut s2i = CONEI;

        if self.az < self.tol {
            self.tiny();
            return;
        }

        let aa = self.az * self.az;
        if aa < self.tol / self.az {
            self.series_end(s1r, s1i, s2r, s2i);
            return;
        }

        let mut trm1r = CONER;
        let mut trm1i = CONEI;
        let mut trm2r = CONER;
        let mut trm2i = CONEI;
        let mut atrm = 1.0;
        let mut str = self.zr * self.zr - self.zi * self.zi;
        let mut sti = self.zr * self.zi + self.zi * self.zr;
        let z3r = str * self.zr - sti * self.zi;
        let z3i = str * self.zi + sti * self.zr;
        let az3 = self.az * aa;

        let mut ak = 2.0 + self.fid;
        let mut bk = 3.0 - self.fid - self.fid;
        let ck = 4.0 - self.fid;
        let dk = 3.0 + self.fid + self.fid;
        let mut d1 = ak * dk;
        let mut d2 = bk * ck;
        let mut ad = d1.min(d2);
        ak = 24.0 + 9.0 * self.fid;
        bk = 30.0 - 9.0 * self.fid;

        for _k in 1..=25 {
            str = (trm1r * z3r - trm1i * z3i) / d1;
            trm1i = (trm1r * z3i + trm1i * z3r) / d1;
            trm1r = str;
            s1r += trm1r;
            s1i += trm1i;

            str = (trm2r * z3r - trm2i * z3i) / d2;
            trm2i = (trm2r * z3i + trm2i * z3r) / d2;
            trm2r = str;
            s2r += trm2r;
            s2i += trm2i;

            atrm = atrm * az3 / ad;
            d1 += ak;
            d2 += bk;
            ad = d1.min(d2);

            if atrm < self.tol * ad {
                break;
            }
            ak += 18.0;
            bk += 18.0;
        }

        self.series_end(s1r, s1i, s2r, s2i);
    }

    // -----------------------------------------------------------------
    // |Z| < TOL branch (Fortran label 170).
    // -----------------------------------------------------------------
    fn tiny(&mut self) {
        let aa = 1.0e+3 * d1mach(1);
        let mut s1r = ZEROR;
        let mut s1i = ZEROI;

        if self.id != 1 {
            if self.az > aa {
                s1r = C2 * self.zr;
                s1i = C2 * self.zi;
            }
            self.air = C1 - s1r;
            self.aii = -s1i;
            return;
        }

        self.air = -C2;
        self.aii = 0.0;
        let aa = aa.sqrt();
        if self.az > aa {
            s1r = 0.5 * (self.zr * self.zr - self.zi * self.zi);
            s1i = self.zr * self.zi;
        }
        self.air += C1 * s1r;
        self.aii += C1 * s1i;
    }

    // -----------------------------------------------------------------
    // Series end (Fortran label 40).
    // -----------------------------------------------------------------
    fn series_end(&mut self, s1r_in: f64, s1i_in: f64, s2r_in: f64, s2i_in: f64) {
        let s1r = s1r_in;
        let s1i = s1i_in;
        let s2r = s2r_in;
        let s2i = s2i_in;
        let mut str;
        let mut sti;

        if self.id != 1 {
            self.air = s1r * C1 - C2 * (self.zr * s2r - self.zi * s2i);
            self.aii = s1i * C1 - C2 * (self.zr * s2i + self.zi * s2r);
            if self.kode == 1 {
                return;
            }
            let (mut csqr, mut csqi) = (0.0f64, 0.0f64);
            zsqrt(self.zr, self.zi, &mut csqr, &mut csqi);
            let ztar = TTH * (self.zr * csqr - self.zi * csqi);
            let ztai = TTH * (self.zr * csqi + self.zi * csqr);
            let (mut er, mut ei) = (0.0f64, 0.0f64);
            zexp(ztar, ztai, &mut er, &mut ei);
            let ptr = self.air * er - self.aii * ei;
            self.aii = self.air * ei + self.aii * er;
            self.air = ptr;
            return;
        }

        self.air = -s2r * C2;
        self.aii = -s2i * C2;
        if self.az > self.tol {
            str = self.zr * s1r - self.zi * s1i;
            sti = self.zr * s1i + self.zi * s1r;
            let cc = C1 / (1.0 + self.fid);
            self.air += cc * (str * self.zr - sti * self.zi);
            self.aii += cc * (str * self.zi + sti * self.zr);
        }

        if self.kode == 1 {
            return;
        }
        let (mut csqr, mut csqi) = (0.0f64, 0.0f64);
        zsqrt(self.zr, self.zi, &mut csqr, &mut csqi);
        let ztar = TTH * (self.zr * csqr - self.zi * csqi);
        let ztai = TTH * (self.zr * csqi + self.zi * csqr);
        let (mut er, mut ei) = (0.0f64, 0.0f64);
        zexp(ztar, ztai, &mut er, &mut ei);
        let ptr = er * self.air - ei * self.aii;
        self.aii = er * self.aii + ei * self.air;
        self.air = ptr;
    }

    // -----------------------------------------------------------------
    // Case for |Z| > 1.0.
    // -----------------------------------------------------------------
    fn large_z(&mut self) {
        let fnu = (1.0 + self.fid) / 3.0;

        // Machine-constant setup.
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
        let rl = 1.2 * dig + 3.0;
        let alaz = self.az.ln();

        // Test for proper range.
        let mut aa = 0.5 / self.tol;
        let bb = (i1mach(9) as f64) * 0.5;
        aa = aa.min(bb);
        aa = aa.powf(TTH);
        if self.az > aa {
            self.ierr = 4;
            self.nz = 0;
            return;
        }
        aa = aa.sqrt();
        if self.az > aa {
            self.ierr = 3;
        }

        let (mut csqr, mut csqi) = (0.0f64, 0.0f64);
        zsqrt(self.zr, self.zi, &mut csqr, &mut csqi);
        let mut ztar = TTH * (self.zr * csqr - self.zi * csqi);
        let mut ztai = TTH * (self.zr * csqi + self.zi * csqr);

        self.csqr = csqr;
        self.csqi = csqi;

        // Re(ZTA) <= 0 when Re(Z) < 0, especially when Im(Z) is small.
        self.iflag = 0;
        self.sfac = 1.0;
        let ak = ztai;
        if self.zr < 0.0 {
            let bk = ztar;
            let ck = -bk.abs();
            ztar = ck;
            ztai = ak;
        }

        if self.zi == 0.0 && self.zr <= 0.0 {
            ztar = 0.0;
            ztai = ak;
        }

        aa = ztar;

        if aa >= 0.0 && self.zr > 0.0 {
            // Underflow test branch (Fortran label 110).
            self.large_underflow(aa, alaz, alim, elim, ztar, ztai, fnu, rl);
            return;
        }

        if self.kode != 2 {
            // Overflow test.
            if aa <= -alim {
                aa = -aa + 0.25 * alaz;
                self.iflag = 1;
                self.sfac = self.tol;
                if aa > elim {
                    self.nz = 0;
                    self.ierr = 2;
                    return;
                }
            }
        }

        // ZACAI returns EXP(ZTA)*K(FNU,ZTA) on KODE=2.
        let mr = if self.zi < 0.0 { -1 } else { 1 };
        let mut cyr = [0.0f64; 1];
        let mut cyi = [0.0f64; 1];
        let mut nn: i32 = 0;
        zacai(
            ztar, ztai, fnu, self.kode, mr, 1, &mut cyr, &mut cyi, &mut nn, rl, self.tol,
            elim, alim,
        );
        if nn < 0 {
            if nn == -1 {
                self.nz = 0;
                self.ierr = 2;
            } else {
                self.nz = 0;
                self.ierr = 5;
            }
            return;
        }
        self.nz += nn;

        self.output(&cyr, &cyi);
    }

    // -----------------------------------------------------------------
    // Underflow-test branch for large |Z| (Fortran label 110).
    // -----------------------------------------------------------------
    #[allow(clippy::too_many_arguments)]
    fn large_underflow(
        &mut self,
        aa_in: f64,
        alaz: f64,
        alim: f64,
        elim: f64,
        ztar: f64,
        ztai: f64,
        fnu: f64,
        rl: f64,
    ) {
        let mut aa = aa_in;

        if self.kode != 2 {
            if aa >= alim {
                aa = -aa - 0.25 * alaz;
                self.iflag = 2;
                self.sfac = 1.0 / self.tol;
                if aa < -elim {
                    self.nz = 1;
                    self.air = ZEROR;
                    self.aii = ZEROI;
                    return;
                }
            }
        }

        let mr = if ztai < 0.0 { -1 } else { 1 };
        let mut cyr = [0.0f64; 1];
        let mut cyi = [0.0f64; 1];
        let mut nn: i32 = 0;
        zacai(
            ztar, ztai, fnu, self.kode, mr, 1, &mut cyr, &mut cyi, &mut nn, rl, self.tol,
            elim, alim,
        );
        if nn < 0 {
            if nn == -1 {
                self.nz = 0;
                self.ierr = 2;
            } else {
                self.nz = 0;
                self.ierr = 5;
            }
            return;
        }
        self.nz = nn;

        self.output(&cyr, &cyi);
    }

    // -----------------------------------------------------------------
    // Output conversion (Fortran labels 130-160).
    // -----------------------------------------------------------------
    fn output(&mut self, cyr: &[f64; 1], cyi: &[f64; 1]) {
        let mut s1r = cyr[0] * COEF;
        let mut s1i = cyi[0] * COEF;

        if self.iflag == 0 {
            if self.id != 1 {
                self.air = self.csqr * s1r - self.csqi * s1i;
                self.aii = self.csqr * s1i + self.csqi * s1r;
                return;
            }
            self.air = -(self.zr * s1r - self.zi * s1i);
            self.aii = -(self.zr * s1i + self.zi * s1r);
            return;
        }

        s1r *= self.sfac;
        s1i *= self.sfac;

        if self.id != 1 {
            let str = s1r * self.csqr - s1i * self.csqi;
            s1i = s1r * self.csqi + s1i * self.csqr;
            s1r = str;
            self.air = s1r / self.sfac;
            self.aii = s1i / self.sfac;
            return;
        }

        let str = -(s1r * self.zr - s1i * self.zi);
        s1i = -(s1r * self.zi + s1i * self.zr);
        s1r = str;
        self.air = s1r / self.sfac;
        self.aii = s1i / self.sfac;
    }
}