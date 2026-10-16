use crate::amos::d1mach::d1mach;
use crate::amos::dgamln::dgamln;
use crate::amos::i1mach::i1mach;
use crate::amos::zabs::zabs;
use crate::amos::zdiv::zdiv;
use crate::amos::zexp::zexp;
use crate::amos::zkscl::zkscl;
use crate::amos::zlog::zlog;
use crate::amos::zmlt::zmlt;
use crate::amos::zshch::zshch;
use crate::amos::zsqrt::zsqrt;
use crate::amos::zuchk::zuchk;

// ---------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------

const CZEROR: f64 = 0.0;
const CZEROI: f64 = 0.0;
const CONER: f64 = 1.0;
const CONEI: f64 = 0.0;
const CTWOR: f64 = 2.0;
const R1: f64 = 2.0;
const DPI: f64 = 3.14159265358979324;
const RTHPI: f64 = 1.25331413731550025;
const SPI: f64 = 1.90985931710274403;
const HPI: f64 = 1.57079632679489662;
const FPI: f64 = 1.89769999331517738;
const TTH: f64 = 6.66666666666666666e-01;
const KMAX: i32 = 30;

const CC: [f64; 8] = [
    5.77215664901532861e-01,
    -4.20026350340952355e-02,
    -4.21977345555443367e-02,
    7.21894324666309954e-03,
    -2.15241674114950973e-04,
    -2.01348547807882387e-05,
    1.13302723198169588e-06,
    6.11609510448141582e-09,
];

// ---------------------------------------------------------------------
// State struct shared across the internal helpers.
//
// This carries every mutable local that in the Fortran survives across
// GO TO jumps between labels. The helpers take &mut self and read/write
// fields directly, matching the original host-association semantics.
// ---------------------------------------------------------------------

struct BknuState<'a> {
    // Parameters (read-only)
    zr: f64,
    zi: f64,
    fnu: f64,
    n: i32,
    tol: f64,
    elim: f64,
    alim: f64,

    // Outputs
    yr: &'a mut [f64],
    yi: &'a mut [f64],
    nz: i32,

    // Control flags / indices
    iflag: i32,
    koded: i32,
    inu: i32,
    kflag: i32,
    inub: i32,
    ic: i32,
    j: i32,
    kk: i32,

    // Persistent recurrence state (Fortran locals)
    dnu: f64,
    dnu2: f64,
    caz: f64,
    rzr: f64,
    rzi: f64,
    smur: f64,
    smui: f64,
    coefr: f64,
    coefi: f64,
    s1r: f64,
    s1i: f64,
    s2r: f64,
    s2i: f64,
    ckr: f64,
    cki: f64,
    fr: f64,
    fi: f64,
    pr: f64,
    pi: f64,
    qr: f64,
    qi: f64,
    ak: f64,
    a1: f64,
    bk: f64,
    zdr: f64,
    zdi: f64,
    cyr: [f64; 2],
    cyi: [f64; 2],
    csrr: [f64; 3],
    cssr: [f64; 3],
    bry: [f64; 3],
}

/// ZBKNU computes the K Bessel function in the right half Z plane.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zbknu(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    n: i32,
    yr: &mut [f64],
    yi: &mut [f64],
    nz: &mut i32,
    tol: f64,
    elim: f64,
    alim: f64,
) {
    let n_usize = n as usize;
    assert!(yr.len() >= n_usize, "yr must have at least n elements");
    assert!(yi.len() >= n_usize, "yi must have at least n elements");

    let caz = zabs(zr, zi);
    let csclr = 1.0 / tol;
    let crscr = tol;
    let cssr = [csclr, 1.0, crscr];
    let csrr = [crscr, 1.0, csclr];
    let bry_inner = 1.0e+3 * d1mach(1) / tol;
    let bry = [bry_inner, 1.0 / bry_inner, d1mach(2)];

    let rcaz = 1.0 / caz;
    let str = zr * rcaz;
    let sti = -zi * rcaz;
    let rzr = (str + str) * rcaz;
    let rzi = (sti + sti) * rcaz;

    let inu = (fnu + 0.5) as i32;
    let dnu = fnu - (inu as f64);

    let mut state = BknuState {
        zr,
        zi,
        fnu,
        n,
        tol,
        elim,
        alim,
        yr,
        yi,
        nz: 0,
        iflag: 0,
        koded: kode,
        inu,
        kflag: 2,
        inub: 0,
        ic: 0,
        j: 0,
        kk: 0,
        dnu,
        dnu2: 0.0,
        caz,
        rzr,
        rzi,
        smur: 0.0,
        smui: 0.0,
        coefr: 0.0,
        coefi: 0.0,
        s1r: 0.0,
        s1i: 0.0,
        s2r: 0.0,
        s2i: 0.0,
        ckr: 0.0,
        cki: 0.0,
        fr: 0.0,
        fi: 0.0,
        pr: 0.0,
        pi: 0.0,
        qr: 0.0,
        qi: 0.0,
        ak: 0.0,
        a1: 0.0,
        bk: 0.0,
        zdr: 0.0,
        zdi: 0.0,
        cyr: [0.0; 2],
        cyi: [0.0; 2],
        csrr,
        cssr,
        bry,
    };

    state.run();
    *nz = state.nz;
}

impl<'a> BknuState<'a> {
    // -----------------------------------------------------------------
    // Main driver: dispatch based on DNU and CAZ.
    // -----------------------------------------------------------------
    fn run(&mut self) {
        // DNU = +/- 0.5 case.
        if self.dnu.abs() == 0.5 {
            self.dnu_half();
            return;
        }

        self.dnu2 = if self.dnu.abs() > self.tol {
            self.dnu * self.dnu
        } else {
            0.0
        };

        // CAZ > R1 case.
        if self.caz > R1 {
            self.large_abs();
            return;
        }

        // -------------------------------------------------------------
        // Series for |Z| <= R1
        // -------------------------------------------------------------
        let mut fc = 1.0;
        let (mut smur, mut smui) = (0.0f64, 0.0f64);
        let mut idum: i32 = 0;
        let _ = zlog(self.rzr, self.rzi, &mut smur, &mut smui);
        let fmur = smur * self.dnu;
        let fmui = smui * self.dnu;
        let (mut cshr, mut cshi, mut cchr, mut cchi) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        zshch(fmur, fmui, &mut cshr, &mut cshi, &mut cchr, &mut cchi);
        if self.dnu != 0.0 {
            fc = self.dnu * DPI;
            fc /= fc.sin();
            smur = cshr / self.dnu;
            smui = cshi / self.dnu;
        }
        self.smur = smur;
        self.smui = smui;

        let a2 = 1.0 + self.dnu;

        // GAM(1-Z)*GAM(1+Z) = PI*Z/SIN(PI*Z)
        let t2 = (-dgamln(a2, &mut idum)).exp();
        let t1 = 1.0 / (t2 * fc);

        let g1;
        if self.dnu.abs() > 0.1 {
            g1 = (t1 - t2) / (self.dnu + self.dnu);
        } else {
            let mut ak = 1.0;
            let mut sum = CC[0];
            for k in 1..8 {
                ak *= self.dnu2;
                let tm = CC[k] * ak;
                sum += tm;
                if tm.abs() < self.tol {
                    break;
                }
            }
            g1 = -sum;
        }

        let g2 = (t1 + t2) * 0.5;
        let mut fr = fc * (cchr * g1 + self.smur * g2);
        let mut fi = fc * (cchi * g1 + self.smui * g2);
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        zexp(fmur, fmui, &mut str, &mut sti);
        let mut pr = 0.5 * str / t2;
        let mut pi_ = 0.5 * sti / t2;
        let (mut ptr, mut pti) = (0.0f64, 0.0f64);
        zdiv(0.5, 0.0, str, sti, &mut ptr, &mut pti);
        let mut qr = ptr / t1;
        let mut qi = pti / t1;

        self.s1r = fr;
        self.s1i = fi;
        self.s2r = pr;
        self.s2i = pi_;

        let mut ak = 1.0;
        let mut a1 = 1.0;
        let mut ckr = CONER;
        let mut cki = CONEI;
        let mut bk = 1.0 - self.dnu2;

        // -------------------------------------------------------------
        // Generate K(FNU,Z), 0 <= FNU < 0.5 and N = 1
        // -------------------------------------------------------------
        if self.inu <= 0 && self.n <= 1 {
            if self.caz >= self.tol {
                let (mut czr, mut czi) = (0.0f64, 0.0f64);
                zmlt(self.zr, self.zi, self.zr, self.zi, &mut czr, &mut czi);
                czr *= 0.25;
                czi *= 0.25;
                let t1 = 0.25 * self.caz * self.caz;
                loop {
                    fr = (fr * ak + pr + qr) / bk;
                    fi = (fi * ak + pi_ + qi) / bk;
                    str = 1.0 / (ak - self.dnu);
                    pr *= str;
                    pi_ *= str;
                    str = 1.0 / (ak + self.dnu);
                    qr *= str;
                    qi *= str;
                    str = ckr * czr - cki * czi;
                    let rak = 1.0 / ak;
                    cki = (ckr * czi + cki * czr) * rak;
                    ckr = str * rak;
                    self.s1r = ckr * fr - cki * fi + self.s1r;
                    self.s1i = ckr * fi + cki * fr + self.s1i;
                    a1 = a1 * t1 * rak;
                    bk = bk + ak + ak + 1.0;
                    ak += 1.0;
                    if a1 <= self.tol {
                        break;
                    }
                }
            }

            self.yr[0] = self.s1r;
            self.yi[0] = self.s1i;
            if self.koded == 1 {
                return;
            }
            zexp(self.zr, self.zi, &mut str, &mut sti);
            let (mut outr, mut outi) = (0.0f64, 0.0f64);
            zmlt(self.s1r, self.s1i, str, sti, &mut outr, &mut outi);
            self.yr[0] = outr;
            self.yi[0] = outi;
            return;
        }

        // -------------------------------------------------------------
        // Generate K(DNU,Z) and K(DNU+1,Z) for forward recurrence
        // -------------------------------------------------------------
        if self.caz >= self.tol {
            let (mut czr, mut czi) = (0.0f64, 0.0f64);
            zmlt(self.zr, self.zi, self.zr, self.zi, &mut czr, &mut czi);
            czr *= 0.25;
            czi *= 0.25;
            let t1 = 0.25 * self.caz * self.caz;
            loop {
                fr = (fr * ak + pr + qr) / bk;
                fi = (fi * ak + pi_ + qi) / bk;
                str = 1.0 / (ak - self.dnu);
                pr *= str;
                pi_ *= str;
                str = 1.0 / (ak + self.dnu);
                qr *= str;
                qi *= str;
                str = ckr * czr - cki * czi;
                let rak = 1.0 / ak;
                cki = (ckr * czi + cki * czr) * rak;
                ckr = str * rak;
                self.s1r = ckr * fr - cki * fi + self.s1r;
                self.s1i = ckr * fi + cki * fr + self.s1i;
                str = pr - fr * ak;
                sti = pi_ - fi * ak;
                self.s2r = ckr * str - cki * sti + self.s2r;
                self.s2i = ckr * sti + cki * str + self.s2i;
                a1 = a1 * t1 * rak;
                bk = bk + ak + ak + 1.0;
                ak += 1.0;
                if a1 <= self.tol {
                    break;
                }
            }
        }

        self.kflag = 2;
        a1 = self.fnu + 1.0;
        ak = a1 * self.smur.abs();
        if ak > self.alim {
            self.kflag = 3;
        }
        str = self.cssr[(self.kflag - 1) as usize];
        let p2r = self.s2r * str;
        let p2i = self.s2i * str;
        zmlt(p2r, p2i, self.rzr, self.rzi, &mut self.s2r, &mut self.s2i);
        self.s1r *= str;
        self.s1i *= str;
        if self.koded != 1 {
            let (mut fr2, mut fi2) = (0.0f64, 0.0f64);
            zexp(self.zr, self.zi, &mut fr2, &mut fi2);
            let (mut a, mut b) = (0.0f64, 0.0f64);
            zmlt(self.s1r, self.s1i, fr2, fi2, &mut a, &mut b);
            self.s1r = a;
            self.s1i = b;
            zmlt(self.s2r, self.s2i, fr2, fi2, &mut a, &mut b);
            self.s2r = a;
            self.s2i = b;
        }
        self.forward_210();
    }

    // -----------------------------------------------------------------
    // DNU = +/- 0.5 branch (Fortran label 110 / 300).
    // -----------------------------------------------------------------
    fn dnu_half(&mut self) {
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        zsqrt(self.zr, self.zi, &mut str, &mut sti);
        zdiv(RTHPI, CZEROI, str, sti, &mut self.coefr, &mut self.coefi);
        self.kflag = 2;

        if self.koded != 2 {
            if self.zr > self.alim {
                // Scale by DEXP(Z), IFLAG = 1.
                self.koded = 2;
                self.iflag = 1;
                self.kflag = 2;
                self.check_half();
                return;
            }
            str = (-self.zr).exp() * self.cssr[(self.kflag - 1) as usize];
            sti = -str * self.zi.sin();
            str *= self.zi.cos();
            let (mut a, mut b) = (0.0f64, 0.0f64);
            zmlt(self.coefr, self.coefi, str, sti, &mut a, &mut b);
            self.coefr = a;
            self.coefi = b;
        }

        self.check_half();
    }

    // -----------------------------------------------------------------
    // |Z| > R1 branch (Fortran label 110 onwards).
    // -----------------------------------------------------------------
    fn large_abs(&mut self) {
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        zsqrt(self.zr, self.zi, &mut str, &mut sti);
        zdiv(RTHPI, CZEROI, str, sti, &mut self.coefr, &mut self.coefi);
        self.kflag = 2;

        if self.koded != 2 {
            if self.zr > self.alim {
                self.koded = 2;
                self.iflag = 1;
                self.kflag = 2;
                self.check_half();
                return;
            }
            str = (-self.zr).exp() * self.cssr[(self.kflag - 1) as usize];
            sti = -str * self.zi.sin();
            str *= self.zi.cos();
            let (mut a, mut b) = (0.0f64, 0.0f64);
            zmlt(self.coefr, self.coefi, str, sti, &mut a, &mut b);
            self.coefr = a;
            self.coefi = b;
        }

        self.check_half();
    }

    fn check_half(&mut self) {
        if self.dnu.abs() == 0.5 {
            // Fortran label 300.
            self.s1r = self.coefr;
            self.s1i = self.coefi;
            self.s2r = self.coefr;
            self.s2i = self.coefi;
            self.forward_210();
            return;
        }
        self.miller();
    }

    // -----------------------------------------------------------------
    // Miller algorithm for |Z| > R1.
    // -----------------------------------------------------------------
    fn miller(&mut self) {
        let mut ak = (DPI * self.dnu).cos().abs();
        if ak == CZEROR {
            self.s1r = self.coefr;
            self.s1i = self.coefi;
            self.s2r = self.coefr;
            self.s2i = self.coefi;
            self.forward_210();
            return;
        }

        let mut fhs = (0.25 - self.dnu2).abs();
        if fhs == CZEROR {
            self.s1r = self.coefr;
            self.s1i = self.coefi;
            self.s2r = self.coefr;
            self.s2i = self.coefi;
            self.forward_210();
            return;
        }

        // Compute R2 = F(E).
        let mut t1 = (i1mach(14) - 1) as f64;
        t1 = t1 * d1mach(5) * 3.321928094;
        t1 = t1.max(12.0).min(60.0);
        let t2 = TTH * t1 - 6.0;

        if self.zr != 0.0 {
            t1 = (self.zi / self.zr).atan().abs();
        } else {
            t1 = HPI;
        }

        let mut fk;
        if t2 <= self.caz {
            // Forward recurrence to determine backward index K.
            let etest = ak / (DPI * self.caz * self.tol);
            fk = CONER;
            if etest >= CONER {
                let mut fks = CTWOR;
                let mut ckr = self.caz + self.caz + CTWOR;
                let mut p1r = CZEROR;
                let mut p2r = CONER;
                let mut found = false;
                let mut i = 0;
                for ii in 1..=KMAX {
                    i = ii;
                    ak = fhs / fks;
                    let cbr = ckr / (fk + CONER);
                    let ptr = p2r;
                    p2r = cbr * p2r - p1r * ak;
                    p1r = ptr;
                    ckr += CTWOR;
                    fks = fks + fk + fk + CTWOR;
                    fhs = fhs + fk + fk;
                    fk += CONER;
                    let str = p2r.abs() * fk;
                    if etest < str {
                        found = true;
                        break;
                    }
                }
                let _ = i;
                if !found {
                    // Fortran label 310: NZ = -2.
                    self.nz = -2;
                    return;
                }
                fk = fk + SPI * t1 * (t2 / self.caz).sqrt();
                fhs = (0.25 - self.dnu2).abs();
            }
        } else {
            // Compute backward index K for |Z| < R2.
            let a2 = self.caz.sqrt();
            ak = FPI * ak / (self.tol * a2.sqrt());
            let aa = 3.0 * t1 / (1.0 + self.caz);
            let bb = 14.7 * t1 / (28.0 + self.caz);
            ak = (ak.ln() + self.caz * aa.cos() / (1.0 + 0.008 * self.caz)) / bb.cos();
            fk = 0.12125 * ak * ak / self.caz + 1.5;
        }

        // -------------------------------------------------------------
        // Backward recurrence loop for Miller algorithm
        // -------------------------------------------------------------
        let k = fk as i32;
        fk = k as f64;
        let mut fks = fk * fk;
        let mut p1r = CZEROR;
        let mut p1i = CZEROI;
        let mut p2r = self.tol;
        let mut p2i = CZEROI;
        let mut csr = p2r;
        let mut csi = p2i;
        for _ in 1..=k {
            let a1 = fks - fk;
            ak = (fks + fk) / (a1 + fhs);
            let rak = 2.0 / (fk + CONER);
            let cbr = (fk + self.zr) * rak;
            let cbi = self.zi * rak;
            let ptr = p2r;
            let pti = p2i;
            p2r = (ptr * cbr - pti * cbi - p1r) * ak;
            p2i = (pti * cbr + ptr * cbi - p1i) * ak;
            p1r = ptr;
            p1i = pti;
            csr += p2r;
            csi += p2i;
            fks = a1 - fk + CONER;
            fk -= CONER;
        }

        // Compute (P2/CS) = (P2/|CS|)*(CONJG(CS)/|CS|)
        let tm = zabs(csr, csi);
        let ptr = 1.0 / tm;
        self.s1r = p2r * ptr;
        self.s1i = p2i * ptr;
        csr *= ptr;
        csi = -csi * ptr;
        let (mut str, mut sti) = (0.0f64, 0.0f64);
        zmlt(self.coefr, self.coefi, self.s1r, self.s1i, &mut str, &mut sti);
        zmlt(str, sti, csr, csi, &mut self.s1r, &mut self.s1i);

        if self.inu <= 0 && self.n <= 1 {
            self.zdr = self.zr;
            self.zdi = self.zi;
            if self.iflag == 1 {
                self.iflag1_270();
                return;
            }
            self.output_240();
            return;
        }

        // Compute P1/P2 = (P1/|P2|)*(CONJG(P2)/|P2|)
        let tm = zabs(p2r, p2i);
        let ptr = 1.0 / tm;
        p1r *= ptr;
        p1i *= ptr;
        p2r *= ptr;
        p2i = -p2i * ptr;
        let (mut ptr2, mut pti2) = (0.0f64, 0.0f64);
        zmlt(p1r, p1i, p2r, p2i, &mut ptr2, &mut pti2);
        let str2 = self.dnu + 0.5 - ptr2;
        let sti2 = -pti2;
        let (mut str3, mut sti3) = (0.0f64, 0.0f64);
        zdiv(str2, sti2, self.zr, self.zi, &mut str3, &mut sti3);
        let str4 = str3 + 1.0;
        zmlt(str4, sti3, self.s1r, self.s1i, &mut self.s2r, &mut self.s2i);

        self.forward_210();
    }

    // -----------------------------------------------------------------
    // Forward recursion (Fortran label 210).
    // -----------------------------------------------------------------
    fn forward_210(&mut self) {
        let str = self.dnu + 1.0;
        self.ckr = str * self.rzr;
        self.cki = str * self.rzi;

        if self.n == 1 {
            self.inu -= 1;
        }

        if self.inu <= 0 {
            if self.n <= 1 {
                self.s1r = self.s2r;
                self.s1i = self.s2i;
            }
            self.zdr = self.zr;
            self.zdi = self.zi;
            if self.iflag == 1 {
                self.iflag1_270();
                return;
            }
            self.output_240();
            return;
        }

        self.inub = 1;
        if self.iflag == 1 {
            self.iflag1_261();
            return;
        }

        self.forward_225();
    }

    // -----------------------------------------------------------------
    // Scaled forward recurrence (Fortran label 225).
    // -----------------------------------------------------------------
    fn forward_225(&mut self) {
        let mut p1r = self.csrr[(self.kflag - 1) as usize];
        let mut ascle = self.bry[(self.kflag - 1) as usize];

        let mut i = self.inub;
        while i <= self.inu {
            let str = self.s2r;
            let sti = self.s2i;
            self.s2r = self.ckr * str - self.cki * sti + self.s1r;
            self.s2i = self.cki * str + self.ckr * sti + self.s1i;
            self.s1r = str;
            self.s1i = sti;
            self.ckr += self.rzr;
            self.cki += self.rzi;

            if self.kflag < 3 {
                let p2r = self.s2r * p1r;
                let p2i = self.s2i * p1r;
                let p2m = p2r.abs().max(p2i.abs());
                if p2m > ascle {
                    self.kflag += 1;
                    ascle = self.bry[(self.kflag - 1) as usize];
                    self.s1r *= p1r;
                    self.s1i *= p1r;
                    self.s2r = p2r;
                    self.s2i = p2i;
                    let s = self.cssr[(self.kflag - 1) as usize];
                    self.s1r *= s;
                    self.s1i *= s;
                    self.s2r *= s;
                    self.s2i *= s;
                    p1r = self.csrr[(self.kflag - 1) as usize];
                }
            }
            i += 1;
        }

        if self.n == 1 {
            self.s1r = self.s2r;
            self.s1i = self.s2i;
        }

        self.output_240();
    }

    // -----------------------------------------------------------------
    // Output section (Fortran label 240).
    // -----------------------------------------------------------------
    fn output_240(&mut self) {
        let str = self.csrr[(self.kflag - 1) as usize];
        self.yr[0] = self.s1r * str;
        self.yi[0] = self.s1i * str;
        if self.n == 1 {
            return;
        }
        self.yr[1] = self.s2r * str;
        self.yi[1] = self.s2i * str;
        if self.n == 2 {
            return;
        }
        self.kk = 2;
        self.output_250();
    }

    fn output_250(&mut self) {
        loop {
            self.kk += 1;
            if self.kk > self.n {
                return;
            }
            let mut p1r = self.csrr[(self.kflag - 1) as usize];
            let mut ascle = self.bry[(self.kflag - 1) as usize];

            let mut i = self.kk;
            while i <= self.n {
                let p2r = self.s2r;
                let p2i = self.s2i;
                self.s2r = self.ckr * p2r - self.cki * p2i + self.s1r;
                self.s2i = self.cki * p2r + self.ckr * p2i + self.s1i;
                self.s1r = p2r;
                self.s1i = p2i;
                self.ckr += self.rzr;
                self.cki += self.rzi;
                let yp2r = self.s2r * p1r;
                let yp2i = self.s2i * p1r;
                self.yr[(i - 1) as usize] = yp2r;
                self.yi[(i - 1) as usize] = yp2i;

                if self.kflag >= 3 {
                    i += 1;
                    continue;
                }
                let p2m = yp2r.abs().max(yp2i.abs());
                if p2m <= ascle {
                    i += 1;
                    continue;
                }

                // KFLAG bump — rescale and break out to the next KK.
                self.kflag += 1;
                ascle = self.bry[(self.kflag - 1) as usize];
                self.s1r *= p1r;
                self.s1i *= p1r;
                self.s2r = yp2r;
                self.s2i = yp2i;
                let s = self.cssr[(self.kflag - 1) as usize];
                self.s1r *= s;
                self.s1i *= s;
                self.s2r *= s;
                self.s2i *= s;
                p1r = self.csrr[(self.kflag - 1) as usize];
                let _ = ascle;
                i += 1;
                break;
            }
        }
    }

    // -----------------------------------------------------------------
    // IFLAG = 1 cases: forward recurrence on scaled values on underflow
    // (Fortran label 261).
    // -----------------------------------------------------------------
    fn iflag1_261(&mut self) {
        let helim = 0.5 * self.elim;
        let elm = (-self.elim).exp();
        let celmr = elm;
        let mut ascle = self.bry[0];
        self.zdr = self.zr;
        self.zdi = self.zi;
        self.ic = -1;
        self.j = 2;

        let mut i = 1;
        while i <= self.inu {
            let str = self.s2r;
            let sti = self.s2i;
            self.s2r = str * self.ckr - sti * self.cki + self.s1r;
            self.s2i = sti * self.ckr + str * self.cki + self.s1i;
            self.s1r = str;
            self.s1i = sti;
            self.ckr += self.rzr;
            self.cki += self.rzi;
            let as_ = zabs(self.s2r, self.s2i);
            let alas = as_.ln();
            let p2r = -self.zdr + alas;

            if p2r >= -self.elim {
                let (mut str2, mut sti2) = (0.0f64, 0.0f64);
                let mut idum: i32 = 0;
                let _ = zlog(self.s2r, self.s2i, &mut str2, &mut sti2);
                let p2r2 = -self.zdr + str2;
                let p2i2 = -self.zdi + sti2;
                let p2m = p2r2.exp() / self.tol;
                let p1r = p2m * p2i2.cos();
                let p1i = p2m * p2i2.sin();
                let mut nw: i32 = 0;
                zuchk(p1r, p1i, &mut nw, ascle, self.tol);
                if nw == 0 {
                    self.j = 3 - self.j;
                    self.cyr[(self.j - 1) as usize] = p1r;
                    self.cyi[(self.j - 1) as usize] = p1i;
                    if self.ic == i - 1 {
                        self.iflag1_264(i);
                        return;
                    }
                    self.ic = i;
                    i += 1;
                    let _ = ascle;
                    continue;
                }
            }

            if alas >= helim {
                self.zdr -= self.elim;
                self.s1r *= celmr;
                self.s1i *= celmr;
                self.s2r *= celmr;
                self.s2i *= celmr;
            }
            ascle = self.bry[0];
            i += 1;
        }

        if self.n != 1 {
            self.iflag1_270();
            return;
        }
        self.s1r = self.s2r;
        self.s1i = self.s2i;
        self.iflag1_270();
    }

    fn iflag1_264(&mut self, i: i32) {
        self.kflag = 1;
        self.inub = i + 1;
        self.s2r = self.cyr[(self.j - 1) as usize];
        self.s2i = self.cyi[(self.j - 1) as usize];
        self.j = 3 - self.j;
        self.s1r = self.cyr[(self.j - 1) as usize];
        self.s1i = self.cyi[(self.j - 1) as usize];
        if self.inub <= self.inu {
            self.forward_225();
            return;
        }
        if self.n != 1 {
            self.output_240();
            return;
        }
        self.s1r = self.s2r;
        self.s1i = self.s2i;
        self.output_240();
    }

    // -----------------------------------------------------------------
    // IFLAG = 1 output path (Fortran label 270).
    // -----------------------------------------------------------------
    fn iflag1_270(&mut self) {
        self.yr[0] = self.s1r;
        self.yi[0] = self.s1i;
        if self.n != 1 {
            self.yr[1] = self.s2r;
            self.yi[1] = self.s2i;
        }
        let ascle = self.bry[0];
        zkscl(
            self.zdr,
            self.zdi,
            self.fnu,
            self.n,
            self.yr,
            self.yi,
            &mut self.nz,
            self.rzr,
            self.rzi,
            ascle,
            self.tol,
            self.elim,
        );
        self.inu = self.n - self.nz;
        if self.inu <= 0 {
            return;
        }
        let mut kk = self.nz + 1;
        self.s1r = self.yr[(kk - 1) as usize];
        self.s1i = self.yi[(kk - 1) as usize];
        self.yr[(kk - 1) as usize] = self.s1r * self.csrr[0];
        self.yi[(kk - 1) as usize] = self.s1i * self.csrr[0];
        if self.inu == 1 {
            return;
        }
        kk = self.nz + 2;
        self.s2r = self.yr[(kk - 1) as usize];
        self.s2i = self.yi[(kk - 1) as usize];
        self.yr[(kk - 1) as usize] = self.s2r * self.csrr[0];
        self.yi[(kk - 1) as usize] = self.s2i * self.csrr[0];
        if self.inu == 2 {
            return;
        }
        let t2 = self.fnu + ((kk - 1) as f64);
        self.ckr = t2 * self.rzr;
        self.cki = t2 * self.rzi;
        self.kflag = 1;
        self.kk = kk - 1;
        self.output_250();
    }
}