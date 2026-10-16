use crate::amos::d1mach::d1mach;
use crate::amos::zabs::zabs;
use crate::amos::zairy::zairy;
use crate::amos::zs1s2::zs1s2;
use crate::amos::zuchk::zuchk;
use crate::amos::zunhj::zunhj;

// ---------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------

const ZEROR: f64 = 0.0;
const ZEROI: f64 = 0.0;
const CONER: f64 = 1.0;
const CR1R: f64 = 1.0;
const CR1I: f64 = 1.73205080756887729;
const CR2R: f64 = -0.5;
const CR2I: f64 = -8.66025403784438647e-01;
const HPI: f64 = 1.57079632679489662;
const PI: f64 = 3.14159265358979324;
const AIC: f64 = 1.26551212348464539;

const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
const CIPI: [f64; 4] = [0.0, -1.0, 0.0, 1.0];

// ---------------------------------------------------------------------
// Flow control.
// ---------------------------------------------------------------------

enum Flow {
    Return,
    Continue,
}

// ---------------------------------------------------------------------
// State struct.
// ---------------------------------------------------------------------

struct Zunk2State<'a> {
    // Inputs (read-only)
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    mr: i32,
    n: i32,
    tol: f64,
    elim: f64,
    alim: f64,

    // Outputs
    yr: &'a mut [f64],
    yi: &'a mut [f64],
    nz: i32,

    // Persistent control state
    kdflg: i32,
    kflag: i32,
    iflag: i32,
    ib: i32,
    ic: i32,
    kk: i32,
    j: i32,

    // Persistent recurrence state
    zrr: f64,
    zri: f64,
    fn_val: f64,
    ckr: f64,
    cki: f64,
    rzr: f64,
    rzi: f64,
    c2r: f64,
    c2i: f64,
    csr: f64,
    csi: f64,
    cspnr: f64,
    cspni: f64,
    csgni: f64,
    sgn: f64,
    fnf: f64,
    inu: i32,
    asc: f64,
    iuf: i32,

    // ZN and associated values
    znr: f64,
    zni: f64,
    zbr: f64,
    zbi: f64,
    yy: f64,
    car: f64,
    sar: f64,

    // Persistent arrays
    bry: [f64; 3],
    cssr: [f64; 3],
    csrr: [f64; 3],
    cyr: [f64; 2],
    cyi: [f64; 2],
    phir: [f64; 2],
    phii: [f64; 2],
    argr: [f64; 2],
    argi: [f64; 2],
    sumr: [f64; 2],
    sumi: [f64; 2],
    zeta1r: [f64; 2],
    zeta1i: [f64; 2],
    zeta2r: [f64; 2],
    zeta2i: [f64; 2],
    asumr: [f64; 2],
    asumi: [f64; 2],
    bsumr: [f64; 2],
    bsumi: [f64; 2],

    // Pending sort results (Fortran labels 172/175)
    pending_phidr: f64,
    pending_phidi: f64,
    pending_argdr: f64,
    pending_argdi: f64,
    pending_zet1dr: f64,
    pending_zet1di: f64,
    pending_zet2dr: f64,
    pending_zet2di: f64,
    pending_asumdr: f64,
    pending_asumdi: f64,
    pending_bsumdr: f64,
    pending_bsumdi: f64,
    pending_called_zunhj: bool,

    // The value of I at exit of the first loop (for IB = I + 1)
    last_i: i32,
}

/// ZUNK2 computes K(FNU,Z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansions for H(KIND,FNU,ZN) and J(FNU,ZN)
/// where ZN is in the right half plane, KIND=(3-MR)/2, MR=+1 or
/// -1. Here ZN=ZR*I or -ZR*I where ZR=Z if Z is in the right
/// half plane or ZR=-Z if Z is in the left half plane. MR indicates
/// the direction of rotation for analytic continuation.
/// NZ=-1 means an overflow will occur.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zunk2(
    zr: f64,
    zi: f64,
    fnu: f64,
    kode: i32,
    mr: i32,
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

    let mut state = Zunk2State {
        zr,
        zi,
        fnu,
        kode,
        mr,
        n,
        tol,
        elim,
        alim,
        yr,
        yi,
        nz: 0,
        kdflg: 1,
        kflag: 2,
        iflag: 0,
        ib: 0,
        ic: 0,
        kk: 0,
        j: 2,
        zrr: 0.0,
        zri: 0.0,
        fn_val: 0.0,
        ckr: 0.0,
        cki: 0.0,
        rzr: 0.0,
        rzi: 0.0,
        c2r: 0.0,
        c2i: 0.0,
        csr: 0.0,
        csi: 0.0,
        cspnr: 0.0,
        cspni: 0.0,
        csgni: 0.0,
        sgn: 0.0,
        fnf: 0.0,
        inu: 0,
        asc: 0.0,
        iuf: 0,
        znr: 0.0,
        zni: 0.0,
        zbr: 0.0,
        zbi: 0.0,
        yy: 0.0,
        car: 0.0,
        sar: 0.0,
        bry: [0.0; 3],
        cssr: [0.0; 3],
        csrr: [0.0; 3],
        cyr: [0.0; 2],
        cyi: [0.0; 2],
        phir: [0.0; 2],
        phii: [0.0; 2],
        argr: [0.0; 2],
        argi: [0.0; 2],
        sumr: [0.0; 2],
        sumi: [0.0; 2],
        zeta1r: [0.0; 2],
        zeta1i: [0.0; 2],
        zeta2r: [0.0; 2],
        zeta2i: [0.0; 2],
        asumr: [0.0; 2],
        asumi: [0.0; 2],
        bsumr: [0.0; 2],
        bsumi: [0.0; 2],
        pending_phidr: 0.0,
        pending_phidi: 0.0,
        pending_argdr: 0.0,
        pending_argdi: 0.0,
        pending_zet1dr: 0.0,
        pending_zet1di: 0.0,
        pending_zet2dr: 0.0,
        pending_zet2di: 0.0,
        pending_asumdr: 0.0,
        pending_asumdi: 0.0,
        pending_bsumdr: 0.0,
        pending_bsumdi: 0.0,
        pending_called_zunhj: false,
        last_i: 0,
    };

    state.run();
    *nz = state.nz;
}

impl<'a> Zunk2State<'a> {
    fn run(&mut self) {
        self.kdflg = 1;
        self.nz = 0;

        // exp(-ALIM) = exp(-ELIM)/TOL
        let cscl = 1.0 / self.tol;
        let crsc = self.tol;
        self.cssr = [cscl, CONER, crsc];
        self.csrr = [crsc, CONER, cscl];
        self.bry[0] = 1.0e+3 * d1mach(1) / self.tol;
        self.bry[1] = 1.0 / self.bry[0];
        self.bry[2] = d1mach(2);

        self.zrr = self.zr;
        self.zri = self.zi;
        if self.zr < 0.0 {
            self.zrr = -self.zr;
            self.zri = -self.zi;
        }

        self.yy = self.zri;
        self.znr = self.zri;
        self.zni = -self.zrr;
        self.zbr = self.zrr;
        self.zbi = self.zri;
        self.inu = self.fnu as i32;
        self.fnf = self.fnu - (self.inu as f64);
        let ang = -HPI * self.fnf;
        self.car = ang.cos();
        self.sar = ang.sin();
        let c2r_0 = HPI * self.sar;
        let c2i_0 = -HPI * self.car;
        let kk_idx = (self.inu % 4 + 4) % 4 + 1;
        let str = c2r_0 * CIPR[(kk_idx - 1) as usize] - c2i_0 * CIPI[(kk_idx - 1) as usize];
        let sti = c2r_0 * CIPI[(kk_idx - 1) as usize] + c2i_0 * CIPR[(kk_idx - 1) as usize];
        self.csr = CR1R * str - CR1I * sti;
        self.csi = CR1R * sti + CR1I * str;

        if self.yy <= 0.0 {
            self.znr = -self.znr;
            self.zbi = -self.zbi;
        }

        self.j = 2;

        // -------------------------------------------------------------
        // First pass: compute K(FNU+I-1,Z) for I=1,N.
        // -------------------------------------------------------------
        if let Flow::Return = self.first_loop() {
            return;
        }

        // -------------------------------------------------------------
        // Forward recurrence setup
        // -------------------------------------------------------------
        let razr = 1.0 / zabs(self.zrr, self.zri);
        let str = self.zrr * razr;
        let sti = -self.zri * razr;
        self.rzr = (str + str) * razr;
        self.rzi = (sti + sti) * razr;
        self.ckr = self.fn_val * self.rzr;
        self.cki = self.fn_val * self.rzi;
        self.ib = self.last_i + 1;

        // -------------------------------------------------------------
        // Test last member for underflow and overflow.
        // -------------------------------------------------------------
        if self.n >= self.ib {
            if let Flow::Return = self.last_member_check() {
                return;
            }
            self.forward_recur();
        }

        if self.mr == 0 {
            return;
        }

        // -------------------------------------------------------------
        // Analytic continuation for Re(Z) < 0
        // -------------------------------------------------------------
        self.nz = 0;
        let fmr = self.mr as f64;
        self.sgn = -PI.copysign(fmr);

        self.csgni = self.sgn;
        if self.yy <= 0.0 {
            self.csgni = -self.csgni;
        }
        let ifn = self.inu + self.n - 1;
        let ang = self.fnf * self.sgn;
        self.cspnr = ang.cos();
        self.cspni = ang.sin();
        if ifn % 2 != 0 {
            self.cspnr = -self.cspnr;
            self.cspni = -self.cspni;
        }

        // CS = coefficient of the J function to get the I function.
        self.csr = self.sar * self.csgni;
        self.csi = self.car * self.csgni;
        let in_ = (ifn % 4 + 4) % 4 + 1;
        let c2r = CIPR[(in_ - 1) as usize];
        let c2i = CIPI[(in_ - 1) as usize];
        let str = self.csr * c2r + self.csi * c2i;
        self.csi = -self.csr * c2i + self.csi * c2r;
        self.csr = str;

        self.asc = self.bry[0];
        self.iuf = 0;
        self.kk = self.n;
        self.kdflg = 1;
        self.ib -= 1;
        self.ic = self.ib - 1;

        if let Flow::Return = self.continuation_loop() {
            return;
        }

        // -------------------------------------------------------------
        // Recur backward for remainder of I sequence and add in the
        // K functions.
        // -------------------------------------------------------------
        let il = self.n - self.kk;
        if il == 0 {
            return;
        }

        self.backward_recur(il);
    }

    // -----------------------------------------------------------------
    // First pass loop (Fortran DO 80 I=1,N).
    // -----------------------------------------------------------------
    fn first_loop(&mut self) -> Flow {
        let n = self.n;

        for i in 1..=n {
            self.j = 3 - self.j;
            let fn_val = self.fnu + ((i - 1) as f64);
            let j = self.j as usize;

            let (mut phir, mut phii) = (0.0f64, 0.0f64);
            let (mut argr, mut argi) = (0.0f64, 0.0f64);
            let (mut zeta1r, mut zeta1i) = (0.0f64, 0.0f64);
            let (mut zeta2r, mut zeta2i) = (0.0f64, 0.0f64);
            let (mut asumr, mut asumi) = (0.0f64, 0.0f64);
            let (mut bsumr, mut bsumi) = (0.0f64, 0.0f64);

            zunhj(
                self.znr, self.zni, fn_val, 0, self.tol, &mut phir, &mut phii, &mut argr,
                &mut argi, &mut zeta1r, &mut zeta1i, &mut zeta2r, &mut zeta2i, &mut asumr,
                &mut asumi, &mut bsumr, &mut bsumi,
            );

            self.phir[j - 1] = phir;
            self.phii[j - 1] = phii;
            self.argr[j - 1] = argr;
            self.argi[j - 1] = argi;
            self.zeta1r[j - 1] = zeta1r;
            self.zeta1i[j - 1] = zeta1i;
            self.zeta2r[j - 1] = zeta2r;
            self.zeta2i[j - 1] = zeta2i;
            self.asumr[j - 1] = asumr;
            self.asumi[j - 1] = asumi;
            self.bsumr[j - 1] = bsumr;
            self.bsumi[j - 1] = bsumi;

            let (s1r, s1i);
            if self.kode != 1 {
                let mut str = self.zbr + zeta2r;
                let mut sti = self.zbi + zeta2i;
                let rast = fn_val / zabs(str, sti);
                str = str * rast * rast;
                sti = -sti * rast * rast;
                s1r = zeta1r - str;
                s1i = zeta1i - sti;
            } else {
                s1r = zeta1r - zeta2r;
                s1i = zeta1i - zeta2i;
            }

            let mut rs1 = s1r;

            if rs1.abs() > self.elim {
                match self.underflow_70(i, rs1) {
                    Flow::Return => return Flow::Return,
                    Flow::Continue => continue,
                }
            }

            if self.kdflg == 1 {
                self.kflag = 2;
            }

            if rs1.abs() >= self.alim {
                let aphi = zabs(phir, phii);
                let aarg = zabs(argr, argi);
                rs1 = rs1 + aphi.ln() - 0.25 * aarg.ln() - AIC;
                if rs1.abs() > self.elim {
                    match self.underflow_70(i, rs1) {
                        Flow::Return => return Flow::Return,
                        Flow::Continue => continue,
                    }
                }
                if self.kdflg == 1 {
                    self.kflag = 1;
                }
                if rs1 >= 0.0 && self.kdflg == 1 {
                    self.kflag = 3;
                }
            }

            // Scale S1.
            let mut c2r = argr * CR2R - argi * CR2I;
            let mut c2i = argr * CR2I + argi * CR2R;
            let mut air = 0.0f64;
            let mut aii = 0.0f64;
            let mut dair = 0.0f64;
            let mut daii = 0.0f64;
            let mut nai: i32 = 0;
            let mut ndai: i32 = 0;
            let mut idum: i32 = 0;
            zairy(c2r, c2i, 0, 2, &mut air, &mut aii, &mut nai, &mut idum);
            zairy(c2r, c2i, 1, 2, &mut dair, &mut daii, &mut ndai, &mut idum);

            let mut str = dair * bsumr - daii * bsumi;
            let mut sti = dair * bsumi + daii * bsumr;
            let ptr = str * CR2R - sti * CR2I;
            let pti = str * CR2I + sti * CR2R;
            str = ptr + (air * asumr - aii * asumi);
            sti = pti + (air * asumi + aii * asumr);
            let ptr = str * phir - sti * phii;
            let pti = str * phii + sti * phir;
            let mut s2r = ptr * self.csr - pti * self.csi;
            let mut s2i = ptr * self.csi + pti * self.csr;

            let str = rs1.exp() * self.cssr[(self.kflag - 1) as usize];
            let s1r_new = str * s1i.cos();
            let s1i_new = str * s1i.sin();
            let str2 = s2r * s1r_new - s2i * s1i_new;
            s2i = s1r_new * s2i + s2r * s1i_new;
            s2r = str2;

            if self.kflag == 1 {
                let mut nw: i32 = 0;
                zuchk(s2r, s2i, &mut nw, self.bry[0], self.tol);
                if nw != 0 {
                    match self.underflow_70(i, rs1) {
                        Flow::Return => return Flow::Return,
                        Flow::Continue => continue,
                    }
                }
            }

            if self.yy <= 0.0 {
                s2i = -s2i;
            }
            self.cyr[(self.kdflg - 1) as usize] = s2r;
            self.cyi[(self.kdflg - 1) as usize] = s2i;
            self.yr[(i - 1) as usize] = s2r * self.csrr[(self.kflag - 1) as usize];
            self.yi[(i - 1) as usize] = s2i * self.csrr[(self.kflag - 1) as usize];

            // Rotate CSR/CSI.
            let str = self.csi;
            self.csi = -self.csr;
            self.csr = str;

            if self.kdflg == 2 {
                self.fn_val = fn_val;
                self.last_i = i;
                return Flow::Continue;
            }
            self.kdflg = 2;
        }

        self.fn_val = self.fnu + ((self.n - 1) as f64);
        self.last_i = self.n;
        Flow::Continue
    }

    // -----------------------------------------------------------------
    // Fortran label-70 underflow/overflow.
    // -----------------------------------------------------------------
    fn underflow_70(&mut self, i: i32, rs1: f64) -> Flow {
        if rs1 > 0.0 {
            self.nz = -1;
            return Flow::Return;
        }
        if self.zr < 0.0 {
            self.nz = -1;
            return Flow::Return;
        }
        self.kdflg = 1;
        self.yr[(i - 1) as usize] = ZEROR;
        self.yi[(i - 1) as usize] = ZEROI;
        self.nz += 1;
        let str = self.csi;
        self.csi = -self.csr;
        self.csr = str;
        if i == 1 {
            return Flow::Continue;
        }
        if self.yr[(i - 2) as usize] == ZEROR && self.yi[(i - 2) as usize] == ZEROI {
            return Flow::Continue;
        }
        self.yr[(i - 2) as usize] = ZEROR;
        self.yi[(i - 2) as usize] = ZEROI;
        self.nz += 1;
        Flow::Continue
    }

    // -----------------------------------------------------------------
    // Test last member for underflow and overflow.
    // -----------------------------------------------------------------
    fn last_member_check(&mut self) -> Flow {
        let fn_val = self.fnu + ((self.n - 1) as f64);
        let ipard = if self.mr != 0 { 0 } else { 1 };

        let (mut phidr, mut phidi) = (0.0f64, 0.0f64);
        let (mut argdr, mut argdi) = (0.0f64, 0.0f64);
        let (mut zet1dr, mut zet1di) = (0.0f64, 0.0f64);
        let (mut zet2dr, mut zet2di) = (0.0f64, 0.0f64);
        let (mut asumdr, mut asumdi) = (0.0f64, 0.0f64);
        let (mut bsumdr, mut bsumdi) = (0.0f64, 0.0f64);

        zunhj(
            self.znr, self.zni, fn_val, ipard, self.tol, &mut phidr, &mut phidi,
            &mut argdr, &mut argdi, &mut zet1dr, &mut zet1di, &mut zet2dr, &mut zet2di,
            &mut asumdr, &mut asumdi, &mut bsumdr, &mut bsumdi,
        );

        let (s1r, s1i);
        if self.kode != 1 {
            let mut str = self.zbr + zet2dr;
            let mut sti = self.zbi + zet2di;
            let rast = fn_val / zabs(str, sti);
            str = str * rast * rast;
            sti = -sti * rast * rast;
            s1r = zet1dr - str;
            s1i = zet1di - sti;
        } else {
            s1r = zet1dr - zet2dr;
            s1i = zet1di - zet2di;
        }

        let mut rs1 = s1r;

        if rs1.abs() > self.elim {
            if rs1 > 0.0 {
                self.nz = -1;
                return Flow::Return;
            }
            if self.zr < 0.0 {
                self.nz = -1;
                return Flow::Return;
            }
            self.nz = self.n;
            for i in 0..(self.n as usize) {
                self.yr[i] = ZEROR;
                self.yi[i] = ZEROI;
            }
            return Flow::Return;
        }

        if rs1.abs() >= self.alim {
            let aphi = zabs(phidr, phidi);
            rs1 += aphi.ln();
            if rs1.abs() >= self.elim {
                if rs1 > 0.0 {
                    self.nz = -1;
                    return Flow::Return;
                }
                if self.zr < 0.0 {
                    self.nz = -1;
                    return Flow::Return;
                }
                self.nz = self.n;
                for i in 0..(self.n as usize) {
                    self.yr[i] = ZEROR;
                    self.yi[i] = ZEROI;
                }
                return Flow::Return;
            }
        }

        Flow::Continue
    }

    // -----------------------------------------------------------------
    // Forward recurrence (Fortran DO 130 I=IB,N).
    // -----------------------------------------------------------------
    fn forward_recur(&mut self) {
        let mut s1r = self.cyr[0];
        let mut s1i = self.cyi[0];
        let mut s2r = self.cyr[1];
        let mut s2i = self.cyi[1];
        let mut c1r = self.csrr[(self.kflag - 1) as usize];
        let mut ascle = self.bry[(self.kflag - 1) as usize];

        let mut ckr = self.ckr;
        let mut cki = self.cki;

        let mut i = self.ib;
        while i <= self.n {
            let c2r = s2r;
            let c2i = s2i;
            s2r = ckr * c2r - cki * c2i + s1r;
            s2i = ckr * c2i + cki * c2r + s1i;
            s1r = c2r;
            s1i = c2i;
            ckr += self.rzr;
            cki += self.rzi;
            let c2r = s2r * c1r;
            let c2i = s2i * c1r;
            self.yr[(i - 1) as usize] = c2r;
            self.yi[(i - 1) as usize] = c2i;

            if self.kflag < 3 {
                let c2m = c2r.abs().max(c2i.abs());
                if c2m > ascle {
                    self.kflag += 1;
                    ascle = self.bry[(self.kflag - 1) as usize];
                    s1r *= c1r;
                    s1i *= c1r;
                    s2r = c2r;
                    s2i = c2i;
                    s1r *= self.cssr[(self.kflag - 1) as usize];
                    s1i *= self.cssr[(self.kflag - 1) as usize];
                    s2r *= self.cssr[(self.kflag - 1) as usize];
                    s2i *= self.cssr[(self.kflag - 1) as usize];
                    c1r = self.csrr[(self.kflag - 1) as usize];
                }
            }
            i += 1;
        }
    }

    // -----------------------------------------------------------------
    // Continuation loop (Fortran DO 290 K=1,N).
    // -----------------------------------------------------------------
    fn continuation_loop(&mut self) -> Flow {
        for _k in 1..=self.n {
            let fn_val = self.fnu + ((self.kk - 1) as f64);

            self.sort_params();

            let mut phidr = self.pending_phidr;
            let mut phidi = self.pending_phidi;
            let mut argdr = self.pending_argdr;
            let mut argdi = self.pending_argdi;
            let mut zet1dr = self.pending_zet1dr;
            let mut zet1di = self.pending_zet1di;
            let mut zet2dr = self.pending_zet2dr;
            let mut zet2di = self.pending_zet2di;
            let mut asumdr = self.pending_asumdr;
            let mut asumdi = self.pending_asumdi;
            let mut bsumdr = self.pending_bsumdr;
            let mut bsumdi = self.pending_bsumdi;
            let called_zunhj = self.pending_called_zunhj;

            if called_zunhj {
                zunhj(
                    self.znr, self.zni, fn_val, 0, self.tol, &mut phidr, &mut phidi,
                    &mut argdr, &mut argdi, &mut zet1dr, &mut zet1di, &mut zet2dr,
                    &mut zet2di, &mut asumdr, &mut asumdi, &mut bsumdr, &mut bsumdi,
                );
            }

            let (s1r_in, s1i_in);
            if self.kode != 1 {
                let mut str = self.zbr + zet2dr;
                let mut sti = self.zbi + zet2di;
                let rast = fn_val / zabs(str, sti);
                str = str * rast * rast;
                sti = -sti * rast * rast;
                s1r_in = -zet1dr + str;
                s1i_in = -zet1di + sti;
            } else {
                s1r_in = -zet1dr + zet2dr;
                s1i_in = -zet1di + zet2di;
            }

            let mut rs1 = s1r_in;
            let mut s1r = s1r_in;
            let mut s1i = s1i_in;

            if rs1.abs() > self.elim {
                if rs1 > 0.0 {
                    self.nz = -1;
                    return Flow::Return;
                }
                let s2r = ZEROR;
                let s2i = ZEROI;
                let (c2r_out, c2i_out) = self.accumulate(s2r, s2i);
                self.kk -= 1;
                let str = self.csi;
                self.csi = -self.csr;
                self.csr = str;
                if c2r_out != 0.0 || c2i_out != 0.0 {
                    if self.kdflg == 2 {
                        return Flow::Continue;
                    }
                    self.kdflg = 2;
                } else {
                    self.kdflg = 1;
                }
                continue;
            }

            if self.kdflg == 1 {
                self.iflag = 2;
            }

            if rs1.abs() >= self.alim {
                let aphi = zabs(phidr, phidi);
                let aarg = zabs(argdr, argdi);
                rs1 = rs1 + aphi.ln() - 0.25 * aarg.ln() - AIC;
                if rs1.abs() > self.elim {
                    if rs1 > 0.0 {
                        self.nz = -1;
                        return Flow::Return;
                    }
                    let s2r = ZEROR;
                    let s2i = ZEROI;
                    let (c2r_out, c2i_out) = self.accumulate(s2r, s2i);
                    self.kk -= 1;
                    let str = self.csi;
                    self.csi = -self.csr;
                    self.csr = str;
                    if c2r_out != 0.0 || c2i_out != 0.0 {
                        if self.kdflg == 2 {
                            return Flow::Continue;
                        }
                        self.kdflg = 2;
                    } else {
                        self.kdflg = 1;
                    }
                    continue;
                }
                if self.kdflg == 1 {
                    self.iflag = 1;
                }
                if rs1 >= 0.0 && self.kdflg == 1 {
                    self.iflag = 3;
                }
            }

            let mut air = 0.0f64;
            let mut aii = 0.0f64;
            let mut dair = 0.0f64;
            let mut daii = 0.0f64;
            let mut nai: i32 = 0;
            let mut ndai: i32 = 0;
            let mut idum: i32 = 0;
            zairy(argdr, argdi, 0, 2, &mut air, &mut aii, &mut nai, &mut idum);
            zairy(argdr, argdi, 1, 2, &mut dair, &mut daii, &mut ndai, &mut idum);

            let mut str = dair * bsumdr - daii * bsumdi;
            let mut sti = dair * bsumdi + daii * bsumdr;
            str += air * asumdr - aii * asumdi;
            sti += air * asumdi + aii * asumdr;
            let ptr = str * phidr - sti * phidi;
            let pti = str * phidi + sti * phidr;
            let mut s2r = ptr * self.csr - pti * self.csi;
            let mut s2i = ptr * self.csi + pti * self.csr;

            let str = rs1.exp() * self.cssr[(self.iflag - 1) as usize];
            s1r = str * s1i.cos();
            s1i = str * s1i.sin();
            let str2 = s2r * s1r - s2i * s1i;
            s2i = s2r * s1i + s2i * s1r;
            s2r = str2;

            if self.iflag == 1 {
                let mut nw: i32 = 0;
                zuchk(s2r, s2i, &mut nw, self.bry[0], self.tol);
                if nw != 0 {
                    s2r = ZEROR;
                    s2i = ZEROI;
                }
            }

            if self.yy <= 0.0 {
                s2i = -s2i;
            }

            let (c2r_out, c2i_out) = self.accumulate(s2r, s2i);
            self.kk -= 1;
            let str = self.csi;
            self.csi = -self.csr;
            self.csr = str;

            if c2r_out != 0.0 || c2i_out != 0.0 {
                if self.kdflg == 2 {
                    return Flow::Continue;
                }
                self.kdflg = 2;
            } else {
                self.kdflg = 1;
            }
        }
        Flow::Continue
    }

    // -----------------------------------------------------------------
    // Sort logic for continuation parameters.
    // -----------------------------------------------------------------
    fn sort_params(&mut self) {
        if self.n > 2 {
            if (self.kk == self.n) && (self.ib < self.n) {
                self.pending_called_zunhj = true;
                return;
            }
            if self.kk == self.ib || self.kk == self.ic {
                self.do_sort_172();
                self.pending_called_zunhj = false;
                return;
            }
            self.pending_called_zunhj = true;
            return;
        }
        self.do_sort_172();
        self.pending_called_zunhj = false;
    }

    fn do_sort_172(&mut self) {
        let j = self.j as usize;
        self.pending_phidr = self.phir[j - 1];
        self.pending_phidi = self.phii[j - 1];
        self.pending_argdr = self.argr[j - 1];
        self.pending_argdi = self.argi[j - 1];
        self.pending_zet1dr = self.zeta1r[j - 1];
        self.pending_zet1di = self.zeta1i[j - 1];
        self.pending_zet2dr = self.zeta2r[j - 1];
        self.pending_zet2di = self.zeta2i[j - 1];
        self.pending_asumdr = self.asumr[j - 1];
        self.pending_asumdi = self.asumi[j - 1];
        self.pending_bsumdr = self.bsumr[j - 1];
        self.pending_bsumdi = self.bsumi[j - 1];
        self.j = 3 - self.j;
    }

    // -----------------------------------------------------------------
    // Accumulate (Fortran label 250).
    // -----------------------------------------------------------------
    fn accumulate(&mut self, s2r_in: f64, s2i_in: f64) -> (f64, f64) {
        let c2r = s2r_in;
        let c2i = s2i_in;

        let kdflg = self.kdflg as usize;
        self.cyr[kdflg - 1] = s2r_in;
        self.cyi[kdflg - 1] = s2i_in;

        let mut s2r = s2r_in * self.csrr[(self.iflag - 1) as usize];
        let mut s2i = s2i_in * self.csrr[(self.iflag - 1) as usize];

        let kk = self.kk as usize;
        let mut s1r = self.yr[kk - 1];
        let mut s1i = self.yi[kk - 1];

        if self.kode != 1 {
            let mut nw: i32 = 0;
            let mut iuf = self.iuf;
            zs1s2(
                self.zrr, self.zri, &mut s1r, &mut s1i, &mut s2r, &mut s2i, &mut nw,
                self.asc, self.alim, &mut iuf,
            );
            self.iuf = iuf;
            self.nz += nw;
        }

        self.yr[kk - 1] = s1r * self.cspnr - s1i * self.cspni + s2r;
        self.yi[kk - 1] = s1r * self.cspni + s1i * self.cspnr + s2i;

        (c2r, c2i)
    }

    // -----------------------------------------------------------------
    // Backward recurrence for remainder of I sequence.
    // -----------------------------------------------------------------
    fn backward_recur(&mut self, il: i32) {
        let mut s1r = self.cyr[0];
        let mut s1i = self.cyi[0];
        let mut s2r = self.cyr[1];
        let mut s2i = self.cyi[1];
        let mut csr = self.csrr[(self.iflag - 1) as usize];
        let mut ascle = self.bry[(self.iflag - 1) as usize];
        let mut fn_val = (self.inu + il) as f64;

        for _i in 1..=il {
            let c2r = s2r;
            let c2i = s2i;
            s2r = s1r + (fn_val + self.fnf) * (self.rzr * c2r - self.rzi * c2i);
            s2i = s1i + (fn_val + self.fnf) * (self.rzr * c2i + self.rzi * c2r);
            s1r = c2r;
            s1i = c2i;
            fn_val -= 1.0;
            let c2r = s2r * csr;
            let c2i = s2i * csr;
            let ckr = c2r;
            let cki = c2i;
            let kk = self.kk as usize;
            let mut c1r = self.yr[kk - 1];
            let mut c1i = self.yi[kk - 1];

            if self.kode != 1 {
                let mut nw: i32 = 0;
                let mut iuf = self.iuf;
                let mut c2r_mut = c2r;
                let mut c2i_mut = c2i;
                zs1s2(
                    self.zrr, self.zri, &mut c1r, &mut c1i, &mut c2r_mut, &mut c2i_mut,
                    &mut nw, self.asc, self.alim, &mut iuf,
                );
                self.iuf = iuf;
                self.nz += nw;
            }

            self.yr[kk - 1] = c1r * self.cspnr - c1i * self.cspni + c2r;
            self.yi[kk - 1] = c1r * self.cspni + c1i * self.cspnr + c2i;
            self.kk -= 1;
            self.cspnr = -self.cspnr;
            self.cspni = -self.cspni;

            if self.iflag < 3 {
                let c2m = ckr.abs().max(cki.abs());
                if c2m > ascle {
                    self.iflag += 1;
                    ascle = self.bry[(self.iflag - 1) as usize];
                    s1r *= csr;
                    s1i *= csr;
                    s2r = ckr;
                    s2i = cki;
                    s1r *= self.cssr[(self.iflag - 1) as usize];
                    s1i *= self.cssr[(self.iflag - 1) as usize];
                    s2r *= self.cssr[(self.iflag - 1) as usize];
                    s2i *= self.cssr[(self.iflag - 1) as usize];
                    csr = self.csrr[(self.iflag - 1) as usize];
                }
            }
        }
    }
}