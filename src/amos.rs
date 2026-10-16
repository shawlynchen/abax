pub(crate) mod i1mach;
use i1mach::i1mach;  // fn i1mach(i: i32) -> i32
pub(crate) mod d1mach;
use d1mach::d1mach;  // fn d1mach(i: i32) -> f64

pub(crate) mod dgamln;
use dgamln::dgamln;  // fn dgamln(z: f64, ierr: &mut i32) -> f64

pub(crate) mod zasyi;
use zasyi::zasyi;  // fn zasyi(zr: f64, zi: f64, fnu: f64, kode: i32, n: i32, yr: &mut [f64], yi: &mut [f64], nz: &mut i32, rl: f64, tol: f64, elim: f64, alim: f64)

pub(crate) mod zabs;
use zabs::zabs;    // fn zabs(zr: f64, zi: f64) -> f64
pub(crate) mod zmlt;
use zmlt::zmlt;    // fn zmlt(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64)
pub(crate) mod zdiv;
use zdiv::zdiv;    // fn zdiv(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64)
pub(crate) mod zsqrt;
use zsqrt::zsqrt;  // fn zsqrt(ar: f64, ai: f64, br: &mut f64, bi: &mut f64)
pub(crate) mod zlog;
use zlog::zlog;    // fn zlog(ar: f64, ai: f64, br: &mut f64, bi: &mut f64) -> i32
pub(crate) mod zexp;
use zexp::zexp;    // fn zexp(ar: f64, ai: f64, br: &mut f64, bi: &mut f64)
pub(crate) mod zshch;
use zshch::zshch;  // fn zshch(zr: f64, zi: f64, cshr: &mut f64, cshi: &mut f64, cchr: &mut f64, cchi: &mut f64)
pub(crate) mod zuchk;
use zuchk::zuchk;  // fn zuchk(yr: f64, yi: f64, nz: &mut i32, ascle: f64, tol: f64)
