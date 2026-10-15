pub(crate) mod i1mach;
use i1mach::i1mach;  // fn i1mach(i: i32) -> i32
pub(crate) mod d1mach;
use d1mach::d1mach;  // fn d1mach(i: i32) -> f64

pub(crate) mod zabs;
use zabs::zabs;    // fn zabs(zr: f64, zi: f64) -> f64
pub(crate) mod zmlt;
use zmlt::zmlt;    // fn zmlt(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64)
pub(crate) mod zdiv;
use zdiv::zdiv;    // fn zdiv(ar: f64, ai: f64, br: f64, bi: f64, cr: &mut f64, ci: &mut f64)
pub(crate) mod zlog;
use zlog::zlog;    // fn zlog(ar: f64, ai: f64, br: &mut f64, bi: &mut f64) -> i32
pub(crate) mod zexp;
use zexp::zexp;    // fn zexp(ar: f64, ai: f64, br: &mut f64, bi: &mut f64)
