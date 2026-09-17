#![allow(non_snake_case, dead_code)]
use crate::consts::LN2;

pub(crate) fn R_forceint(x: f64) -> f64 {
    f64::round_ties_even(x)
}
pub(crate) fn fmax2(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        f64::NAN
    } else {
        f64::max(x, y)
    }
}
pub(crate) fn fmin2(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        f64::NAN
    } else {
        f64::min(x, y)
    }
}
pub(crate) fn R_nonint(x: f64) -> bool {
    f64::abs(x) - R_forceint(x) > 1.0e-7 * fmax2(1.0, f64::abs(x))
}

pub(crate) fn ISNAN(x: f64) -> bool {
    x.is_nan()
}
pub(crate) fn R_FINITE(x: f64) -> bool {
    x.is_finite()
}
pub(crate) fn R_D__0(log_p: bool) -> f64 {
    if log_p { f64::NEG_INFINITY } else { 0.0 }
}
pub(crate) fn R_D__1(log_p: bool) -> f64 {
    if log_p { 0.0 } else { 1.0 }
}
pub(crate) fn R_DT_0(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D__0(log_p)
    } else {
        R_D__1(log_p)
    }
}
pub(crate) fn R_DT_1(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D__1(log_p)
    } else {
        R_D__0(log_p)
    }
}
pub(crate) fn R_D_half(log_p: bool) -> f64 {
    if log_p { -LN2 } else { 0.5 }
}

pub(crate) fn R_D_Lval(p: f64, lower_tail: bool) -> f64 {
    if lower_tail { p } else { 0.5 - p + 0.5 }
}
pub(crate) fn R_D_Cval(p: f64, lower_tail: bool) -> f64 {
    if lower_tail { 0.5 - p + 0.5 } else { p }
}

pub(crate) fn R_D_val(x: f64, log_p: bool) -> f64 {
    if log_p { f64::ln(x) } else { x }
}
pub(crate) fn R_D_qIv(p: f64, log_p: bool) -> f64 {
    if log_p { f64::exp(p) } else { p }
}
pub(crate) fn R_D_exp(x: f64, log_p: bool) -> f64 {
    if log_p { x } else { f64::exp(x) }
}
pub(crate) fn R_D_log(p: f64, log_p: bool) -> f64 {
    if log_p { p } else { f64::ln(p) }
}
pub(crate) fn R_D_Clog(p: f64, log_p: bool) -> f64 {
    if log_p { f64::ln_1p(-p) } else { 0.5 - p + 0.5 }
}

// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
pub(crate) fn R_Log1_Exp(x: f64) -> f64 {
    if x > -LN2 {
        f64::ln(-f64::exp_m1(x))
    } else {
        f64::ln_1p(-f64::exp(x))
    }
}

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
pub(crate) fn R_D_LExp(x: f64, log_p: bool) -> f64 {
    if log_p { R_Log1_Exp(x) } else { f64::ln_1p(-x) }
}

pub(crate) fn R_DT_val(x: f64, lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D_val(x, log_p)
    } else {
        R_D_Clog(x, log_p)
    }
}
pub(crate) fn R_DT_Cval(x: f64, lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D_Clog(x, log_p)
    } else {
        R_D_val(x, log_p)
    }
}

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
pub(crate) fn R_DT_qIv(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if log_p {
        if lower_tail {
            f64::exp(p)
        } else {
            -f64::exp_m1(p)
        }
    } else {
        R_D_Lval(p, lower_tail)
    }
}

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
pub(crate) fn R_DT_CIv(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if log_p {
        if lower_tail {
            -f64::exp_m1(p)
        } else {
            f64::exp(p)
        }
    } else {
        R_D_Cval(p, lower_tail)
    }
}

/* exp(x) */
pub(crate) fn R_DT_exp(x: f64, lower_tail: bool, log_p: bool) -> f64 {
    R_D_exp(R_D_Lval(x, lower_tail), log_p)
}
/* exp(1 - x) */
pub(crate) fn R_DT_Cexp(x: f64, lower_tail: bool, log_p: bool) -> f64 {
    R_D_exp(R_D_Cval(x, lower_tail), log_p)
}
/* log(p) in qF */
pub(crate) fn R_DT_log(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D_log(p, log_p)
    } else {
        R_D_LExp(p, log_p)
    }
}
/* log(1-p) in qF*/
pub(crate) fn R_DT_Clog(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        R_D_LExp(p, log_p)
    } else {
        R_D_log(p, log_p)
    }
}
pub(crate) fn R_DT_Log(p: f64, lower_tail: bool) -> f64 {
    if lower_tail { p } else { R_Log1_Exp(p) }
}
// ==   R_DT_log when we already "know" log_p == TRUE

// ML_WARN_return_NAN
pub(crate) fn R_Q_P01_check(p: f64, log_p: bool) -> Option<f64> {
    if (log_p && p > 0.0) || (!log_p && (p < 0.0 || p > 1.0)) {
        return Some(f64::NAN);
    }
    None
}

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
pub(crate) fn R_Q_P01_boundaries(
    p: f64,
    _LEFT_: f64,
    _RIGHT_: f64,
    lower_tail: bool,
    log_p: bool,
) -> Option<f64> {
    if log_p {
        if p > 0.0 {
            return Some(f64::NAN);
        }
        if p == 0.0
        /* upper bound*/
        {
            return Some(if lower_tail { _RIGHT_ } else { _LEFT_ });
        }
        if p == f64::NEG_INFINITY {
            return Some(if lower_tail { _LEFT_ } else { _RIGHT_ });
        }
    } else {
        /* !log_p */
        if p < 0.0 || p > 1.0 {
            return Some(f64::NAN);
        }
        if p == 0.0 {
            return Some(if lower_tail { _LEFT_ } else { _RIGHT_ });
        }
        if p == 1.0 {
            return Some(if lower_tail { _RIGHT_ } else { _LEFT_ });
        }
    }
    None
}

pub(crate) fn R_P_bounds_01(
    x: f64,
    x_min: f64,
    x_max: f64,
    lower_tail: bool,
    log_p: bool,
) -> Option<f64> {
    if x <= x_min {
        return Some(R_DT_0(lower_tail, log_p));
    }
    if x >= x_max {
        return Some(R_DT_1(lower_tail, log_p));
    }
    None
}
/* is typically not quite optimal for (-Inf,Inf) where
 * you'd rather have */
pub(crate) fn R_P_bounds_Inf_01(x: f64, lower_tail: bool, log_p: bool) -> Option<f64> {
    if !R_FINITE(x) {
        if x > 0.0 {
            return Some(R_DT_1(lower_tail, log_p));
        }
        /* x < 0 */
        return Some(R_DT_0(lower_tail, log_p));
    }
    None
}

/* additions for density functions (C.Loader) */
pub(crate) fn R_D_fexp(f: f64, x: f64, give_log: bool) -> f64 {
    if give_log {
        -0.5 * f64::ln(f) + x
    } else {
        f64::exp(x) / f64::sqrt(f)
    }
}
// version working with rf := sqrt(f) [avoiding overflow in computation of f in the caller]
pub(crate) fn R_D_rtxp(rf: f64, x: f64, give_log: bool) -> f64 {
    if give_log {
        -f64::ln(rf) + x
    } else {
        f64::exp(x) / rf
    }
}

/* [neg]ative or [non int]eger : */
pub(crate) fn R_D_negInonint(x: f64) -> bool {
    x < 0. || R_nonint(x)
}

// for discrete d<distr>(x, ...) :
pub(crate) fn R_D_nonint_check(x: f64, log_p: bool) -> Option<f64> {
    if R_nonint(x) {
        return Some(R_D__0(log_p));
    }
    return None;
}
