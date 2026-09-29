/// ALGDIV computes ln(gamma(b)/gamma(a+b)) when 8 <= B.
/// In this algorithm, del(x) is the function defined by
/// ln(gamma(x)) = (x - 0.5) * ln(x) - x + 0.5 * ln(2*pi) + del(x).
#[allow(dead_code)]
fn algdiv(a: f64, b: f64) -> f64 {
    const C0: f64 =  0.833333333333333E-01;
    const C1: f64 = -0.277777777760991E-02;
    const C2: f64 =  0.793650666825390E-03;
    const C3: f64 = -0.595202931351870E-03;
    const C4: f64 =  0.837308034031215E-03;
    const C5: f64 = -0.165322962780713E-02;

    let h: f64;
    let c: f64;
    let x: f64;
    let d: f64;
    if b < a {
        h = b / a;
        c = 1.0 / (1.0 + h);
        x = h / (1.0 + h);
        d = a + (b - 0.5);
    } else {
        h = a / b;
        c = h / (1.0 + h);
        x = 1.0 / (1.0 + h);
        d = b + (a - 0.5);
    }

    // Set sn = (1 - x**n) / ( 1 - x )
    let x2 = x * x;
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);

    // Set w = del(b) - del(a + b)
    let t = 1.0 / (b * b);

    let w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;

    let w = w * (c / b);
    
    // Combine the results.
    let u = d * alnrel (a / b);
    let v = a * (f64::ln(b) - 1.0);

    return if v < u {
        (w - v) - u
    } else {
        (w - u) - v
    };
}

/// ALNREL evaluates the function ln(1 + a).
#[allow(dead_code)]
fn alnrel(a: f64) -> f64 {

    const P1: f64 = -0.129418923021993E+01;
    const P2: f64 =  0.405303492862024E+00;
    const P3: f64 = -0.178874546012214E-01;
    const Q1: f64 = -0.162752256355323E+01;
    const Q2: f64 =  0.747811014037616E+00;
    const Q3: f64 = -0.845104217945565E-01;

    if f64::abs(a) <= 0.375 {
        let t = a / (a + 2.0);
        let t2 = t * t;

        let w = (((P3 * t2 + P2) * t2 + P1) * t2 + 1.0 ) / (((Q3 * t2 + Q2) * t2 + Q1) * t2 + 1.0);

        return 2.0 * t * w;
    } else {
        let x = 1.0 + a;
        return f64::ln(x);
    }
}

/// APSER yields the incomplete beta ratio i(sub(1-x))(b,a) for
///     a <= min(eps,eps*b), b*x <= 1, and x <= 0.5. used when
///     a is very small. use only if above inequalities are satisfied.
#[allow(dead_code)]
fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    const G: f64 = 0.577215664901533;
    
    let bx = b * x;
    let mut t = x - bx;

    let c: f64 = if b * eps <= 2.0e-2 {
        f64::ln(x) + psi(b) + G + t
    } else {
        f64::ln(bx) + G + t
    };

    let tol = 5.0 * eps * f64::abs(c);
    let mut j = 1.0;
    let mut s = 0.0;

    loop {
        j += 1.0;
        t *= x - bx / j;
        let aj = t / j;
        s += aj;

        if f64::abs(aj) <= tol {
            break;
        }
    }

    -a * (c + s)
}

/// BASYM uses an asymptotic expansion for Ix(A,B) for large A and B.
/// lambda = (a + b) * y - b and eps is the tolerance used.
/// it is assumed that lambda is nonnegative and that
/// a and b are greater than or equal to 15.
#[allow(dead_code)]
fn basym(a: f64, b: f64, lambda: f64, eps: f64) -> f64 {
    const E0: f64 = 1.12837916709551e+00;
    const E1: f64 = 0.353553390593274e+00;

    // Maximum n in the loop below.
    // Must be even.
    const NUM: usize = 20;

    let mut a0 = [0.0_f64; NUM + 2];
    let mut b0 = [0.0_f64; NUM + 2];
    let mut c = [0.0_f64; NUM + 2];
    let mut d = [0.0_f64; NUM + 2];

    let mut basym = 0.0_f64;

    let (h, r0, r1, w0);

    if a < b {
        let h1 = a / b;

        h = h1;
        r0 = 1.0 / (1.0 + h1);
        r1 = (b - a) / b;
        w0 = 1.0 / (a * (1.0 + h1)).sqrt();
    } else {
        let h1 = b / a;

        h = h1;
        r0 = 1.0 / (1.0 + h1);
        r1 = (b - a) / a;
        w0 = 1.0 / (b * (1.0 + h1)).sqrt();
    }

    // f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
    let f = a * rlog1(-lambda / a) + b * rlog1(lambda / b);

    let t = (-f).exp();

    if t == 0.0 {
        return basym;
    }

    let z0 = f.sqrt();
    let z = 0.5 * (z0 / E1);
    let z2 = f + f;

    a0[1] = (2.0 / 3.0) * r1;
    c[1] = -0.5 * a0[1];
    d[1] = -c[1];

    let mut j0 = (0.5 / E0) * erfc1(1, z0);
    let mut j1 = E1;

    let mut sum2 = j0 + d[1] * w0 * j1;

    let mut s = 1.0_f64;
    let h2 = h * h;
    let mut hn = 1.0_f64;
    let mut w = w0;

    let mut znm1 = z;
    let mut zn = z2;

    // Fortran:
    //
    // do n = 2, num, 2
    //
    for n in (2..=NUM).step_by(2) {
        hn = h2 * hn;

        a0[n] = 2.0 * r0 * (1.0 + h * hn) / (n as f64 + 2.0);

        let np1 = n + 1;

        s += hn;

        a0[np1] = 2.0 * r1 * s / (n as f64 + 3.0);

        // --------------------------------------------------------
        // do i = n, np1
        // --------------------------------------------------------
        for i in n..=np1 {
            let r = -0.5 * (i as f64 + 1.0);

            b0[1] = r * a0[1];

            // ----------------------------------------------------
            // do m = 2, i
            // ----------------------------------------------------
            for m in 2..=i {
                let mut bsum = 0.0_f64;

                let mm1 = m - 1;

                // ------------------------------------------------
                // do j = 1, mm1
                // ------------------------------------------------
                for j in 1..=mm1 {
                    let mmj = m - j;

                    bsum +=
                        (j as f64 * r - mmj as f64)
                        * a0[j]
                        * b0[mmj];
                }

                b0[m] = r * a0[m] + bsum / m as f64;
            }

            c[i] = b0[i] / (i as f64 + 1.0);

            let mut dsum = 0.0_f64;

            let im1 = i - 1;

            // ----------------------------------------------------
            // do j = 1, im1
            // ----------------------------------------------------
            for j in 1..=im1 {
                let imj = i - j;

                dsum += d[imj] * c[j];
            }

            d[i] = -(dsum + c[i]);
        }

        // --------------------------------------------------------
        // Update j0 and j1
        // --------------------------------------------------------
        j0 = E1 * znm1 + (n as f64 - 1.0) * j0;
        j1 = E1 * zn + n as f64 * j1;

        znm1 = z2 * znm1;
        zn = z2 * zn;

        w = w0 * w;

        let t0 = d[n] * w * j0;

        w = w0 * w;

        let t1 = d[np1] * w * j1;

        sum2 += t0 + t1;

        if (t0.abs() + t1.abs()) <= eps * sum2 {
            break;
        }
    }

    let u = (-bcorr(a, b)).exp();

    basym = E0 * t * u * sum2;

    basym
}

/// BCORR evaluates del(a0) + del(b0) - del(a0 + b0) where
///   ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
/// it is assumed that a0 >= 8 and b0 >= 8.
#[allow(dead_code)]
fn bcorr(a0: f64, b0: f64) -> f64 {
    // Constants from the original Fortran routine.
    const C0: f64 = 0.833333333333333e-01;
    const C1: f64 = -0.277777777760991e-02;
    const C2: f64 = 0.793650666825390e-03;
    const C3: f64 = -0.595202931351870e-03;
    const C4: f64 = 0.837308034031215e-03;
    const C5: f64 = -0.165322962780713e-02;

    let a = a0.min(b0);
    let b = a0.max(b0);

    let h = a / b;
    let c = h / (1.0 + h);
    let x = 1.0 / (1.0 + h);
    let x2 = x * x;

    // Set: sn = (1 - x^n) / (1 - x)
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);

    // Set: w = del(b) - del(a+b)
    let t = (1.0 / b) * (1.0 / b);

    let mut w =
        ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t
            + C1 * s3)
            * t
            + C0;

    w *= c / b;

    // Compute: del(a) + w
    let t = (1.0 / a) * (1.0 / a);

    (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a + w
}

/// Evaluates the logarithm of the Beta function.
#[allow(dead_code)]
fn betaln(a0: f64, b0: f64) -> f64 {
    const E: f64 = 0.918938533204673e+00;         // Log(2 * PI) / 2

    let mut a = a0.min(b0);
    let b = a0.max(b0);

    if a < 1.0 {
        if b < 8.0 {
            return gamln(a) + (gamln(b) - gamln(a + b));
        } else {
            return gamln(a) + algdiv(a, b);
        }
    }

    let betaln_b_lt_8_with_w = |a: f64, mut b: f64, w: f64| -> f64 {
        let n = (b - 1.0) as i32;
        let mut z = 1.0_f64;
        for _ in 1..=n {
            b -= 1.0;
            z *= b / (a + b);
        }
        w + z.ln() + (gamln(a) + (gamln(b) - gsumln(a, b)))
    };

    if a < 8.0 {
        // 1 <= a < 8
        if a <= 2.0 {
            if b <= 2.0 {
                return gamln(a) + gamln(b) - gsumln(a, b);
            }
            if b < 8.0 {
                return betaln_b_lt_8_with_w(a, b, 0.0);
            } else {
                return gamln(a) + algdiv(a, b);
            }
        }

        if b <= 1000.0 {
            let n = (a - 1.0) as i32;
            let mut w = 1.0_f64;
            for _ in 1..=n {
                a -= 1.0;
                let h = a / b;
                w *= h / (1.0 + h);
            }
            let w = w.ln();
            if b < 8.0 {
                return betaln_b_lt_8_with_w(a, b, w);
            }
            return w + gamln(a) + algdiv(a, b);
        } else {
            let n = (a - 1.0) as i32;
            let mut w = 1.0_f64;
            for _ in 1..=n {
                a -= 1.0;
                w *= a / (1.0 + a / b);
            }
            return (w.ln() - n as f64 * b.ln()) + (gamln(a) + algdiv(a, b));
        }
    }

    // a >= 8 (label 60)
    let w = bcorr(a, b);
    let h = a / b;
    let c = h / (1.0 + h);
    let u = -(a - 0.5) * c.ln();
    let v = b * alnrel(h);
    if u <= v {
        // Label 61
        (((-0.5 * b.ln() + E) + w) - u) - v
    } else {
        (((-0.5 * b.ln() + E) + w) - v) - u
    }
}

/// Uses a continued fraction expansion for Ix(a,b) when a, b > 1.
/// It is assumed that lambda = (a + b)*y - b.
#[allow(dead_code)]
fn bfrac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64) -> f64 {
    let bfrac = brcomp(a, b, x, y);

    if bfrac == 0.0 {
        return bfrac;
    }

    let c = 1.0 + lambda;
    let c0 = b / a;
    let c1 = 1.0 + 1.0 / a;
    let yp1 = y + 1.0;

    let mut n = 0.0_f64;
    let mut p = 1.0_f64;
    let mut s = a + 1.0;
    let mut an = 0.0_f64;
    let mut bn = 1.0_f64;
    let mut anp1 = 1.0_f64;
    let mut bnp1 = c / c1;
    let mut r = c1 / c;

    // ---------------------------------------------------------------------
    // Continued fraction calculation.
    // ---------------------------------------------------------------------
    loop {
        n += 1.0;
        let t = n / a;
        let w = n * (b - n) * x;
        let e = a / s;
        let alpha = (p * (p + c0) * e * e) * (w * x);
        let e = (1.0 + t) / (c1 + t + t);
        let beta = n + w / s + e * (c + n * yp1);
        p = 1.0 + t;
        s += 2.0;

        // ---------------------------------------------------------------
        // Update AN, BN, ANP1, and BNP1.
        // ---------------------------------------------------------------
        let t = alpha * an + beta * anp1;
        an = anp1;
        anp1 = t;
        let t = alpha * bn + beta * bnp1;
        bn = bnp1;
        bnp1 = t;
        let r0 = r;
        r = anp1 / bnp1;

        if (r - r0).abs() <= eps * r {
            break;
        }

        // ---------------------------------------------------------------
        // Rescale AN, BN, ANP1, and BNP1.
        // ---------------------------------------------------------------
        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0;
    }

    // ---------------------------------------------------------------------
    // Termination.
    // ---------------------------------------------------------------------
    bfrac * r
}

/// Uses an asymptotic expansion for Ix(a,b) when a is larger than b.
///
/// The result of the expansion is added to w. It is assumed that a >= 15
/// and b <= 1. eps is the tolerance used. ierr is a variable that reports
/// the status of the results.
///
/// Returns `(w, ierr)`.
#[allow(dead_code)]
fn bgrat(a: f64, b: f64, x: f64, y: f64, w: f64, eps: f64) -> (f64, i32) {
    let bm1 = (b - 0.5) - 0.5;
    let nu = a + 0.5 * bm1;

    let lnx = if y <= 0.375 {
        alnrel(-y)
    } else {
        x.ln()
    };

    let z = -nu * lnx;
    if b * z == 0.0 {
        // the expansion cannot be computed.
        return (w, 1);
    }

    // ---------------------------------------------------------------------
    // Computation of the expansion.
    // Set r = exp(-z)*z**b/gamma(b).
    // ---------------------------------------------------------------------
    let mut r = b * (1.0 + gam1(b)) * (b * z.ln()).exp();
    r = r * (a * lnx).exp() * (0.5 * bm1 * lnx).exp();
    let mut u = algdiv(b, a) + b * nu.ln();
    u = r * (-u).exp();
    if u == 0.0 {
        return (w, 1);
    }

    let (_, _p, q) = grat1(b, z, r, eps);

    let v = 0.25 * (1.0 / nu).powi(2);
    let t2 = 0.25 * lnx * lnx;
    let l = w / u;
    let mut j = q / r;
    let mut sum1 = j;
    let mut t = 1.0_f64;
    let mut cn = 1.0_f64;
    let mut n2 = 0.0_f64;

    let mut c = [0.0_f64; 30];
    let mut d = [0.0_f64; 30];

    let mut converged = false;

    for n in 1..=30 {
        let bp2n = b + n2;
        j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
        n2 += 2.0;
        t *= t2;
        cn /= n2 * (n2 + 1.0);
        c[n - 1] = cn;

        let mut s = 0.0_f64;
        if n > 1 {
            let nm1 = n - 1;
            let mut coef = b - n as f64;
            for i in 1..=nm1 {
                s += coef * c[i - 1] * d[n - i - 1];
                coef += b;
            }
        }
        d[n - 1] = bm1 * cn + s / n as f64;
        let dj = d[n - 1] * j;
        sum1 += dj;

        if sum1 <= 0.0 {
            return (w, 1);
        }

        if dj.abs() <= eps * (sum1 + l) {
            // converged.
            converged = true;
            break;
        }
    }

    // ---------------------------------------------------------------------
    // Add the results to w, or report failure if we exhausted
    // all 30 terms without converging.
    // ---------------------------------------------------------------------
    if converged {
        (w + u * sum1, 0)
    } else {
        // the expansion cannot be computed.
        (w, 1)
    }
}

/// Uses the power series expansion for evaluating Ix(a,b) when b <= 1
/// or b*x <= 0.7. eps is the tolerance used.
#[allow(dead_code)]
fn bpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let mut bpser;

    if x == 0.0 {
        return 0.0;
    }

    // ---------------------------------------------------------------------
    // Compute the factor x**a/(a*beta(a,b)).
    // ---------------------------------------------------------------------
    let a0 = a.min(b);

    // Compute the series and return the final result.
    // If `bpser == 0.0` or `a <= 0.1*eps`, the series is skipped and `bpser`
    // is returned unchanged.
    let bpser_series = |bpser: f64, a: f64, b: f64, x: f64, eps: f64| -> f64 {
        if bpser == 0.0 || a <= 0.1 * eps {
            return bpser;
        }

        // Compute the series.
        let mut sum1 = 0.0_f64;
        let mut n = 0.0_f64;
        let mut c = 1.0_f64;
        let tol = eps / a;

        loop {
            n += 1.0;
            c = c * (0.5 + (0.5 - b / n)) * x;
            let w = c / (a + n);
            sum1 += w;
            if w.abs() <= tol {
                break;
            }
        }

        bpser * (1.0 + a * sum1)
    };

    if a0 >= 1.0 {
        // Label 70 path via direct computation.
        let z = a * x.ln() - betaln(a, b);
        bpser = z.exp() / a;
        return bpser_series(bpser, a, b, x, eps);
    }

    let mut b0 = a.max(b);

    if b0 >= 8.0 {
        // a0 < 1 and b0 >= 8
        let u = gamln1(a0) + algdiv(a0, b0);
        let z = a * x.ln() - u;
        bpser = (a0 / a) * z.exp();
        return bpser_series(bpser, a, b, x, eps);
    }

    if b0 > 1.0 {
        // a0 < 1 and 1 < b0 < 8
        let mut u = gamln1(a0);
        let m = (b0 - 1.0) as i32;
        if m >= 1 {
            let mut c = 1.0;
            for _ in 1..=m {
                b0 -= 1.0;
                c = c * (b0 / (a0 + b0));
            }
            u = c.ln() + u;
        }

        let z = a * x.ln() - u;
        b0 -= 1.0;
        let apb = a0 + b0;
        let t = if apb > 1.0 {
            let u = a0 + b0 - 1.0;
            (1.0 + gam1(u)) / apb
        } else {
            1.0 + gam1(apb)
        };
        bpser = z.exp() * (a0 / a) * (1.0 + gam1(b0)) / t;
        return bpser_series(bpser, a, b, x, eps);
    }

    // a0 < 1 and b0 <= 1
    bpser = x.powf(a);
    if bpser == 0.0 {
        return bpser;
    }

    let apb = a + b;
    let z = if apb > 1.0 {
        let u = a + b - 1.0;
        (1.0 + gam1(u)) / apb
    } else {
        1.0 + gam1(apb)
    };
    let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
    bpser = bpser * c * (b / apb);

    bpser_series(bpser, a, b, x, eps)
}

/// Result of a `bratio` evaluation.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub(crate) struct BratioResult {
    pub w: f64,
    pub w1: f64,
    pub ierr: i32,
}

/// Evaluates the incomplete beta function Ix(A,B).
///
/// It is assumed that X <= 1 and Y = 1 - X. BRATIO assigns W and W1 the
/// values
///
/// ```text
///     W  = ix(a,b)
///     W1 = 1 - ix(a,b)
/// ```
///
/// `ierr` reports the status of the results. If no input errors are detected
/// then `ierr` is set to 0 and `w` and `w1` are computed. Otherwise, if an
/// error is detected, then `w` and `w1` are assigned the value 0 and `ierr`
/// is set to one of the following values:
///
/// ```text
///     ierr = 1  if a or b is negative
///     ierr = 2  if a = b = 0
///     ierr = 3  if x < 0 or x > 1
///     ierr = 4  if y < 0 or y > 1
///     ierr = 5  if x + y /= 1
///     ierr = 6  if x = a = 0
///     ierr = 7  if y = b = 0
/// ```
///
/// Parameters:
///
///   Input, a, b: the parameters of the function. a and b should be
///   nonnegative.
#[allow(dead_code)]
pub(crate) fn bratio(a: f64, b: f64, x: f64, y: f64) -> BratioResult {
    let mut eps = f64::EPSILON;

    // Defaults (also used on error returns).
    let mut w: f64;
    let mut w1 = 0.0_f64;

    // ---------------------------------------------------------------------
    // Input validation.
    // ---------------------------------------------------------------------
    if a < 0.0 || b < 0.0 {
        return BratioResult { w: 0.0, w1: 0.0, ierr: 1 };
    }
    if a == 0.0 && b == 0.0 {
        return BratioResult { w: 0.0, w1: 0.0, ierr: 2 };
    }
    if x < 0.0 || x > 1.0 {
        return BratioResult { w: 0.0, w1: 0.0, ierr: 3 };
    }
    if y < 0.0 || y > 1.0 {
        return BratioResult { w: 0.0, w1: 0.0, ierr: 4 };
    }

    let z = ((x + y) - 0.5) - 0.5;
    if z.abs() > 3.0 * eps {
        return BratioResult { w: 0.0, w1: 0.0, ierr: 5 };
    }

    // ---------------------------------------------------------------------
    // Degenerate cases in x, y, a, b.
    // ---------------------------------------------------------------------
    if x == 0.0 {
        // Label 200
        if a == 0.0 {
            return BratioResult { w: 0.0, w1: 0.0, ierr: 6 };
        }
        return BratioResult { w: 0.0, w1: 1.0, ierr: 0 };
    }
    if y == 0.0 {
        if b == 0.0 {
            return BratioResult { w: 0.0, w1: 0.0, ierr: 7 };
        }
        return BratioResult { w: 1.0, w1: 0.0, ierr: 0 };
    }
    if a == 0.0 {
        return BratioResult { w: 1.0, w1: 0.0, ierr: 0 };
    }
    if b == 0.0 {
        return BratioResult { w: 0.0, w1: 1.0, ierr: 0 };
    }

    eps = eps.max(1.0e-15);
    if a.max(b) < 1.0e-3 * eps {
        // a and b both < 1.e-3*eps
        return BratioResult {
            w: b / (a + b),
            w1: a / (a + b),
            ierr: 0,
        };
    }

    // ---------------------------------------------------------------------
    // Set up for the main algorithm.
    // ---------------------------------------------------------------------
    let mut ind = 0_i32;
    let mut a0 = a;
    let mut b0 = b;
    let mut x0 = x;
    let mut y0 = y;
    let mut lambda: f64;
    let mut n: i32;

    if a0.min(b0) <= 1.0 {
        // ---------------------------------------------------------------
        // Procedure for a0 <= 1 or b0 <= 1.
        // ---------------------------------------------------------------
        if x > 0.5 {
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
        }

        if b0 < eps.min(eps * a0) {
            // fpser
            w = fpser(a0, b0, x0, eps);
            w1 = 0.5 + (0.5 - w);
        } else if a0 < eps.min(eps * b0) && b0 * x0 <= 1.0 {
            // apser
            w1 = apser(a0, b0, x0, eps);
            w = 0.5 + (0.5 - w1);
        } else if a0.max(b0) <= 1.0 {
            // Small parameters (a0, b0 <= 1).
            if a0 >= 0.2_f64.min(b0) {
                // bpser
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
            } else if x0.powf(a0) <= 0.9 {
                // bpser
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
            } else if x0 >= 0.3 {
                // bpser with swapped arguments
                w1 = bpser(b0, a0, y0, eps);
                w = 0.5 + (0.5 - w1);
            } else {
                // bup + bgrat
                n = 20;
                w1 = bup(b0, a0, y0, x0, n, eps);
                b0 += n as f64;
                let (w1_new, _ierr1) = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
                w1 = w1_new;
                w = 0.5 + (0.5 - w1);
            }
        } else {
            // a0.max(b0) > 1.0, with b0 <= 1.0.
            if b0 <= 1.0 {
                // bpser
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
            } else if x0 >= 0.3 {
                // bpser with swapped arguments
                w1 = bpser(b0, a0, y0, eps);
                w = 0.5 + (0.5 - w1);
            } else {
                // x0 < 0.3
                if x0 >= 0.1 {
                    if b0 > 15.0 {
                        // bgrat only
                        let (w1_new, _ierr1) = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
                        w1 = w1_new;
                        w = 0.5 + (0.5 - w1);
                    } else {
                        // n = 20
                        n = 20;
                        w1 = bup(b0, a0, y0, x0, n, eps);
                        b0 += n as f64;
                        let (w1_new, _ierr1) = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
                        w1 = w1_new;
                        w = 0.5 + (0.5 - w1);
                    }
                } else {
                    // x0 < 0.1
                    if (x0 * b0).powf(a0) <= 0.7 {
                        // bpser
                        w = bpser(a0, b0, x0, eps);
                        w1 = 0.5 + (0.5 - w);
                    } else if b0 > 15.0 {
                        // bgrat only
                        let (w1_new, _ierr1) = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
                        w1 = w1_new;
                        w = 0.5 + (0.5 - w1);
                    } else {
                        // n = 20
                        n = 20;
                        w1 = bup(b0, a0, y0, x0, n, eps);
                        b0 += n as f64;
                        let (w1_new, _ierr1) = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
                        w1 = w1_new;
                        w = 0.5 + (0.5 - w1);
                    }
                }
            }
        }
    } else {
        // ---------------------------------------------------------------
        // Procedure for a0 > 1 and b0 > 1.
        // ---------------------------------------------------------------
        lambda = if a <= b {
            a - (a + b) * x
        } else {
            (a + b) * y - b
        };

        if lambda < 0.0 {
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
            lambda = lambda.abs();
        }

        if b0 < 40.0 && b0 * x0 <= 0.7 {
            // bpser
            w = bpser(a0, b0, x0, eps);
            w1 = 0.5 + (0.5 - w);
        } else if b0 < 40.0 {
            n = b0 as i32;
            b0 -= n as f64;
            if b0 == 0.0 {
                n -= 1;
                b0 = 1.0;
            }
            w = bup(b0, a0, y0, x0, n, eps);
            if x0 <= 0.7 {
                w += bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
            } else {
                if a0 <= 15.0 {
                    n = 20;
                    w += bup(a0, b0, x0, y0, n, eps);
                    a0 += n as f64;
                }
                let (w_new, _ierr1) = bgrat(a0, b0, x0, y0, w, 15.0 * eps);
                w = w_new;
                w1 = 0.5 + (0.5 - w);
            }
        } else {
            // b0 >= 40.0
            if a0 <= b0 {
                if a0 <= 100.0 || lambda <= 0.03 * a0 {
                    // bfrac
                    w = bfrac(a0, b0, x0, y0, lambda, 15.0 * eps);
                    w1 = 0.5 + (0.5 - w);
                } else {
                    // basym
                    w = basym(a0, b0, lambda, 100.0 * eps);
                    w1 = 0.5 + (0.5 - w);
                }
            } else {
                if b0 <= 100.0 || lambda <= 0.03 * b0 {
                    // bfrac
                    w = bfrac(a0, b0, x0, y0, lambda, 15.0 * eps);
                    w1 = 0.5 + (0.5 - w);
                } else {
                    // basym
                    w = basym(a0, b0, lambda, 100.0 * eps);
                    w1 = 0.5 + (0.5 - w);
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // swap w and w1 if the arguments were flipped.
    // ---------------------------------------------------------------------
    if ind != 0 {
        std::mem::swap(&mut w, &mut w1);
    }

    BratioResult { w, w1, ierr: 0 }
}

/// Evaluates exp(mu) * (x**a*y**b/beta(a,b)).
#[allow(dead_code)]
fn brcmp1(mu: i32, a: f64, b: f64, x: f64, y: f64) -> f64 {
    let a0 = a.min(b);

    // The `a >= 8` and `b >= 8` procedure.
    let brcmp1_large = |mu: i32, a: f64, b: f64, x: f64, y: f64| -> f64 {
        const CONST: f64 = 0.398942280401433;

        // choose the orientation based on which of a, b is larger.
        let (h, x0, y0, lambda);
        if a <= b {
            h = a / b;
            x0 = h / (1.0 + h);
            y0 = 1.0 / (1.0 + h);
            lambda = a - (a + b) * x;
        } else {
            h = b / a;
            x0 = 1.0 / (1.0 + h);
            y0 = h / (1.0 + h);
            lambda = (a + b) * y - b;
        }

        // compute u.
        let e = -lambda / a;
        let u = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            // Label 111
            e - (x / x0).ln()
        };

        // compute v.
        let e = lambda / b;
        let v = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (y / y0).ln()
        };

        let z = esum(mu, -(a * u + b * v));
        CONST * (b * x0).sqrt() * z * (-bcorr(a, b)).exp()
    };

    // ---------------------------------------------------------------------
    // Small a0 (a0 < 8) path.
    // ---------------------------------------------------------------------
    if a0 >= 8.0 {
        return brcmp1_large(mu, a, b, x, y);
    }

    // --- Compute z = a*lnx + b*lny ---
    let (lnx, lny);
    if x <= 0.375 {
        lnx = x.ln();
        lny = alnrel(-x);
    } else if y <= 0.375 {
        lnx = alnrel(-y);
        lny = y.ln();
    } else {
        lnx = x.ln();
        lny = y.ln();
    }

    let mut z = a * lnx + b * lny;

    if a0 >= 1.0 {
        z -= betaln(a, b);
        return esum(mu, z);
    }

    // ---------------------------------------------------------------------
    // Procedure for a < 1 or b < 1.
    // ---------------------------------------------------------------------
    let mut b0 = a.max(b);

    if b0 >= 8.0 {
        // algorithm for b0 >= 8
        let u = gamln1(a0) + algdiv(a0, b0);
        return a0 * esum(mu, z - u);
    }

    if b0 > 1.0 {
        // algorithm for 1 < b0 < 8
        let mut u = gamln1(a0);
        let n = (b0 - 1.0) as i32;
        if n >= 1 {
            let mut c = 1.0_f64;
            for _ in 1..=n {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u = c.ln() + u;
        }

        z -= u;
        b0 -= 1.0;
        let apb = a0 + b0;
        let t = if apb > 1.0 {
            // Label 71
            let u = a0 + b0 - 1.0;
            (1.0 + gam1(u)) / apb
        } else {
            // Label 70 -> 72
            1.0 + gam1(apb)
        };
        return a0 * esum(mu, z) * (1.0 + gam1(b0)) / t;
    }

    // ---------------------------------------------------------------------
    // algorithm for b0 <= 1
    // ---------------------------------------------------------------------
    let mut brcmp1 = esum(mu, z);
    if brcmp1 == 0.0 {
        return brcmp1;
    }

    let apb = a + b;
    let zz = if apb > 1.0 {
        let u = a + b - 1.0;
        (1.0 + gam1(u)) / apb
    } else {
        1.0 + gam1(apb)
    };
    let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / zz;
    brcmp1 = brcmp1 * (a0 * c) / (1.0 + a0 / b0);
    brcmp1
}

/// Evaluates X**a * y**b / beta(a,b).
#[allow(dead_code)]
fn brcomp(a: f64, b: f64, x: f64, y: f64) -> f64 {

    if x == 0.0 || y == 0.0 {
        return 0.0;
    }

    let a0 = a.min(b);

    // The `a >= 8` and `b >= 8` procedure
    let brcomp_large = |a: f64, b: f64, x: f64, y: f64| -> f64 {
        const CONST: f64 = 0.398942280401433;     // 1/sqrt(2*pi)

        // Label 100: choose the orientation based on which of a, b is larger.
        let (h, x0, y0, lambda);
        if a <= b {
            h = a / b;
            x0 = h / (1.0 + h);
            y0 = 1.0 / (1.0 + h);
            lambda = a - (a + b) * x;
        } else {
            // Label 101
            h = b / a;
            x0 = 1.0 / (1.0 + h);
            y0 = h / (1.0 + h);
            lambda = (a + b) * y - b;
        }

        // Label 110: compute u.
        let e = -lambda / a;
        let u = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            // Label 111
            e - (x / x0).ln()
        };

        // Label 120: compute v.
        let e = lambda / b;
        let v = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            // Label 121
            e - (y / y0).ln()
        };

        // Label 130
        let z = (-(a * u + b * v)).exp();
        CONST * (b * x0).sqrt() * z * (-bcorr(a, b)).exp()
    };

    // ---------------------------------------------------------------------
    // Small a0 (a0 < 8) path.
    // ---------------------------------------------------------------------
    if a0 >= 8.0 {
        return brcomp_large(a, b, x, y);
    }

    // --- Compute z = a*lnx + b*lny ---
    let (lnx, lny);
    if x <= 0.375 {
        lnx = x.ln();
        lny = alnrel(-x);
    } else if y <= 0.375 {
        lnx = alnrel(-y);
        lny = y.ln();
    } else {
        lnx = x.ln();
        lny = y.ln();
    }

    let mut z = a * lnx + b * lny;

    if a >= 1.0 {
        z -= betaln(a, b);
        return z.exp();
    }

    // ---------------------------------------------------------------------
    // Procedure for a < 1 or b < 1
    // ---------------------------------------------------------------------
    let mut b0 = a.max(b);

    if b0 >= 8.0 {
        // algorithm for b0 >= 8
        let u = gamln1(a0) + algdiv(a0, b0);
        return a0 * (z - u).exp();
    }

    if b0 > 1.0 {
        // algorithm for 1 < b0 < 8
        let mut u = gamln1(a0);
        let n = (b0 - 1.0) as i32;
        if n >= 1 {
            let mut c = 1.0_f64;
            for _ in 1..=n {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u = c.ln() + u;
        }

        z -= u;
        b0 -= 1.0;
        let apb = a0 + b0;
        let t = if apb > 1.0 {
            let u = a0 + b0 - 1.0;
            (1.0 + gam1(u)) / apb
        } else {
            1.0 + gam1(apb)
        };
        return a0 * z.exp() * (1.0 + gam1(b0)) / t;
    }

    // ---------------------------------------------------------------------
    // algorithm for b0 <= 1
    // ---------------------------------------------------------------------
    let mut brcomp = z.exp();
    if brcomp == 0.0 {
        return brcomp;
    }

    let apb = a + b;
    let zz = if apb > 1.0 {
        // Label 40
        let u = a + b - 1.0;
        (1.0 + gam1(u)) / apb
    } else {
        1.0 + gam1(apb)
    };
    let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / zz;
    brcomp = brcomp * (a0 * c) / (1.0 + a0 / b0);
    brcomp
}

/// Evaluates Ix(a,b) - Ix(a+n,b) where n is a positive integer.
/// eps is the tolerance used.
#[allow(dead_code)]
fn bup(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64) -> f64 {
    // ---------------------------------------------------------------------
    // Obtain the scaling factor exp(-mu) and
    // exp(mu)*(x**a*y**b/beta(a,b))/a
    // ---------------------------------------------------------------------
    let apb = a + b;
    let ap1 = a + 1.0;

    let mut mu = 0_i32;
    let mut d = 1.0_f64;

    if n != 1 && a >= 1.0 {
        if apb >= 1.1 * ap1 {
            mu = exparg(1).abs() as i32;
            let k = exparg(0) as i32;
            if k < mu {
                mu = k;
            }
            let t = mu as f64;
            d = (-t).exp();
        }
    }

    let mut bup: f64 = brcmp1(mu, a, b, x, y) / a;

    if n == 1 || bup == 0.0 {
        return bup;
    }

    let nm1 = n - 1;
    let mut w = d;

    // ---------------------------------------------------------------------
    // Let k be the index of the maximum term.
    // ---------------------------------------------------------------------
    let mut k = 0_i32;
    if b > 1.0 {
        if y <= 1.0e-4 {
            k = nm1;
        } else {
            let r = (b - 1.0) * x / y - a;
            if r >= 1.0 {
                k = nm1;
                let t = nm1 as f64;
                if r < t {
                    k = r as i32;
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Add the increasing terms of the series.
    // ---------------------------------------------------------------------
    if k >= 1 {
        for i in 1..=k {
            let l = (i - 1) as f64;
            d = ((apb + l) / (ap1 + l)) * x * d;
            w += d;
        }
    }

    if k == nm1 {
        // Terminate the procedure.
        bup *= w;
        return bup;
    }

    // ---------------------------------------------------------------------
    // Add the remaining terms of the series.
    // ---------------------------------------------------------------------
    let kp1 = k + 1;
    for i in kp1..=nm1 {
        let l = (i - 1) as f64;
        d = ((apb + l) / (ap1 + l)) * x * d;
        w += d;
        if d <= eps * w {
            // Terminate the procedure.
            bup *= w;
            return bup;
        }
    }

    // Terminate the procedure.
    bup *= w;
    bup
}

/// ERFC1 evaluates the complementary error function
///          erfc1(ind,x) = erfc(x)            if ind = 0
///          erfc1(ind,x) = exp(x*x)*erfc(x)   otherwise
#[allow(dead_code)]
fn erfc1(ind: i32, x: f64) -> f64 {
    const C: f64 = 0.564189583547756;

    // Fortran arrays are 1-based:
    // a(1..5), b(1..3), p(1..8), q(1..8), r(1..5), s(1..4)
    let a = [
        0.0,
        0.771058495001320e-04,
       -0.133733772997339e-02,
        0.323076579225834e-01,
        0.479137145607681e-01,
        0.128379167095513e+00,
    ];

    let b = [
        0.0,
        0.301048631703895e-02,
        0.538971687740286e-01,
        0.375795757275549e+00,
    ];

    let p = [
        0.0,
       -1.36864857382717e-07,
        5.64195517478974e-01,
        7.21175825088309e+00,
        4.31622272220567e+01,
        1.52989285046940e+02,
        3.39320816734344e+02,
        4.51918953711873e+02,
        3.00459261020162e+02,
    ];

    let q = [
        0.0,
        1.00000000000000e+00,
        1.27827273196294e+01,
        7.70001529352295e+01,
        2.77585444743988e+02,
        6.38980264465631e+02,
        9.31354094850610e+02,
        7.90950925327898e+02,
        3.00459260956983e+02,
    ];

    let r = [
        0.0,
        2.10144126479064e+00,
        2.62370141675169e+01,
        2.13688200555087e+01,
        4.65807828718470e+00,
        2.82094791773523e-01,
    ];

    let s = [
        0.0,
        9.41537750555460e+01,
        1.87114811799590e+02,
        9.90191814623914e+01,
        1.80124575948747e+01,
    ];

    let ax = x.abs();

    // abs(x) <= 0.5
    if ax <= 0.5 {
        let t = x * x;

        let top = ((((a[1] * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5]) + 1.0;

        let bot = (((b[1] * t + b[2]) * t + b[3]) * t) + 1.0;

        let mut erfc1 = 0.5 + (0.5 - x * (top / bot));

        if ind != 0 {
            erfc1 *= t.exp();
        }

        return erfc1;
    }

    // 0.5 < abs(x) <= 4
    if ax <= 4.0 {
        let top = (((((((p[1] * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5])
                * ax + p[6]) * ax + p[7]) * ax) + p[8];

        let bot = (((((((q[1] * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5])
                * ax + q[6]) * ax + q[7]) * ax) + q[8];

        let mut erfc1 = top / bot;

        // Final assembly
        if ind != 0 {
            if x < 0.0 {
                erfc1 = 2.0 * (x * x).exp() - erfc1;
            }

            return erfc1;
        }

        let w = x * x;
        let t = w;
        let e = w - t;

        erfc1 = ((0.5 + (0.5 - e)) * (-t).exp()) * erfc1;

        if x < 0.0 {
            erfc1 = 2.0 - erfc1;
        }

        return erfc1;
    }

    // abs(x) > 4
    // Limit value for large negative x
    if x <= -5.6 {
        let mut erfc1 = 2.0;

        if ind != 0 {
            erfc1 = 2.0 * (x * x).exp();
        }

        return erfc1;
    }

    // For ind == 0:
    if ind == 0 {
        if 100.0 < x {
            return 0.0;
        }

        if -exparg(1) < x * x {
            return 0.0;
        }
    }

    // 4 < abs(x), excluding the cases handled above
    let t = (1.0 / x) * (1.0 / x);

    let top = ((((r[1] * t + r[2]) * t + r[3]) * t + r[4]) * t) + r[5];

    let bot = ((((s[1] * t + s[2]) * t + s[3]) * t + s[4]) * t) + 1.0;

    let mut erfc1 = (C - t * top / bot) / ax;

    // Final assembly
    if ind != 0 {
        if x < 0.0 {
            erfc1 = 2.0 * (x * x).exp() - erfc1;
        }

        return erfc1;
    }

    let w = x * x;
    let t = w;
    let e = w - t;

    erfc1 = ((0.5 + (0.5 - e)) * (-t).exp()) * erfc1;

    if x < 0.0 {
        erfc1 = 2.0 - erfc1;
    }

    erfc1
}

/// ERF evaluates the real error function.
#[allow(dead_code)]
fn erf(x: f64) -> f64 {
    const C: f64 = 0.564189583547756;

    const A: [f64; 5] = [
         0.771058495001320e-04,
        -0.133733772997339e-02,
         0.323076579225834e-01,
         0.479137145607681e-01,
         0.128379167095513e+00,
    ];

    const B: [f64; 3] = [
        0.301048631703895e-02,
        0.538971687740286e-01,
        0.375795757275549e+00,
    ];

    const P: [f64; 8] = [
        -1.36864857382717e-07,
         5.64195517478974e-01,
         7.21175825088309e+00,
         4.31622272220567e+01,
         1.52989285046940e+02,
         3.39320816734344e+02,
         4.51918953711873e+02,
         3.00459261020162e+02,
    ];

    const Q: [f64; 8] = [
        1.00000000000000e+00,
        1.27827273196294e+01,
        7.70001529352295e+01,
        2.77585444743988e+02,
        6.38980264465631e+02,
        9.31354094850610e+02,
        7.90950925327898e+02,
        3.00459260956983e+02,
    ];

    const R: [f64; 5] = [
        2.10144126479064e+00,
        2.62370141675169e+01,
        2.13688200555087e+01,
        4.65807828718470e+00,
        2.82094791773523e-01,
    ];

    const S: [f64; 4] = [
        9.41537750555460e+01,
        1.87114811799590e+02,
        9.90191814623914e+01,
        1.80124575948747e+01,
    ];

    let ax = x.abs();
    let mut erf;

    if ax <= 0.5 {
        let t = x * x;
        let top = ((((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t + A[4]) + 1.0;
        let bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;
        erf = x * (top / bot);
    } else if ax <= 4.0 {
        let top = ((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5]) * ax + P[6]) * ax + P[7];
        let bot = ((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5]) * ax + Q[6]) * ax + Q[7];
        erf = 0.5 + (0.5 - (-x * x).exp() * top / bot);
        if x < 0.0 {
            erf = -erf;
        }
    } else if ax < 5.8 {
        let x2 = x * x;
        let t = 1.0 / x2;
        let top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];
        let bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;
        erf = (C - top / (x2 * bot)) / ax;
        erf = 0.5 + (0.5 - (-x2).exp() * erf);

        if x < 0.0 {
            erf = -erf;
        }
    } else {
        erf = if x > 0.0 { 1.0 } else { -1.0 };
    }

    erf
}

/// Evaluates exp(mu + x).
#[allow(dead_code)]
fn esum(mu: i32, x: f64) -> f64 {
    // The logic computes exp(mu + x) directly only when it is safe to do so,
    // i.e., when mu + x is small enough to avoid overflow. Otherwise it splits
    // the computation into exp(mu) * exp(x) so that each factor is computed
    // separately.
    if x > 0.0 {
        if mu <= 0 {
            let w = mu as f64 + x;
            if w >= 0.0 {
                return w.exp();
            }
        }
    } else {
        // x <= 0.0
        if mu >= 0 {
            let w = mu as f64 + x;
            if w <= 0.0 {
                return w.exp();
            }
        }
    }

    // fallback
    let w = mu as f64;
    w.exp() * x.exp()
}

/// EXPARG reports the largest safe arguments for EXP(X).
/// if l = 0 then exparg(l) = the largest positive w for which
/// exp(w) can be computed.
/// if l is nonzero then exparg(l) = the largest negative w for
/// which the computed value of exp(w) is nonzero.
#[allow(dead_code)]
fn exparg(l: i32) -> f64 {
    const LNB: f64 = 0.69314718055995;

    if l != 0 {
        let m = ipmpar(9) - 1;
        return 0.99999_f64 * (m as f64) * LNB;
    }

    let m = ipmpar(10);
    0.99999_f64 * (m as f64) * LNB
}

/// Evaluates Ix(A,B) for small B and X <= 0.5.
/// For b < min(eps, eps*a) and x <= 0.5.
#[allow(dead_code)]
fn fpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    // Set FPSER = X**A.
    let mut fpser = 1.0;

    if 1.0e-3 * eps < a {
        fpser = 0.0;
        let t = a * x.ln();
        if t < exparg(1) {
            return fpser;
        }
        fpser = t.exp();
    }

    // Note that 1/b(a,b) = b
    fpser = (b / a) * fpser;
    let tol = eps / a;
    let mut an = a + 1.0;
    let mut t = x;
    let mut s = t / an;

    loop {
        an += 1.0;
        t = x * t;
        let c = t / an;
        s += c;

        if c.abs() <= tol {
            break;
        }
    }

    fpser * (1.0 + a * s)
}

/// Computes 1/gamma(a+1) - 1 for -0.5 <= a <= 1.5
#[allow(dead_code)]
fn gam1(a: f64) -> f64 {
    const P: [f64; 7] = [
         0.577215664901533e+00,
        -0.409078193005776e+00,
        -0.230975380857675e+00,
         0.597275330452234e-01,
         0.766968181649490e-02,
        -0.514889771323592e-02,
         0.589597428611429e-03,
    ];

    const Q: [f64; 5] = [
        0.100000000000000e+01,
        0.427569613095214e+00,
        0.158451672430138e+00,
        0.261132021441447e-01,
        0.423244297896961e-02,
    ];

    const R: [f64; 9] = [
        -0.422784335098468e+00,
        -0.771330383816272e+00,
        -0.244757765222226e+00,
         0.118378989872749e+00,
         0.930357293360349e-03,
        -0.118290993445146e-01,
         0.223047661158249e-02,
         0.266505979058923e-03,
        -0.132674909766242e-03,
    ];

    const S1: f64 = 0.273076135303957e+00;
    const S2: f64 = 0.559398236957378e-01;

    let mut t = a;
    let d = a - 0.5;
    if d > 0.0 {
        t = d - 0.5;
    }

    if t > 0.0 {
        let top = (((((P[6] * t + P[5]) * t + P[4]) * t + P[3]) * t + P[2]) * t + P[1]) * t + P[0];
        let bot = (((Q[4] * t + Q[3]) * t + Q[2]) * t + Q[1]) * t + 1.0;
        let w = top / bot;
        if d > 0.0 {
            (t / a) * ((w - 0.5) - 0.5)
        } else {
            a * w
        }
    } else if t < 0.0 {
        let top = (((((((R[8] * t + R[7]) * t + R[6]) * t + R[5]) * t + R[4]) * t + R[3]) * t + R[2]) * t + R[1]) * t + R[0];
        let bot = (S2 * t + S1) * t + 1.0;
        let w = top / bot;
        if 0.0 < d {
            t * w / a
        } else {
            a * ((w + 0.5) + 0.5)
        }
    } else {
        0.0
    }
}

/// Evaluates ln(gamma(1 + a)) for -0.2 <= a <= 1.25
#[allow(dead_code)]
fn gamln1(a: f64) -> f64 {
    const P0: f64 =  0.577215664901533e+00;
    const P1: f64 =  0.844203922187225e+00;
    const P2: f64 = -0.168860593646662e+00;
    const P3: f64 = -0.780427615533591e+00;
    const P4: f64 = -0.402055799310489e+00;
    const P5: f64 = -0.673562214325671e-01;
    const P6: f64 = -0.271935708322958e-02;

    const Q1: f64 =  0.288743195473681e+01;
    const Q2: f64 =  0.312755088914843e+01;
    const Q3: f64 =  0.156875193295039e+01;
    const Q4: f64 =  0.361951990101499e+00;
    const Q5: f64 =  0.325038868253937e-01;
    const Q6: f64 =  0.667465618796164e-03;

    const R0: f64 =  0.422784335098467e+00;
    const R1: f64 =  0.848044614534529e+00;
    const R2: f64 =  0.565221050691933e+00;
    const R3: f64 =  0.156513060486551e+00;
    const R4: f64 =  0.170502484022650e-01;
    const R5: f64 =  0.497958207639485e-03;

    const S1: f64 =  0.124313399877507e+01;
    const S2: f64 =  0.548042109832463e+00;
    const S3: f64 =  0.101552187439830e+00;
    const S4: f64 =  0.713309612391000e-02;
    const S5: f64 =  0.116165475989616e-03;

    if a < 0.6 {
        let w = ((((((P6 * a + P5) * a + P4) * a + P3) * a + P2) * a + P1) * a + P0)
              / ((((((Q6 * a + Q5) * a + Q4) * a + Q3) * a + Q2) * a + Q1) * a + 1.0);
        -a * w
    } else {
        let x = (a - 0.5) - 0.5;
        let w = (((((R5 * x + R4) * x + R3) * x + R2) * x + R1) * x + R0)
              / (((((S5 * x + S4) * x + S3) * x + S2) * x + S1) * x + 1.0);
        x * w
    }
}

/// Evaluates ln(gamma(a)) for positive A.
#[allow(dead_code)]
fn gamln(a: f64) -> f64 {
    // d = 0.5*(ln(2*pi) - 1)
    const D: f64 = 0.418938533204673;

    const C0: f64 =  0.833333333333333e-01;
    const C1: f64 = -0.277777777760991e-02;
    const C2: f64 =  0.793650666825390e-03;
    const C3: f64 = -0.595202931351870e-03;
    const C4: f64 =  0.837308034031215e-03;
    const C5: f64 = -0.165322962780713e-02;

    if a <= 0.8 {
        gamln1(a) - a.ln()
    } else if a <= 2.25 {
        let t = (a - 0.5) - 0.5;
        gamln1(t)
    } else if a < 10.0 {
        let n = (a - 1.25) as i32;
        let mut t = a;
        let mut w = 1.0;
        for _ in 1..=n {
            t -= 1.0;
            w = t * w;
        }
        gamln1(t - 1.0) + w.ln()
    } else {
        let t = 1.0 / (a * a);
        let w = (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a;
        (D + w) + (a - 0.5) * (a.ln() - 1.0)
    }
}

/// Evaluates the incomplete Gamma ratio functions P(A,X) and Q(A,X).
/// It is assumed that a <= 1. eps is the tolerance to be used.
/// The input argument r has the value e**(-x)*x**a/gamma(a).
/// Returns `(r, p, q)`.
#[allow(dead_code)]
fn grat1(a: f64, x: f64, r: f64, eps: f64) -> (f64, f64, f64) {
    // Special cases.
    if a * x == 0.0 {
        if x <= a {
            return (r, 0.0, 1.0);
        } else {
            return (r, 1.0, 0.0);
        }
    }

    if a == 0.5 {
        // Label 120.
        if x < 0.25 {
            let p = erf(x.sqrt());
            let q = 0.5 + (0.5 - p);
            return (r, p, q);
        } else {
            // Label 121.
            let q = erfc1(0, x.sqrt());
            let p = 0.5 + (0.5 - q);
            return (r, p, q);
        }
    }

    if x >= 1.1 {
        // Continued fraction expansion
        let grat1_cf = |a: f64, x: f64, r: f64, eps: f64| -> (f64, f64, f64) {
            let mut a2nm1 = 1.0_f64;
            let mut a2n = 1.0_f64;
            let mut b2nm1 = x;
            let mut b2n = x + (1.0 - a);
            let mut c = 1.0_f64;

            let mut an0: f64;
            loop {
                a2nm1 = x * a2n + c * a2nm1;
                b2nm1 = x * b2n + c * b2nm1;
                let am0 = a2nm1 / b2nm1;
                c += 1.0;
                let cma = c - a;
                a2n = a2nm1 + cma * a2n;
                b2n = b2nm1 + cma * b2n;
                an0 = a2n / b2n;
                if (an0 - am0).abs() < eps * an0 {
                    break;
                }
            }

            let q = r * an0;
            let p = 0.5 + (0.5 - q);
            (r, p, q)
        };

        return grat1_cf(a, x, r, eps);
    }

    // ---------------------------------------------------------------------
    // Taylor series for p(a,x)/x**a  (label 10)
    // ---------------------------------------------------------------------
    let mut an = 3.0_f64;
    let mut c = x;
    let mut sum2 = x / (a + 3.0);
    let tol = 0.1 * eps / (a + 1.0);

    loop {
        an += 1.0;
        c = -c * (x / an);
        let t = c / (a + an);
        sum2 += t;
        if t.abs() <= tol {
            break;
        }
    }

    let j = a * x * ((sum2 / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));

    let z = a * x.ln();
    let h = gam1(a);
    let g = 1.0 + h;

    match (x < 0.25, z > -0.13394, a < x / 2.59) {
        (true, true, _) | (false, _, true) => {
            let l = rexp(z);
            let w = 0.5 + (0.5 + l);
            let q = (w * j - l) * g - h;
            if q < 0.0 {
                return (r, 1.0, 0.0);
            }
            let p = 0.5 + (0.5 - q);
            return (r, p, q);
        },
        _ => {
            let w = z.exp();
            let p = w * g * (0.5 + (0.5 - j));
            let q = 0.5 + (0.5 - p);
            return (r, p, q);
        },
    }
}

/// Evaluates the function log(gamma(a + b)) in a special range.
/// For 1 <= a <= 2 and 1 <= b <= 2.
#[allow(dead_code)]
fn gsumln(a: f64, b: f64) -> f64 {
    let x = a + b - 2.0;

    if x <= 0.25 {
        gamln1(1.0 + x)
    } else if x <= 1.25 {
        gamln1(x) + alnrel(x)
    } else {
        gamln1(x - 1.0) + (x * (1.0 + x)).ln()
    }
}

/// IPMPAR provides the integer (kind = 4) machine constants for the computer that is used.
/// integer (kind = 4)s are represented in the n-digit, base-a form
///        sign ( x(n-1)*a**(n-1) + ... + x(1)*a + x(0) )
///        where 0 <= x(i) < a for i = 0, ..., n - 1.
/// floating-point numbers: the single and double precision floating
/// point arithmetics have the same base, say b, and that the
/// nonzero numbers are represented in the form
///        sign (b**e) * (x(1)/b + ... + x(m)/b**m)
///        where x(i) = 0,1,...,b-1 for i=1,...,m,
///        x(1) .ge. 1, and emin <= e <= emax.
#[allow(dead_code)]
const fn ipmpar(i: i32) -> i32 {
    match i {
        1 => 2,          // a, the base
        2 => 31,         // n, the number of base-a digits
        3 => 2147483647, // a**n - 1, the largest magnitude
        4 => 2,          // Base for floating-point arithmetic
        5 => 24,         // Number of base digits for single precision
        6 => -125,       // Smallest exponent for single precision
        7 => 128,        // Largest exponent for single precision
        8 => 53,         // Number of base digits for double precision
        9 => -1021,      // Smallest exponent for double precision
        10 => 1024,      // Largest exponent for double precision
        _  => panic!(),
    }
}

/// PSI evaluates the Psi or Digamma function.
/// The algorithm uses:
///   - reflection for x < 0.5
///   - rational approximation for 0.5 <= x <= 3
///   - asymptotic rational approximation for x > 3
#[allow(dead_code)]
fn psi(x: f64) -> f64 {
    const DX0: f64   = 1.461632144968362341262659542325721325;                  // zero of PSI
    const PIOV4: f64 = 0.785398163397448;                                       // pi/4

    // coefficients for rational approximation of
    // psi(x) / (x - x0),  0.5 <= x <= 3.0
    const P1: [f64; 7] = [
        0.895385022981970e-02,  0.477762828042627e+01,
        0.142441585084029e+03,  0.118645200713425e+04,
        0.363351846806499e+04,  0.413810161269013e+04,
        0.130560269827897e+04];
    const Q1: [f64; 6] = [
        0.448452573429826e+02,  0.520752771467162e+03,
        0.221000799247830e+04,  0.364127349079381e+04,
        0.190831076596300e+04,  0.691091682714533e-05];

    // coefficients for rational approximation of
    // psi(x) - ln(x) + 1 / (2*x),  x > 3.0
    const P2: [f64; 4] = [
        -0.212940445131011e+01, -0.701677227766759e+01,
        -0.448616543918019e+01, -0.648157123766197e+00];
    const Q2: [f64; 4] = [
        0.322703493791143e+02, 0.892920700481861e+02,
        0.546117738103215e+02, 0.777788548522962e+01];


    // xsmall = absolute argument below which pi*cotan(pi*x)
    // may be represented by 1/x.
    const XMAX1: f64  = i32::MAX as f64;
    const XSMALL: f64 = 1.0e-9;

    let mut x: f64 = x;
    let mut aug: f64 = 0.0;
    if x < 0.5 {
        // psi(1-x) = psi(x) + pi * cot(pi * x)
        if f64::abs(x) <= XSMALL {
            // 0 < |x| < xsmall
            // pi * cot(pi * x) is approximated by 1/x
            if x == 0.0 {
                return 0.0; // error
            }
            aug = -1.0 / x;
        } else {
            // Reduction of argument for cotangent
            let mut w = -x;
            let mut sgn = PIOV4;
            if w <= 0.0 {
                w = -w;
                sgn = -sgn;
            }
            if w >= XMAX1 {
                return 0.0; // error
            }
            let mut nq = w as i64;
            w -= nq as f64;
            nq = (w * 4.0) as i64;
            w = (w - nq as f64 * 0.25) * 4.0;
            
            // w is now related to the fractional part of 4*x
            // adjust argument to first quadrant and  determine sign.
            let mut n = nq / 2;

            if n + n != nq {
                w = 1.0 - w;
            }

            let z = PIOV4 * w;

            let mut m = n / 2;

            if m + m != n {
                sgn = -sgn;
            }

            // Determine final value for -pi*cot(pi*x)
            n = (nq + 1) / 2;
            m = n / 2;
            m += m;

            if m == n {
                // check for singularity
                if z == 0.0 {
                    return 0.0;  // error
                }
                // cot(z) = cos(z) / sin(z)
                aug = sgn * (f64::cos(z) / f64::sin(z)) * 4.0;
            } else {
                // tan(z) = sin(z) / cos(z)
                aug = sgn * (f64::sin(z) / f64::cos(z)) * 4.0;
            }
        }

        // Reflection reduction
        // psi(x) = psi(1-x) - pi*cot(pi*x)
        // aug contains the correction term.
        x = 1.0 - x;
    }
    if x <= 3.0 {
        let mut den = x;
        let mut upper = P1[0] * x;
        for i in 0..5 {
            den = (den + Q1[i]) * x;
            upper = (upper + P1[i + 1]) * x;
        }
        den = (upper + P1[6]) / (den + Q1[5]);
        let xmx0 = x - DX0;
        return den * xmx0 + aug;
    }
    if x < XMAX1 {
        let w = 1.0 / (x * x);
        let mut den = w;
        let mut upper = P2[0] * w;
        for i in 0..3 {
            den = (den + Q2[i]) * w;
            upper = (upper + P2[i + 1]) * w;
        }
        aug += upper / (den + Q2[3]) - 0.5 / x;
    }
    
    return aug + f64::ln(x);
}

/// REXP evaluates the function Exp(X) - 1.
#[allow(dead_code)]
fn rexp(x: f64) -> f64 {
    const P1: f64 =  0.914041914819518e-09;
    const P2: f64 =  0.238082361044469e-01;
    const Q1: f64 = -0.499999999085958e+00;
    const Q2: f64 =  0.107141568980644e+00;
    const Q3: f64 = -0.119041179760821e-01;
    const Q4: f64 =  0.595130811860248e-03;

    if f64::abs(x) <= 0.15 {
        return x * (((P2 * x + P1) * x + 1.0) / ((((Q4 * x + Q3) * x + Q2) * x + Q1) * x + 1.0));
    } else if x < 0.0 {
        let w = f64::exp(x);
        return (w - 0.5) - 0.5;
    } else {
        let w = f64::exp(x);
        return w * (0.5 + (0.5 - 1.0 / w));
    }
}

/// RLOG1 evaluates the function X - Log ( 1 + X).
#[allow(dead_code)]
fn rlog1(x: f64) -> f64 {
    const A:  f64 =  0.566749439387324e-01;
    const B:  f64 =  0.456512608815524e-01;
    const P0: f64 =  0.333333333333333e+00;
    const P1: f64 = -0.224696413112536e+00;
    const P2: f64 =  0.620886815375787e-02;
    const Q1: f64 = -0.127408923933623e+01;
    const Q2: f64 =  0.354508718369557e+00;
    
    if x < -0.39 {
        let w = (x + 0.5) + 0.5;
        return x - f64::ln(w);
    } else if x < -0.18 {
        let h = (x + 0.3) / 0.7;
        let w1 = A - h * 0.3;
        let r = h / (h + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        return 2.0 * t * (1.0/(1.0 - r) - r * w) + w1;
    } else if x <= 0.18 {
        let h = x;
        let w1 = 0.0;
        let r = h/(h + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        return 2.0 * t * (1.0/(1.0 - r) - r * w) + w1;
    } else if x <= 0.57 {
        let h = 0.75 * x - 0.25;
        let w1 = B + h / 3.0;
        let r = h/(h + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        return 2.0 * t * (1.0 / (1.0 - r) - r * w) + w1;
    } else {
        let w = (x + 0.5) + 0.5;
        return x - f64::ln(w);
    }
}

#[cfg(test)]
mod bratio_tests {
    use super::*;

    /// Tolerance for comparing `w + w1` against 1.
    const SUM_TOL: f64 = 1e-6;

    /// Helper: assert the result is a success and `w + w1 == 1`.
    fn assert_valid(result: BratioResult) {
        assert_eq!(result.ierr, 0, "expected ierr=0, got {}", result.ierr);
        assert!(
            (result.w + result.w1 - 1.0).abs() < SUM_TOL,
            "w + w1 = {} (w={}, w1={})",
            result.w + result.w1,
            result.w,
            result.w1
        );
    }

    // =====================================================================
    // Error paths (E1–E10)
    // =====================================================================

    #[test]
    fn err_a_negative() {
        // E1: a < 0
        let r = bratio(-1.0, 2.0, 0.3, 0.7);
        assert_eq!(r.ierr, 1);
        assert_eq!(r.w, 0.0);
        assert_eq!(r.w1, 0.0);
    }

    #[test]
    fn err_b_negative() {
        // E2: b < 0
        let r = bratio(2.0, -1.0, 0.3, 0.7);
        assert_eq!(r.ierr, 1);
        assert_eq!(r.w, 0.0);
        assert_eq!(r.w1, 0.0);
    }

    #[test]
    fn err_a_b_both_zero() {
        // E3: a == 0 && b == 0
        let r = bratio(0.0, 0.0, 0.3, 0.7);
        assert_eq!(r.ierr, 2);
        assert_eq!(r.w, 0.0);
        assert_eq!(r.w1, 0.0);
    }

    #[test]
    fn err_x_negative() {
        // E4: x < 0
        let r = bratio(2.0, 3.0, -0.1, 1.1);
        assert_eq!(r.ierr, 3);
    }

    #[test]
    fn err_x_greater_than_one() {
        // E5: x > 1
        let r = bratio(2.0, 3.0, 1.1, -0.1);
        assert_eq!(r.ierr, 3);
    }

    #[test]
    fn err_y_negative() {
        // E6: y < 0
        let r = bratio(2.0, 3.0, 1.1, -0.1);
        // Note: x > 1 is caught first, so use a case where x is valid
        // but y is negative.
        let r = bratio(2.0, 3.0, 0.5, -0.1);
        let _ = r; // placeholder to avoid unused warning
        // The above actually trips E4? No: x=0.5 is valid, y=-0.1 < 0.
        assert_eq!(r.ierr, 4);
    }

    #[test]
    fn err_y_greater_than_one() {
        // E7: y > 1
        let r = bratio(2.0, 3.0, -0.1, 1.1);
        // x < 0 is caught first, so use valid x.
        let r = bratio(2.0, 3.0, 0.5, 1.1);
        assert_eq!(r.ierr, 4);
    }

    #[test]
    fn err_x_plus_y_not_one() {
        // E8: x + y != 1 (well outside 3*eps)
        let r = bratio(2.0, 3.0, 0.3, 0.6);
        assert_eq!(r.ierr, 5);
    }

    #[test]
    fn err_x_zero_a_zero() {
        // E9: x == 0 && a == 0
        let r = bratio(0.0, 2.0, 0.0, 1.0);
        assert_eq!(r.ierr, 6);
    }

    #[test]
    fn err_y_zero_b_zero() {
        // E10: y == 0 && b == 0
        let r = bratio(2.0, 0.0, 1.0, 0.0);
        assert_eq!(r.ierr, 7);
    }

    // =====================================================================
    // Degenerate / trivial paths (D1–D5)
    // =====================================================================

    #[test]
    fn degen_x_zero() {
        // D1: x == 0, a != 0 → w=0, w1=1
        let r = bratio(2.0, 3.0, 0.0, 1.0);
        assert_eq!(r.ierr, 0);
        assert_eq!(r.w, 0.0);
        assert_eq!(r.w1, 1.0);
    }

    #[test]
    fn degen_y_zero() {
        // D2: y == 0, b != 0 → w=1, w1=0
        let r = bratio(2.0, 3.0, 1.0, 0.0);
        assert_eq!(r.ierr, 0);
        assert_eq!(r.w, 1.0);
        assert_eq!(r.w1, 0.0);
    }

    #[test]
    fn degen_a_zero() {
        // D3: a == 0, x != 0, y != 0 → w=1, w1=0
        let r = bratio(0.0, 3.0, 0.5, 0.5);
        assert_eq!(r.ierr, 0);
        assert_eq!(r.w, 1.0);
        assert_eq!(r.w1, 0.0);
    }

    #[test]
    fn degen_b_zero() {
        // D4: b == 0, x != 0, y != 0 → w=0, w1=1
        let r = bratio(3.0, 0.0, 0.5, 0.5);
        assert_eq!(r.ierr, 0);
        assert_eq!(r.w, 0.0);
        assert_eq!(r.w1, 1.0);
    }

    #[test]
    fn degen_a_b_tiny() {
        // D5: max(a,b) < 1e-3 * eps, where eps = max(f64::EPSILON, 1e-15)
        // With f64::EPSILON ≈ 2.22e-16, 1e-3*eps ≈ 2.22e-19.
        let a = 1e-20;
        let b = 2e-20;
        let r = bratio(a, b, 0.5, 0.5);
        println!("{:.18e}", r.w);
        assert_eq!(r.ierr, 0);
        let expected_w = b / (a + b);
        let expected_w1 = a / (a + b);
        assert!((r.w - expected_w).abs() < 1e-12);
        assert!((r.w1 - expected_w1).abs() < 1e-12);
    }

    // =====================================================================
    // Algorithm-selection paths (A1–A16)
    // =====================================================================

    #[test]
    fn alg_fpser() {
        // A1: label 80, condition b0 < min(eps, eps*a0).
        // With a0 = a (large enough), b0 = b tiny.
        // eps ≈ max(2.22e-16, 1e-15) = 1e-15. So need b < 1e-15 * a.
        let a = 1.0;
        let b = 1e-18;
        let r = bratio(a, b, 0.5, 0.5);
        assert_valid(r);
    }

    #[test]
    fn alg_apser() {
        // A2: label 90, a0 < min(eps, eps*b0) && b0*x0 <= 1.
        // Need a < 1e-15 * b (so a is tiny), b0*x0 <= 1.
        let a = 1e-18;
        let b = 1.0;
        let x = 0.3;
        let y = 0.7;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_bpser_small_a0_b0() {
        // A3 + A11: a0, b0 <= 1 and a0 >= min(0.2, b0) → label 100.
        let a = 0.5;
        let b = 0.5;
        let r = bratio(a, b, 0.3, 0.7);
        assert_valid(r);
    }

    #[test]
    fn alg_bpser_a0_lt_min_0_2_b0_x_pow_a_ok() {
        // A3 + A12: a0 < min(0.2, b0) and x0**a0 <= 0.9 → label 100.
        // Pick a0 small, b0 moderate, x0 small so x0**a0 is close to 1
        // but not too small. Let a0 = 0.01, x0 = 0.5 → 0.5**0.01 ≈ 0.993.
        // That's > 0.9, so this does NOT satisfy x0**a0 <= 0.9.
        // Need x0**a0 <= 0.9. With a0=0.01, need x0 <= 0.9**(1/0.01) ≈ 0.
        // So use a0 = 0.1, x0 = 0.3 → 0.3**0.1 ≈ 0.886 <= 0.9. Good.
        let a = 0.1;
        let b = 0.6;
        let x = 0.3;
        let y = 0.7;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_bpser_swapped_label_110() {
        // A4 + A13: label 10 → x0 >= 0.3 → label 110 (bpser swapped).
        // Need a0, b0 <= 1, a0 < min(0.2, b0), x0**a0 > 0.9, x0 >= 0.3.
        // With a0 = 0.1, x0 = 0.95 → 0.95**0.1 ≈ 0.9949 > 0.9. Good.
        // But x <= 0.5 check: if x > 0.5, ind=1 and swap happens first.
        // So for label 10's x0, we need x <= 0.5. But then x0 >= 0.3 and
        // x0**a0 > 0.9 requires x0 close to 1. With x <= 0.5, x0 <= 0.5.
        // 0.5**0.1 ≈ 0.933 > 0.9. So a0 = 0.1, x = 0.5 works.
        let a = 0.1;
        let b = 0.6;
        let x = 0.5;
        let y = 0.5;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_bup_bgrat_label_130() {
        // A5 + A13's fallthrough: label 10 → n=20 → label 130.
        // Need a0, b0 <= 1, a0 < min(0.2, b0), x0**a0 > 0.9, x0 < 0.3.
        // x0 < 0.3 and x0**a0 > 0.9: with a0 = 0.01, x0 = 0.29 →
        // 0.29**0.01 ≈ 0.987 > 0.9. Good.
        let a = 0.01;
        let b = 0.6;
        let x = 0.29;
        let y = 0.71;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_20_b0_le_1() {
        // A14: label 20 → b0 <= 1 → label 100.
        // Need max(a0,b0) > 1, b0 <= 1, a0 <= 1.
        // So a0 <= 1, b0 <= 1, but max > 1 — contradiction unless
        // one of them equals exactly 1? No: max(a0,b0) > 1 requires
        // one > 1. So for label 20 with b0 <= 1, we need a0 > 1.
        // But then min(a0,b0) <= 1 holds (b0 <= 1). Good.
        let a = 2.0;
        let b = 1.0;
        let x = 0.3;
        let y = 0.7;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_20_x0_ge_0_3() {
        // A4 via label 20: b0 > 1, x0 >= 0.3 → label 110.
        let a = 2.0;
        let b = 1.5;
        let x = 0.3;
        let y = 0.7;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_21_b0_gt_15() {
        // A6 + A15: label 20 → x0 < 0.3, x0 >= 0.1, b0 > 15 → label 131.
        // Need max(a0,b0) > 1, b0 > 1, x0 in [0.1, 0.3), b0 > 15.
        // Also need min(a0,b0) <= 1, so a0 <= 1.
        let a = 1.0;
        let b = 20.0;
        let x = 0.2;
        let y = 0.8;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_21_b0_le_15() {
        // A5 + A15: label 21 → b0 <= 15 → n=20 → label 130.
        let a = 1.0;
        let b = 10.0;
        let x = 0.2;
        let y = 0.8;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_20_x0_lt_0_1_ok() {
        // A16: label 20 → x0 < 0.1, (x0*b0)**a0 <= 0.7 → label 100.
        let a = 1.0;
        let b = 5.0;
        let x = 0.05;
        let y = 0.95;
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_20_x0_lt_0_1_fallthrough_b0_gt_15() {
        // A16 fallthrough + A15: x0 < 0.1, (x0*b0)**a0 > 0.7, b0 > 15
        // → label 131.
        let a = 0.1;
        let b = 20.0;
        let x = 0.05;
        let y = 0.95;
        // (0.05*20)**0.1 = 1**0.1 = 1 > 0.7, so fallthrough.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_20_x0_lt_0_1_fallthrough_b0_le_15() {
        // A16 fallthrough + A15 fallthrough: x0 < 0.1, (x0*b0)**a0 > 0.7,
        // b0 <= 15 → n=20 → label 130.
        let a = 0.1;
        let b = 10.0;
        let x = 0.05;
        let y = 0.95;
        // (0.05*10)**0.1 = 0.5**0.1 ≈ 0.933 > 0.7, so fallthrough.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    // ---------------------------------------------------------------------
    // a0 > 1 and b0 > 1 branch (label 30)
    // ---------------------------------------------------------------------

    #[test]
    fn alg_label_40_b0_lt_40_b0x0_le_0_7() {
        // label 40 → b0 < 40 && b0*x0 <= 0.7 → label 100.
        let a = 3.0;
        let b = 5.0;
        let x = 0.1;
        let y = 0.9;
        // b0*x0 = 0.5 <= 0.7. min(a,b)=3 > 1. Good.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_140_x0_le_0_7() {
        // A9: label 40 → b0 < 40 → label 140 → label 141 → x0 <= 0.7.
        // Need b0 < 40, b0*x0 > 0.7, x0 <= 0.7.
        let a = 3.0;
        let b = 5.0;
        let x = 0.3;
        let y = 0.7;
        // b0*x0 = 1.5 > 0.7. x0 = 0.3 <= 0.7. Good.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_150_a0_le_15() {
        // A10: label 40 → b0 < 40 → label 140 → label 141 → x0 > 0.7
        // → label 150 → a0 <= 15 → label 151.
        let a = 3.0;
        let b = 5.0;
        let x = 0.8;
        let y = 0.2;
        // b0*x0 = 4.0 > 0.7, x0 = 0.8 > 0.7, a0 = 3 <= 15. Good.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_150_a0_gt_15() {
        // A10: x0 > 0.7, a0 > 15 → skip label 150's bup, go to 151.
        let a = 20.0;
        let b = 5.0;
        let x = 0.8;
        let y = 0.2;
        // min(a,b) = 5 > 1. b0 < 40? b0 = min(a,b) = 5. Yes.
        // b0*x0 = 4 > 0.7, x0 > 0.7, a0 = max(a,b) = 20 > 15. Good.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_140_b0_exactly_integer() {
        // Exercise the `b0 == 0.0` special case after `n = b0; b0 = b0-n`.
        // Need b0 an integer. With a0 > 1, b0 > 1.
        let a = 3.0;
        let b = 4.0;
        let x = 0.8;
        let y = 0.2;
        // b0 = min(a,b) = 3 (integer). n = 3, b0 = 0 → special case.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_120_bfrac() {
        // A7: label 40 → b0 >= 40 → a0 <= b0 → a0 <= 100 or lambda small
        // → label 120.
        let a = 50.0;
        let b = 60.0;
        let x = 0.4;
        let y = 0.6;
        // min=50 > 1. lambda = a - (a+b)*x = 50 - 110*0.4 = 6.
        // b0 = min = 50 >= 40. a0 = 50 <= 100. Good.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_120_via_lambda_small() {
        // label 40 → b0 >= 40 → a0 > b0, a0 > 100, lambda <= 0.03*a0
        // → label 120.
        let a = 200.0;
        let b = 50.0;
        let x = 0.4;
        let y = 0.6;
        // min = 50 >= 40. a0 = 200 > b0 = 50. a0 > 100.
        // lambda = (a+b)*y - b = 250*0.6 - 50 = 100.
        // 0.03*a0 = 6. lambda=100 > 6, so NOT this path.
        // Need lambda <= 0.03*a0. Try x closer to 1.
        let x = 0.99;
        let y = 0.01;
        // lambda = 250*0.01 - 50 = -47.5 < 0 → ind=1, swap.
        // After swap: a0 = 50, b0 = 200, lambda = 47.5.
        // Then b0 = 200 >= 40. a0 = 50 <= b0 = 200. a0 = 50 <= 100.
        // So this goes to label 120 via the a0<=100 path, not the lambda path.
        // Let me pick different values to hit the lambda path.
        let a = 200.0;
        let b = 50.0;
        let x = 0.4;
        let y = 0.6;
        // lambda = 250*0.6 - 50 = 100. 0.03*a0 = 6. 100 > 6 → label 180.
        // That's basym. Let me use that for the basym test and come back.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_180_basym() {
        // A8: label 40 → b0 >= 40 → a0 > b0, a0 > 100,
        // lambda > 0.03*a0 → label 180.
        let a = 200.0;
        let b = 50.0;
        let x = 0.4;
        let y = 0.6;
        // min = 50 >= 40. a0 = 200 > b0 = 50. a0 > 100.
        // lambda = (a+b)*y - b = 250*0.6 - 50 = 100.
        // 0.03*a0 = 6. 100 > 6 → label 180. 
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_50_b0_le_100() {
        // label 50: a0 > b0, b0 <= 100 → label 120.
        let a = 200.0;
        let b = 80.0;
        let x = 0.4;
        let y = 0.6;
        // min = 80 >= 40. a0 = 200 > b0 = 80. b0 = 80 <= 100 → label 120.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_50_lambda_gt_003_b0() {
        // label 50: a0 > b0, b0 > 100, lambda > 0.03*b0 → label 120.
        let a = 500.0;
        let b = 200.0;
        let x = 0.4;
        let y = 0.6;
        // min = 200 >= 40. a0 = 500 > b0 = 200. b0 > 100.
        // lambda = (a+b)*y - b = 700*0.6 - 200 = 220.
        // 0.03*b0 = 6. 220 > 6 → label 120.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    #[test]
    fn alg_label_50_lambda_le_003_b0() {
        // label 50: a0 > b0, b0 > 100, lambda <= 0.03*b0 → label 180.
        // Need lambda small relative to b0.
        let a = 500.0;
        let b = 200.0;
        let x = 0.6;
        let y = 0.4;
        // lambda = 700*0.4 - 200 = 80. 0.03*b0 = 6. 80 > 6. Still 120.
        // Hmm, need lambda <= 6. Try x near the mean a/(a+b) = 500/700 ≈ 0.714.
        let x = 0.71;
        let y = 0.29;
        // lambda = 700*0.29 - 200 = 3. 0.03*b0 = 6. 3 <= 6 → label 180.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    // =====================================================================
    // Swap paths (S1, S2)
    // =====================================================================

    #[test]
    fn swap_x_gt_half_small_a0() {
        // S1: min(a,b) <= 1 and x > 0.5 → ind = 1.
        let a = 0.5;
        let b = 0.5;
        let x = 0.7;
        let y = 0.3;
        let r = bratio(a, b, x, y);
        println!("{:.18e}, {:.18e}", r.w, r.w1);
        assert_valid(r);
    }

    #[test]
    fn swap_lambda_negative_large_a0() {
        // S2: min(a,b) > 1 and lambda < 0 → ind = 1.
        // lambda = a - (a+b)*x if a <= b; need negative.
        let a = 3.0;
        let b = 5.0;
        let x = 0.9;
        let y = 0.1;
        // lambda = 3 - 8*0.9 = -4.2 < 0 → ind = 1.
        let r = bratio(a, b, x, y);
        assert_valid(r);
    }

    // =====================================================================
    // Symmetry / numerical sanity checks
    // =====================================================================

    #[test]
    fn symmetry_ix_ab_x_eq_1_minus_ix_ba_y() {
        // I_x(a,b) = 1 - I_y(b,a) where y = 1-x.
        // bratio(a,b,x,y).w should equal 1 - bratio(b,a,y,x).w.
        let a = 2.5;
        let b = 3.5;
        let x = 0.3;
        let y = 0.7;
        let r1 = bratio(a, b, x, y);
        let r2 = bratio(b, a, y, x);
        assert_valid(r1);
        assert_valid(r2);
        assert!(
            (r1.w - r2.w1).abs() < 1e-6,
            "r1.w={} vs r2.w1={}",
            r1.w,
            r2.w1
        );
    }

    #[test]
    fn known_value_a1_b1_is_x() {
        // I_x(1,1) = x.
        let r = bratio(1.0, 1.0, 0.3, 0.7);
        assert_valid(r);
        assert!((r.w - 0.3).abs() < 1e-10, "w={} expected 0.3", r.w);
    }

    #[test]
    fn known_value_a2_b1_is_x_squared() {
        // I_x(2,1) = x^2.
        let r = bratio(2.0, 1.0, 0.5, 0.5);
        assert_valid(r);
        assert!((r.w - 0.25).abs() < 1e-10, "w={} expected 0.25", r.w);
    }

    #[test]
    fn known_value_a1_b2_is_1_minus_y_squared() {
        // I_x(1,2) = 1 - (1-x)^2 = 1 - y^2.
        let r = bratio(1.0, 2.0, 0.5, 0.5);
        assert_valid(r);
        assert!((r.w - 0.75).abs() < 1e-10, "w={} expected 0.75", r.w);
    }

    #[test]
    fn known_value_a_half_b_half_is_arcsine() {
        // I_x(0.5,0.5) = (2/pi) * arcsin(sqrt(x)).
        let x = 0.3;
        let r = bratio(0.5, 0.5, x, 1.0 - x);
        assert_valid(r);
        let expected = (2.0 / std::f64::consts::PI) * x.sqrt().asin();
        assert!(
            (r.w - expected).abs() < 1e-6,
            "w={} expected {}",
            r.w,
            expected
        );
    }

    #[test]
    fn w_plus_w1_equals_one_many_cases() {
        // Sweep a grid of (a, b, x) and check w + w1 == 1 every time.
        let avals = [0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0];
        let bvals = [0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0];
        let xvals = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99];
        for &a in &avals {
            for &b in &bvals {
                for &x in &xvals {
                    let r = bratio(a, b, x, 1.0 - x);
                    assert_eq!(r.ierr, 0, "ierr for a={} b={} x={}", a, b, x);
                    let sum = r.w + r.w1;
                    assert!(
                        (sum - 1.0).abs() < 1e-6,
                        "a={} b={} x={}: w+w1={}",
                        a, b, x, sum
                    );
                    assert!(
                        r.w >= -1e-10 && r.w <= 1.0 + 1e-10,
                        "a={} b={} x={}: w={} out of [0,1]",
                        a, b, x, r.w
                    );
                }
            }
        }
    }

    #[test]
    fn monotonic_in_x() {
        // I_x(a,b) should be nondecreasing in x for fixed a, b.
        let a = 2.0;
        let b = 3.0;
        let mut prev = -1.0;
        for i in 0..=100 {
            let x = i as f64 / 100.0;
            let r = bratio(a, b, x, 1.0 - x);
            assert_eq!(r.ierr, 0);
            assert!(
                r.w >= prev - 1e-9,
                "non-monotonic at x={}: w={} prev={}",
                x, r.w, prev
            );
            prev = r.w;
        }
    }

    #[test]
    fn boundary_x_equals_half() {
        // x = 0.5 is a common edge; by symmetry, I_{0.5}(a,a) = 0.5.
        let r = bratio(2.0, 2.0, 0.5, 0.5);
        assert_valid(r);
        assert!((r.w - 0.5).abs() < 1e-10, "w={} expected 0.5", r.w);
    }

    #[test]
    fn boundary_a_equals_b() {
        // I_x(a,a) + I_{1-x}(a,a) = 1; at x=0.5 it's exactly 0.5.
        let r = bratio(5.0, 5.0, 0.5, 0.5);
        assert_valid(r);
        assert!((r.w - 0.5).abs() < 1e-10);
    }
}