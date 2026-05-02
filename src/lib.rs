//! # abax
//!
//! A lightweight Rust library providing high-precision mathematical constants and special functions.
//! This library includes Bernoulli numbers (<math><msub><mi>B</mi><mrow><mn>2</mn><mi>n</mi></mrow></msub></math>), Riemann Zeta values (<math><mi>ζ</mi><mo>(</mo><mi>s</mi><mo>)</mo></math>), and Stirling series coefficients.

mod consts;
mod digamma;
mod gamma;
mod gammaln;
mod trigamma;

pub use digamma::digamma;
pub use gamma::gamma;
pub use gammaln::gammaln;
pub use trigamma::trigamma;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gammaln_export() {
        assert_eq!(gammaln(1.0), 0.0);
    }

    #[test]
    fn test_gamma_export() {
        assert_eq!(gamma(5.0), 24.0);
    }

    #[test]
    fn test_digamma_export() {
        assert!((digamma(1.0) + 0.5772156649015329).abs() < 1e-14);
    }

    #[test]
    fn test_trigamma_export() {
        assert!((trigamma(1.0) - 1.6449340668482264).abs() < 1e-14);
    }
}
