use crate::erfc;
use std::f64::consts::SQRT_2;

/*
function [varargout] = normcdf(x,varargin)
%NORMCDF Normal cumulative distribution function (cdf).
%   P = NORMCDF(X,MU,SIGMA) returns the cdf of the normal distribution with
%   mean MU and standard deviation SIGMA, evaluated at the values in X.
%   The size of P is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [P,PLO,PUP] = NORMCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  PLO and PUP are arrays of the same size as P containing the lower
%   and upper confidence bounds.
%
%   [...] = NORMCDF(...,'upper') computes the upper tail probability of the 
%   normal distribution. This can be used to compute a right-tailed p-value. 
%   To compute a two-tailed p-value, use 2*NORMCDF(-ABS(X),MU,SIGMA).
%
%   See also ERF, ERFC, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2024 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<1
    error(message('stats:normcdf:TooFewInputsX'));
end

if nargin>1 && strcmpi(varargin{end},'upper')
    % Compute upper tail and remove 'upper' flag
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;  
end

[varargout{1:max(1,nargout)}] = localnormcdf(uflag,x,varargin{:});
end

function [p,plo,pup] = localnormcdf(uflag,x,mu,sigma,pcov,alpha)

if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% Inputs need to match in size. This is guaranteed for the common case of
% scalar mu and sigma.
if ~isscalar(mu) || ~isscalar(sigma)
    [errorcode, x,mu,sigma] = distchck(3,x,mu,sigma);
    if errorcode > 0
        error(message('stats:normcdf:InputSizeMismatch'));
    end
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<5
      error(message('stats:normcdf:TooFewInputsCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:normcdf:BadCovarianceSize'));
   end
   if nargin<6
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:normcdf:BadAlpha'));
   end
end

z = (x-mu) ./ sigma;
if uflag==true
    z = -z;
end

% Prepare output
p = NaN(size(z),"like",z);
if nargout>=2
    plo = NaN(size(z),"like",z);
    pup = NaN(size(z),"like",z);
end

% Set edge case sigma=0
if any(sigma==0, "all")
    isSigmaZero = sigma==0;
    xLtMu = isSigmaZero & x<mu;
    xGtMu = isSigmaZero & x>=mu;

    if uflag==true
        p(xLtMu) = 1;
        p(xGtMu) = 0;
        if nargout>=2
            plo(xLtMu) = 1;
            plo(xGtMu) = 0;
            pup(xLtMu) = 1;
            pup(xGtMu) = 0;
        end
    else
        p(xLtMu) = 0;
        p(xGtMu) = 1;
        if nargout>=2
            plo(xLtMu) = 0;
            plo(xGtMu) = 1;
            pup(xLtMu) = 0;
            pup(xGtMu) = 1;
        end
    end
end

% Normal cases
if isscalar(sigma)
    if sigma > 0
        p = iErf(z);
        % Compute confidence bounds if requested.
        if nargout > 1
            [plo(:), pup(:)] = iComputeConfidenceBounds(z, sigma, alpha, pcov);
        end
    end
    return
end

todo = sigma>0;
sigma = sigma(todo);
z = z(todo);
p(todo) = iErf(z);
% Compute confidence bounds if requested.
if nargout > 1
    [plo(todo), pup(todo)] = iComputeConfidenceBounds(z, sigma, alpha, pcov);
end
end

function p = iErf(z)
% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p = 0.5 * erfc(-z ./ sqrt(2));
end

function [plo, pup] = iComputeConfidenceBounds(z,sigma,alpha,pcov)
z = z(:);
zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma(:).^2);
if any(zvar<0)
    error(message('stats:normcdf:BadCovarianceSymPos'));
end
normz = -norminv(alpha/2);
halfwidth = normz * sqrt(zvar);
zlo = z - halfwidth;
zup = z + halfwidth;

plo = iErf(zlo);
pup = iErf(zup);
end

 */

/// Normal cumulative distribution function (CDF).
///
/// Given a value `x`, a mean `mu`, and a standard deviation `sigma`,
/// this function returns the probability that a normal random variable
/// is less than or equal to `x`.
///
/// # Mathematical Definition
/// For a normal distribution with mean <math><mi>μ</mi></math> and standard deviation <math><mi>σ</mi></math>:
/// - Lower tail (`upper = false`): <math><mi>Φ</mi><mo>(</mo><mi>x</mi><mo>;</mo><mi>μ</mi><mo>,</mo><mi>σ</mi><mo>)</mo><mo>=</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mi>erfc</mi><mo>(</mo><mo>-</mo><mfrac><mrow><mi>x</mi><mo>-</mo><mi>μ</mi></mrow><mrow><mi>σ</mi><msqrt><mn>2</mn></msqrt></mrow></mfrac><mo>)</mo></math>
/// - Upper tail (`upper = true`): <math><mn>1</mn><mo>-</mo><mi>Φ</mi><mo>(</mo><mi>x</mi><mo>;</mo><mi>μ</mi><mo>,</mo><mi>σ</mi><mo>)</mo><mo>=</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><mi>erfc</mi><mo>(</mo><mfrac><mrow><mi>x</mi><mo>-</mo><mi>μ</mi></mrow><mrow><mi>σ</mi><msqrt><mn>2</mn></msqrt></mrow></mfrac><mo>)</mo></math>
///
/// # Domain
/// - `sigma >= 0`
/// - Returns `NaN` if `sigma` is negative or any input is `NaN`.
/// - If `sigma` is 0, the distribution is treated as a Dirac delta concentrated at `mu`.
///   - Returns 0 (lower) / 1 (upper) if `x < mu`.
///   - Returns 1 (lower) / 0 (upper) if `x >= mu`.
///
/// # Examples
/// ```
/// use abax::normcdf;
///
/// // Standard normal median is 0.5
/// assert!((normcdf(0.0, 0.0, 1.0, false) - 0.5).abs() < 1e-15);
/// // One sigma upper bound
/// assert!((normcdf(1.0, 0.0, 1.0, false) - 0.8413447460685429).abs() < 1e-15);
/// // Upper tail one sigma
/// assert!((normcdf(1.0, 0.0, 1.0, true) - 0.15865525393145705).abs() < 1e-15);
/// ```
pub fn normcdf(x: f64, mu: f64, sigma: f64, upper: bool) -> f64 {
    if x.is_nan() || mu.is_nan() || sigma.is_nan() || sigma < 0.0 {
        return f64::NAN;
    }
 
    if sigma == 0.0 {
        return if upper {
            if x < mu { 1.0 } else { 0.0 }
        } else {
            if x < mu { 0.0 } else { 1.0 }
        };
    }
 
    let z = (x - mu) / sigma;
    if upper {
        0.5 * erfc(z / SQRT_2)
    } else {
        0.5 * erfc(-z / SQRT_2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normcdf_standard() {
        let tol = 1e-14;
        assert!((normcdf(0.0, 0.0, 1.0, false) - 0.5).abs() < tol);
        assert!((normcdf(1.959963984540054, 0.0, 1.0, false) - 0.975).abs() < tol);
        assert!((normcdf(-1.959963984540054, 0.0, 1.0, false) - 0.025).abs() < tol);
    }
 
    #[test]
    fn test_normcdf_upper() {
        let tol = 1e-14;
        assert!((normcdf(0.0, 0.0, 1.0, true) - 0.5).abs() < tol);
        assert!((normcdf(1.959963984540054, 0.0, 1.0, true) - 0.025).abs() < tol);
        assert!((normcdf(-1.959963984540054, 0.0, 1.0, true) - 0.975).abs() < tol);
    }

    #[test]
    fn test_normcdf_zero_sigma() {
        assert_eq!(normcdf(5.0, 5.0, 0.0, false), 1.0);
        assert_eq!(normcdf(4.99, 5.0, 0.0, false), 0.0);
        assert_eq!(normcdf(4.99, 5.0, 0.0, true), 1.0);
    }

    #[test]
    fn test_normcdf_invalid() {
        assert!(normcdf(0.5, 0.0, -1.0, false).is_nan());
    }
}
