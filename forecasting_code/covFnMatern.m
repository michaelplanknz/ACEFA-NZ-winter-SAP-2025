function Sigma = covFnMatern(d, s0, L, v)

% Matern covariance function Sigma(d) with parameters s0 (s.d. at d=0), L
% (correlation length scale) and v

z = sqrt(2*v)*abs(d)/L;
Sigma = s0^2 * 2^(1-v)/gamma(v) * z.^v .* besselk(v, z);
Sigma(d == 0) = s0^2;

