function [a, b] = gamShapeScale(m, sd)

% Return the shape and scale parameters for a gamma dsitribution with mean
% m and std. dev. sd

a = m.^2./sd.^2;
b = m./a;
