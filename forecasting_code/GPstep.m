function F =  GPstep(xs, Fs, xq, covFn, SD, L, sNoise)

% Function to generate values of a zero-mean Gaussian process at query points xq conditional on the known (noisy) values at sample points xs
%
% USAGE: F =  GPstep(xs, Fs, xq, covFn, SD, L, sNoise) 
%
% INPUTS: xs - row vector of sample points xs
%         Fs - observed value of the GP at the sample points (may be a matrix with nr rows where each row is an independent sample)
%         xq - row vector of nq query points xq
%         covFn - function that returns the covariance for any distance d =
%         |x-x'|, assumed to be isotropic
%         SD - std. dev. for covariance funcion
%         L - length scale L of covariance function
%         sNoise - standard deviation of the IID observation noise (needs to be nonzero to avoid numerical instability)
%
% OUTPUTS: F - nr x nq matrix containing generated values of the GP at the query points - each column of F is a query point and each row of F is the generated data for the corresponding row of Fs

% See p16 of Rasmussen


SMALL = 1e-12;

K_XsX = covFn(xs-xq.', SD, L);
K_XX  = covFn(xs-xs.', SD, L) + sNoise^2*eye(length(xs)); 
m = (K_XsX * (K_XX \ Fs.') ).';
Sigma = covFn(xq-xq.', SD, L) - K_XsX * (K_XX \ K_XsX.');
%Sigma = 0.5*(Sigma+Sigma.');                % force covariance matrix to be symmetric
    
domEV = eigs(Sigma, [], 1);
Sigma = Sigma + max(0, -domEV+SMALL)*eye(size(Sigma));        % force covariance matrix to be positive definite

 
F = mvnrnd(m, Sigma);

