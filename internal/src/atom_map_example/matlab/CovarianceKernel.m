%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance kernel.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cov = CovarianceKernel(p1, p2, gamma)
    cov = exp(-gamma * norm(p1 - p2)^2);
end