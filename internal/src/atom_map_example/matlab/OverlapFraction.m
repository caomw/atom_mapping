%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute overlap fraction between atoms separated by distance
% 'd' and of radius 'r'.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = OverlapFraction(d, r)
w = 1 - (0.25 / (r*r*r)) * d * (3.0 * r*r - 0.25 * d*d);
end