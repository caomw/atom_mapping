%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to evaluate the GP regression at a particular 3D point.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sdf, var] = InterpolateGP(x, y, z, kdtree, data)
    % kNN search.
    neighbors_idx = knnsearch(kdtree, [x, y, z], 'k', 30);
    num_neighbors = length(neighbors_idx);
    
    % Set gamma and noise variance parameters.
    GAMMA = 0.7;
    NOISE = 0.5;
    
    % Compute training covariance.
    K11 = zeros(num_neighbors, num_neighbors);
    for ii = 1:num_neighbors
       for jj = 1:ii
           var = CovarianceKernel(data(neighbors_idx(ii), 1:3), ...
                                  data(neighbors_idx(jj), 1:3), GAMMA);
           K11(ii, jj) = var;
           K11(jj, ii) = var;
       end
    end
    
    K11 = K11 + NOISE * eye(size(K11));
    
    % Compute cross covariance.
    K12 = zeros(num_neighbors, 1);
    for ii = 1:num_neighbors
       K12(ii) = ...
          CovarianceKernel(data(neighbors_idx(ii), 1:3), [x, y, z], GAMMA);       
    end
    
    % Compute expected sdf.
    sdf = K12' * (K11 \ data(neighbors_idx, 4));
    
    % Compute variance.
    var = 1.0 - K12' * (K11 \ K12);
end