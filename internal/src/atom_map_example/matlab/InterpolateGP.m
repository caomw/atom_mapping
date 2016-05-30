%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to evaluate the GP regression at a particular 3D point.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sdf = InterpolateGP(x, y, z, kdtree, data)
    % kNN search.
    neighbors_idx = knnsearch(kdtree, [x, y, z], 'k', 30);
    num_neighbors = length(neighbors_idx);
    
    % Set gamma parameter.
    GAMMA = 1.0;
    
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
    
    % Compute cross covariance.
    K12 = zeros(num_neighbors, 1);
    for ii = 1:num_neighbors
       K12(ii) = ...
          CovarianceKernel(data(neighbors_idx(ii), 1:3), [x, y, z], GAMMA);       
    end
    
    % Compute expected sdf.
    sdf = K12' * (K11 \ data(neighbors_idx, 4));
end