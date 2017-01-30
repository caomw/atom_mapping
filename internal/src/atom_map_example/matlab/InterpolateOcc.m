%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to interpolate log odds of occupancy at a particular point,
% using multi-linear interpolation of log odds.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function occ = InterpolateOcc(x, y, z, kdtree, data, r)

    % Radius search.
    [neighbors_idx, distances] = rangesearch(kdtree, [x, y, z], 2.0*r);
    
    if numel(neighbors_idx{1,1}) == 0
        % Default to occupied.
        occ = 10.0;
    else        
        % Compute weighted average.
        total_weight = 0.0;
        weighted_sum = 0.0;
        for ii = 1:numel(neighbors_idx{1,1})
            idx = neighbors_idx{1,1}(ii);
            d = distances{1,1}(ii);
            weight = OverlapFraction(d, r);

            total_weight = total_weight + weight;
            weighted_sum = weighted_sum + weight * data(idx, 4);
        end

        occ = weighted_sum / total_weight;    
    end
end