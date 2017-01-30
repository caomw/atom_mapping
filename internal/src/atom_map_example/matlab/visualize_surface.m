%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization script. Load up an AtomMap, allow user to select a box, and
% compute the implicit surface in that region.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%%
MODE = 0; % 0 = occ, 1 = sdf
RADIUS = 0.3;
SDF_MIN = -1.2;
SDF_MAX = 1.2;

%% Load an AtomMap.
data = csvread('../saved_maps/nsh_300cm_occ.csv');

%% Display a 2D scatterplot of the data projected onto the x-y plane.
%figure;
%scatter(data(:, 1), data(:, 2), ...
%       'MarkerEdgeColor', [0 .5 .5], ...
%       'MarkerFaceColor', [0 .7 .7], ...
%       'LineWidth', 1.5);

%% Let the user select two points (to define a bounding box -- LL then TR).
%fprintf('Please select two points and then press ENTER...\n');
%[pts_x, pts_y] = getpts
pts_x = [77.6728, 81]; %86.3825];
pts_y = [-56.7347, -49.7376];

subset = data(data(:, 1) > pts_x(1) & data(:, 1) < pts_x(2) & ...
              data(:, 2) > pts_y(1) & data(:, 2) < pts_y(2), :);

z_min = min(subset(:, 3))-1;
z_max = max(subset(:, 3))+1;

%% Create a kdtree.
kdtree = KDTreeSearcher(subset(:, 1:3));

%% Allocate a mesh over the specified area.
resolution = 0.25;
[X, Y, Z] = meshgrid(pts_x(1):resolution:pts_x(2), ...
                     pts_y(1):resolution:pts_y(2), ...
                     z_min:resolution:z_max);

size(X)

%% Evaluate GP at each voxel in mesh.
if MODE == 1
    [sdfs, vars] = arrayfun(@(x, y, z) ...
        InterpolateGP(x, y, z, kdtree, subset), X, Y, Z);
else
    [occs, vars] = arrayfun(@(x, y, z) ...
        InterpolateGP(x, y, z, kdtree, subset), X, Y, Z);
end

%% Plot
figure; hold on; set(gca, 'fontsize', 16);

if MODE == 1
    isosurface(X, Y, Z, sdfs, 0, vars); colormap bone;
    freezeColors;
    %scatter3(subset(:, 1), subset(:, 2), subset(:, 3))
    % figure; hold on;

    %pcshow(subset(:, 1:3), 'MarkerSize', 18);
    colormap cool;
    min_sdf = min(subset(:, 4)) - 0.1;
    max_sdf = max(subset(:, 4)) + 0.1;
    [atom_x, atom_y, atom_z] = sphere;
    for ii = 1:size(subset, 1)
        color = ones(size(atom_x)) * (subset(ii, 4) - min_sdf) / (max_sdf - min_sdf);
    %    color(:, 3) = 1.0 - color(:, 1);
        surf(RADIUS * atom_x + subset(ii, 1), ...
             RADIUS * atom_y + subset(ii, 2), ...
             RADIUS * atom_z + subset(ii, 3), ...
             color, 'facealpha', 0.9, 'edgealpha', 0.0);
    end
else
    isosurface(X, Y, Z, occs, 0, vars); colormap bone;
    freezeColors;
    %scatter3(subset(:, 1), subset(:, 2), subset(:, 3))
    % figure; hold on;

    %pcshow(subset(:, 1:3), 'MarkerSize', 18);
    colormap cool;
    min_occ = min(subset(:, 4)) - 0.1;
    max_occ = max(subset(:, 4)) + 0.1;
    [atom_x, atom_y, atom_z] = sphere;
    for ii = 1:size(subset, 1)
        color = ones(size(atom_x)) * (subset(ii, 4) - min_occ) /...
            (max_occ - min_occ);
    %    color(:, 3) = 1.0 - color(:, 1);
        surf(RADIUS * atom_x + subset(ii, 1), ...
             RADIUS * atom_y + subset(ii, 2), ...
             RADIUS * atom_z + subset(ii, 3), ...
             color, 'facealpha', 0.9, 'edgealpha', 0.0);
    end

end

axis equal;
axis tight;
axis off;