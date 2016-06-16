%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization script. Load up an AtomMap, allow user to select a box, and
% compute the implicit surface in that region.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Load an AtomMap.
data = csvread('../saved_maps/map.csv');

%% Display a 2D scatterplot of the data projected onto the x-y plane.
figure;
scatter(data(:, 1), data(:, 2), ...
       'MarkerEdgeColor', [0 .5 .5], ...
       'MarkerFaceColor', [0 .7 .7], ...
       'LineWidth', 1.5);

%% Let the user select two points (to define a bounding box -- LL then TR).
fprintf('Please select two points and then press ENTER...\n');
[pts_x, pts_y] = getpts

subset = data(data(:, 1) > pts_x(1) & data(:, 1) < pts_x(2) & ...
              data(:, 2) > pts_y(1) & data(:, 2) < pts_y(2), :);

z_min = min(subset(:, 3)) - 1;
z_max = max(subset(:, 3)) + 1;

%% Create a kdtree.
kdtree = KDTreeSearcher(subset(:, 1:3));

if length(pts_x) == 2
%% Allocate a mesh over the specified area.
    resolution = 0.5;
    [X, Y, Z] = meshgrid(pts_x(1):resolution:pts_x(2), ...
                         pts_y(1):resolution:pts_y(2), ...
                         z_min:resolution:z_max);

    size(X)
%% Evaluate GP at each voxel in mesh.
    [sdfs, vars] = arrayfun(@(x, y, z) InterpolateGP(x, y, z, kdtree, subset), X, Y, Z);
    figure; set(gca, 'fontsize', 16);
    isosurface(X, Y, Z, sdfs, 0, vars);
end

axis equal;