%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization script. Load up an AtomMap, allow user to select a box, and
% compute the implicit surface in that region.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% (1) Load an AtomMap.
data = csvread('../saved_maps/map.csv');
z_min = min(data(:, 3));
z_max = max(data(:, 3));

%% (2) Create a kdtree.
kdtree = KDTreeSearcher(data(:, 1:3));

%% (3) Display a 2D scatterplot of the data projected onto the x-y plane.
figure;
scatter(data(:, 1), data(:, 2), ...
       'MarkerEdgeColor', [0 .5 .5], ...
       'MarkerFaceColor', [0 .7 .7], ...
       'LineWidth', 1.5);

%% (4) Let the user select two points (to define a bounding box).
fprintf('Please select two points and then press ENTER...\n');
[pts_x, pts_y] = getpts

if length(pts_x) == 2
%% (5) Allocate a mesh over the specified area.
    resolution = 0.3;
    [X, Y, Z] = meshgrid(pts_x(1):resolution:pts_x(2), ...
                         pts_y(1):resolution:pts_y(2), ...
                         z_min:resolution:z_max);

    size(X)
%% (6) Evaluate GP at each voxel in mesh.
    sdfs = arrayfun(@(x, y, z) InterpolateGP(x, y, z, kdtree, data), X, Y, Z);
    figure;
    isosurface(X, Y, Z, sdfs, 0);
end