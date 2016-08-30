%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to extract map dimensions.
% 
% Author: David Fridovich-Keil (dfk@eecs.berkeley.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Load an AtomMap.
data = csvread('../saved_maps/lbl_500cm.csv');

%% Extract dimensions.
max_x = max(data(:, 1));
min_x = min(data(:, 1));
max_y = max(data(:, 2));
min_y = min(data(:, 2));
max_z = max(data(:, 3));
min_z = min(data(:, 3));

fprintf('Map dimensions: %f x %f x %f\n', ...
    max_x - min_x, max_y - min_y, max_z - min_z);
fprintf('Volume = %f\n', (max_x - min_x) * (max_y - min_y) * (max_z - min_z));