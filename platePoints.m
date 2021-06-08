clear all
close all
clc
% IMPORTANT: THIS SCRIPT REQUIRES THE Partial Differential Equation Toolbox.
%% Import the stl file

 x = stlread('NAH_ESM_mesh_light v15.stl');
 trimesh(x);