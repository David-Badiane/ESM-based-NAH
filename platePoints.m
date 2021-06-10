clear all
close all
clc

%% Import the stl file

 x = stlread('NAH_ESM_mesh_light v15.stl');
 trimesh(x);