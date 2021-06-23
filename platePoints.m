clear all
close all
clc

%% Import the stl file

 x = stlread('NAH_ESM_meshA.stl');
 trimesh(x);