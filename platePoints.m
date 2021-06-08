clear all
close all
clc
% IMPORTANT: THIS SCRIPT REQUIRES THE Partial Differential Equation Toolbox.
%% Import the stl file

figure
model = createpde(3);
plateStl = importGeometry(model, 'NAH_ESM_mesh_light v15.stl');
pdegplot(plateStl)