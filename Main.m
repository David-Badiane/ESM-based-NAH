close all
clearvars
clc
addpath 'CSV\pressure' 
addpath 'CSV\velocity' 
addpath 'functions'

%% read the csv

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';
[violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = importData(velocityFileName, pressureFileName);


%% Green's functions matrix

% Points on the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints);

% Green's matrices in a cell array ( for each eigenfrequency)
[G_w] = Green_matrix(hologramPoints , virtualPoints , eigenFreqz);

