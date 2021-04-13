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

% Green's matrices in a cell array (for each eigenfrequency)
[G_w] = Green_matrix(hologramPoints , virtualPoints , eigenFreqz);

%% Inverse problem

% choose the mode 

mode = 2;

% measurements of pressure vector

p = pressureFields{mode};
mes_size = numel(p);
p = reshape(p , [mes_size,1]); % convert the measurement matrix into an array

p = abs(p); % the magnitude of pressure is needed

p_n = whiteNoise(p); % add white gaussian noise to the mesurement

G = G_w{mode}; % take the Green's function matrix of the chosen mode




