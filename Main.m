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

% Points on the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints);

% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p] = Green_matrix(hologramPoints , virtualPoints , eigenFreqz);

% Green's matrices of the plate surface - equivalent sources in a cell array (for each eigenfrequency)
[G_v] = Green_matrix(platePoints , virtualPoints , eigenFreqz);

%% Inverse problem

% choose the mode 

mode = 2;

omega = eigenFreqz(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector

p = pressureFields{mode};
mes_size = numel(p);
p = reshape(p , [mes_size,1]); % convert the measurement matrix into an array

p = abs(p); % the magnitude of pressure is needed

p_n = whiteNoise(p); % add white gaussian noise to the mesurement

G = G_p{mode}; % take the Green's function matrix of the chosen mode

q_TSVD = (1/1i*omega*rho).*TSVD(G, p_n , 10); % perform the TSVD -> estimate the source strength

q_TIK = (1/1i*omega*rho).*Tikhonov_SVD(G , p_n , 10);  % perform the Tikhonov SVD -> estimate the source strength

%% direct problem 

ext_v = velocityFields{mode}; % exact velocity 

vel_size = numel(ext_v);




