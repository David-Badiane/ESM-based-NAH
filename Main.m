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
eigenFreqzRad = 2*pi*eigenFreqz;

%% Green's functions matrix

% Points on the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Points on the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints);

% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p] = Green_matrix(hologramPoints , virtualPoints , eigenFreqzRad);

% [G_components] = Green_matrixComponents(hologramPoints , virtualPoints ,
% eigenFreqzRad); NOT USEFUL - cancelliamo la funzione ?

%% Inverse problem

% choose the mode 

mode = 2;

omega = eigenFreqzRad(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector

p = pressureFields{mode};
mes_size = numel(p);
p = reshape(abs(p) , [mes_size,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed

p_n = whiteNoise(p); % add white gaussian noise to the mesurement

G = G_p{mode}; % take the Green's function matrix of the chosen mode

q_TSVD = (1/1i*omega*rho).*TSVD(G, p_n , 10); % perform the TSVD -> estimate the source strength

q_TIK = (1/1i*omega*rho).*Tikhonov_SVD(G , p_n , 10);  % perform the Tikhonov SVD -> estimate the source strength

%% direct problem 

ext_v = velocityFields{mode}; % exact velocity 

vel_size = numel(ext_v);

% Get the normal vectors on the plate surface

gridX = length(unique(platePoints(:,1)));
gridY = length(unique(platePoints(:, 2)));
X =  reshape(platePoints(:,1), [gridY, gridX]).'; 
Y =  reshape(platePoints(:,2), [gridY, gridX]).'; 
Z =  reshape(platePoints(:,3), [gridY, gridX]).';

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface

normalPoints = [reshape(nx, [1024,1]), reshape(ny, [1024,1]), reshape(nz, [1024,1]) ];

% Calculate the derivative of the Green function along the normal direction
[G_v] = normalGradient(virtualPoints, platePoints , eigenFreqzRad, normalPoints);

% Green's matrices of the plate surface - equivalent sources in a cell array (for each eigenfrequency)
G_v_omega = G_v{mode}; % Green's function equivalent points - surface  for the mode









