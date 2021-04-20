close all
clearvars
clc
addpath 'CSV\pressure' 
addpath 'CSV\velocity' 
addpath 'functions'

%% read the csv

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';
[violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz, hologramInfos8] = importData(velocityFileName, pressureFileName);
eigenFreqzRad = 2*pi*eigenFreqz; % convert in [rad/s]

%% Green's functions matrix

% Points on the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Points on the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints);

% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(1)]);

% [G_components] = Green_matrixComponents(hologramPoints , virtualPoints ,
% eigenFreqzRad); NOT USEFUL - cancelliamo la funzione ?

%% Inverse problem

% choose the mode 

mode = 1;

omega = eigenFreqzRad(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector

p = pressureFields{mode};
mes_size = numel(p);
p = reshape(abs(p) , [mes_size,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed

p_n = whiteNoise(p); % add white gaussian noise to the mesurement

G_p_omega = G_p{mode}; % take the Green's function matrix of the chosen mode

% TO DO: REMEBER TO ADD THE NOISE INTO THE FOLLOWING
q_TSVD = (1/1i*omega*rho).*TSVD(G_p_omega, p_n , 10); % perform the TSVD -> estimate the source strength

q_TIK = (1/1i*omega*rho).*Tikhonov_SVD(G_p_omega , p_n , 100);  % perform the Tikhonov SVD -> estimate the source strength

%% direct problem - green function computation

% Get the normal vectors on the plate surface

gridX = length(unique(platePoints(:,1)));
gridY = length(unique(platePoints(:, 2)));

X =  reshape(platePoints(:,1), [gridY, gridX]).'; 
Y =  reshape(platePoints(:,2), [gridY, gridX]).'; 
Z =  reshape(platePoints(:,3), [gridY, gridX]).';

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface

normalPoints = [reshape(nx, [1024,1]), reshape(ny, [1024,1]), reshape(nz, [1024,1]) ];

normalPoints(isnan(normalPoints)) = 0; % surfnorm generates NaN and we need to eliminate such entries to perform the next operations

% Calculate the derivative of the Green function along the normal direction
[G_v] = normalGradient(virtualPoints, platePoints , [eigenFreqzRad(1)], normalPoints);

%% direct problem - reconstruction

% Green's matrices of the plate surface - equivalent sources in a cell array (for each eigenfrequency)
G_v_omega = G_v{mode}; % Green's function equivalent points - surface  for the mode

G_v_omega(isnan(G_v_omega)) = 0;
v_TSVD = G_v_omega*q_TSVD; % TO DO: check the zeros, they might result from NaN, so better convert back in NaN
v_TIK = G_v_omega*q_TIK;

%% error evalutation

v_ex = velocityFields{mode}; % exact velocity IMPORTANT!!!!!!!! ASK THE UNITY OF MEASURE HERE BECAUSE THE VALUES ARE TOO SMALL

vel_size = numel(v_ex);

v_ex_vector = reshape( v_ex, [vel_size, 1]);

v_ex_vector(isnan(v_ex_vector))=0;

NCC = (v_TSVD'*v_ex_vector)/(norm(v_ex_vector)*norm(v_TSVD));

% for the L curve we need to compute the reconstructed pressure

%% plot of the reconstruced velocity field vs. exact velocity field

surfVelRecTSVD = reshape( v_TSVD , [gridY, gridX]).'; 
surfVelRecTIK = reshape( v_TIK , [gridY, gridX]).'; 

figure(600) % TO DO! applicare una maschera al contorno per rimuovere i punti non del violino
subplot 131
surf(X,Y,abs(v_ex))
title('Exact velocity')
subplot 132
surf(X, Y, abs(surfVelRecTSVD))
title('TSVD velocity')
subplot 133
surf(X, Y, abs(surfVelRecTIK))
title('Tik velocity')

%% to do per martedì: ricostruire la pressione

recP = G_p_omega*q_TIK;
surfRecP = reshape( recP , [8, 8]); 
figure(6543)
surf(hologramInfos8{1},hologramInfos8{2},abs(surfRecP))
%%

figure(565)

surf(X,Y,Z);
zlim([0,100]);