close all
clearvars
clc
addpath 'CSV\pressure' 
addpath 'CSV\velocity' 
addpath 'functions'
addpath 'functions\regu'

%% read the csv

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';
resampleX = 8;
resampleY = 8;
[violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = importData(velocityFileName, pressureFileName, resampleX, resampleY);
eigenFreqzRad = 2*pi*eigenFreqz; % convert in [rad/s]

%conversion from [mm] to [m]
for ii =1:4
hologramInfos{ii} = hologramInfos{ii}*0.001;
violinInfos{ii} = violinInfos{ii}*0.001;
end

%% Setup of global variables 

% choose the mode 
nModes = 10;
alphas = zeros(8,nModes);
rowsNames = { 'NMSE TIK' 'NCC TIK' 'NMSE TSVD' 'NCC TSVD' 'alpha NMSE TIK' 'alpha NCC TIK' 'k NMSE TSVD' 'k NCC TSVD'};
freqzNames = cell(nModes,1);

for ii = 1:nModes
    freqzNames{ii} = ['f', int2str(ii)];
end

rho = 1.2; % [Kg/m3] air density 

% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Coordinates of the plate surface
platePoints = violinInfos{4};

% 1) Get the normal vectors on the plate surface
gridX = length(unique(platePoints(:, 1)));
gridY = length(unique(platePoints(:, 2)));

X =  violinInfos{1};  Y =  violinInfos{2};  Z =  violinInfos{3}; 
Z(isnan(Z)) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
nNormPoints = gridX*gridY;                 
normalPoints = [reshape(nx', [nNormPoints,1]),...
                reshape(ny', [nNormPoints,1]),...
                reshape(nz', [nNormPoints,1]) ];

%% START
for mode = 1:nModes

% Setup of local variables
omega = eigenFreqzRad(mode);

% Virtual sources points
virtnames = {'scale' 'offset' 'cutX' 'cutY' 'borderX' 'borderY' 'active' 'plotData'};
params = {1.15,5,1,3,3,true};
[virtualPoints, lattice,deleteIndexes] = getVirtualPoints(violinInfos,hologramPoints, params, true)

% pressure vector setup
measuredPressure = pressureFields{mode};
meshSize = numel(measuredPressure);
measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed
measuredPressureN = whiteNoise(measuredPressure,-10); % add white gaussian noise to the mesurement

% velocity vector setup
v_ex = velocityFields{mode};
vel_size = numel(v_ex);
v_ex_vector = reshape( v_ex.', [vel_size, 1]);
v_ex_vector(deleteIndexes,:) = [];
check = find(abs(v_ex_vector) == 0);
v_ex_vector(check) = v_ex_vector(check-1);


%% Inverse problem (individuation of ES weights)

% 1) Green's functions matrix
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(mode)]);
G_p_omega = G_p{1}; % take the Green's function matrix of the chosen mode

% 2) individuation of regularization parameter (lambda)
[U,s,V] = csvd (G_p_omega);
figure(2)
lambda_l = l_curve (U,s,measuredPressure);
figure(3)
k_l = l_curve (U,s,measuredPressure,'tsvd');

% 3) calculation of equivalent sources weights
q_TIK =  (1/(1i*omega*rho)) * tikhonov (U,s,V,measuredPressure,lambda_l);
q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V,measuredPressure,k_l); % perform the TSVD -> estimate the source strength

%% direct problem - reconstruction

% 1) Calculate the gradient of G along the normal vectors
[G_v] = normalGradient(virtualPoints, platePoints , [eigenFreqzRad(mode)], normalPoints);
G_v_omega = G_v{1};

% 2) Reconstruction %%%%%% INUTILE POSSIAMO TOGLIERLO %%%%%%%
v_TSVD = G_v_omega*q_TSVD; % reconstructed velocity with truncated SVD
v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov

%% compute metrics to choose best solution 
% the L curve computed with the reconstructed pressure

rangeTIK = [0,100]; % range of value for the regularization parameter
rangeTSVD = [1,32 ]; % range of value for the regularization parameter
numParamsTIK = 1e2;
numParamsTSVD = 64;
% points =  L_Curve(G_p_omega, measuredPressure, rangeTIK, numParamsTIK, rho, omega);
% points =  L_CurveV(G_p_omega, measuredPressure, G_v_omega, v_ex_vector, rangeTIK, numParamsTIK, rho, omega);
[velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexes);
%[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );

%% Recalculation of solution with best metrics %%%% CALCULATION OF SOLUTION %%%%%

%[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );
 q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressure  ,desiredAlpha(4,2)); % perform the TSVD -> estimate the source strength
 q_TIK= (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressure,desiredAlpha(2,2));  % perform the Tikhonov SVD -> estimate the source strength

 v_TSVD = G_v_omega*q_TSVD; % reconstructed velocity with truncated SVD
 v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov

%% plot of the reconstruced velocity field vs. exact velocity field

v_TSVD_Fin = addNans(platePoints, v_TSVD);
v_TIK_Fin = addNans(platePoints, v_TIK);
v_ex_Fin = addNans(platePoints, v_ex_vector);

surfVelRecTSVD = reshape( v_TSVD_Fin , [gridY, gridX]).'; 
surfVelRecTIK = reshape( v_TIK_Fin , [gridY, gridX]).'; 
surfVel = reshape( v_ex_Fin , [gridY, gridX]).';

figure(600) 
subplot 131
surf(X, Y, abs(surfVel))
title('Exact velocity')
subplot 132
surf(X, Y, abs(surfVelRecTSVD))
title('TSVD velocity')
subplot 133
surf(X, Y, abs(surfVelRecTIK))
title('Tik velocity')

%% plot of the reconstruced pressure field vs. exact pressure field
p_TIK = 1i*omega*rho*G_p_omega*q_TIK;

surfRecP = reshape( p_TIK , [resampleX, resampleY]); 

figure(6543)
subplot(121)
surf(hologramInfos{1},hologramInfos{2},abs(surfRecP))
title('reconstructed pressure')
subplot(122)
surf(hologramInfos{1},hologramInfos{2},abs(pressureFields{mode}))
title('actual pressure')

% alphas(:,mode) = desiredAlpha(:);
end

alphasTable = array2table(alphas, 'rowNames', rowsNames,'variableNames', freqzNames)


%% SEE Virtual Points grids
figure(150)
for ii = 1:20
  virtualPoints = table2array(readtable(['VP_', int2str(ii),'.csv']));  
  plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 10);
  pause(3);
end
