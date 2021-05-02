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
eigenFreqzRad = 2*pi*eigenFreqz; % convert in [rad/s]

for ii =1:4
hologramInfos{ii} = hologramInfos{ii}*0.001;
violinInfos{ii} = violinInfos{ii}*0.001;
end

%% Green's functions matrix

% choose the mode 
nModes = 2;
alphas = zeros(8,nModes);
rowsNames = { 'NMSE TIK' 'NCC TIK' 'NMSE TSVD' 'NCC TSVD' 'alpha NMSE TIK' 'alpha NCC TIK' 'k NMSE TSVD' 'k NCC TSVD'};
freqzNames = cell(nModes,1);

for ii = 1:nModes
    freqzNames{ii} = ['f', int2str(ii)];
end

for mode = 1:nModes
%conversion from [mm] to [m]


% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Coordinates of the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice, deleteIndexes] = getVirtualPoints(violinInfos,hologramPoints, true); 

% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(mode)]);

%% Inverse problem

omega = eigenFreqzRad(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector
measuredPressure = pressureFields{mode};
meshSize = numel(measuredPressure);
measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed

measuredPressureN = whiteNoise(measuredPressure, 10); % add white gaussian noise to the mesurement

G_p_omega = G_p{1}; % take the Green's function matrix of the chosen mode

q_TSVD = (1/(1i*omega*rho)).*TSVD(G_p_omega, measuredPressure  , 1); % perform the TSVD -> estimate the source strength

q_TIK= (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega, measuredPressure  , 0.5);  % perform the Tikhonov SVD -> estimate the source strength

%% direct problem - green function computation

% Get the normal vectors on the plate surface

gridX = length(unique(platePoints(:, 1)));
gridY = length(unique(platePoints(:, 2)));

X =  violinInfos{1}; 
Y =  violinInfos{2}; 
Z =  violinInfos{3}; 

Z(isnan(Z)) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
%quiver3(X,Y,Z,nx,ny,nz)         % to plot the normal vectors

normalPoints = [reshape(nx', [1024,1]), reshape(ny', [1024,1]), reshape(nz', [1024,1]) ];

% Calculate the derivative of the Green function along the normal direction
[G_v] = normalGradient(virtualPoints, platePoints , [eigenFreqzRad(mode)], normalPoints);

%% direct problem - reconstruction

% Green's matrices of the plate surface - equivalent sources in a cell array (for each eigenfrequency)
G_v_omega = G_v{1}; % Green's function (equivalent points - surface) for the mode

v_TSVD = G_v_omega.'*q_TSVD; % reconstructed velocity with truncated SVD

v_TIK = G_v_omega.'*q_TIK; % reconstructed velocity with Tikhonov

%% error evalutation - velocity 

v_ex = velocityFields{mode};

vel_size = numel(v_ex);

v_ex_vector = reshape( v_ex.', [vel_size, 1]);
v_ex_vector(deleteIndexes,:) = [];
%[NCCv , NMSEv] = errorEvaluation(v_ex_vector, v_TIK.');

%% error evaluation - pressure

p_TIK = 1i*omega*rho*G_p_omega*q_TIK;

[NCCp , NMSEp] = errorEvaluation(measuredPressure, p_TIK);


%% L curve (Tikhonov) or similar 
% the L curve computed with the reconstructed pressure

rangeTIK = [1e-5,100]; % range of value for the regularization parameter
rangeTSVD = [1,64 ]; % range of value for the regularization parameter
numParamsTIK = 1e2;
numParamsTSVD = 64;
%points =  L_Curve(G_p_omega, measuredPressure, rangeTIK, numParamsTIK, rho, omega);
%points =  L_CurveV(G_p_omega, measuredPressure, G_v_omega, v_ex_vector, rangeTIK, numParamsTIK, rho, omega);
[velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexes);

%[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );

%% Recalculation
%[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );
q_TSVD = (1/(1i*omega*rho)).*TSVD(G_p_omega, measuredPressure  ,desiredAlpha(4,2)); % perform the TSVD -> estimate the source strength
q_TIK= (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega, measuredPressure  ,desiredAlpha(2,2));  % perform the Tikhonov SVD -> estimate the source strength
v_TSVD = G_v_omega.'*q_TSVD; % reconstructed velocity with truncated SVD
v_TIK = G_v_omega.'*q_TIK; % reconstructed velocity with Tikhonov


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


surfRecP = reshape( p_TIK , [8, 8]); 

figure(6543)
surf(hologramInfos{1},hologramInfos{2},abs(surfRecP))
title('reconstructed pressure')

figure(6544)
surf(hologramInfos{1},hologramInfos{2},abs(pressureFields{mode}))
title('actual pressure')

alphas(:,mode) = desiredAlpha(:);
end

alphasTable = array2table(alphas, 'rowNames', rowsNames,'variableNames', freqzNames)
