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

%% Green's functions matrix

% choose the mode 

mode = 2;

%conversion from [mm] to [m]
for ii =1:4
hologramInfos{ii} = hologramInfos{ii}*0.001;
violinInfos{ii} = violinInfos{ii}*0.001;
end

% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Coordinates of the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints, true); 


% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(mode)]);

%% Inverse problem

omega = eigenFreqzRad(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector

measuredPressure = pressureFields{mode};
meshSize = numel(measuredPressure);
measuredPressure = reshape(abs(measuredPressure) , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed

measuredPressureN = whiteNoise(measuredPressure, 20); % add white gaussian noise to the mesurement

G_p_omega = G_p{1}; % take the Green's function matrix of the chosen mode

q_TSVD = (1/(1i*omega*rho)).*TSVD(G_p_omega, measuredPressure  , 60); % perform the TSVD -> estimate the source strength

q_TIK= (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega, measuredPressure  , 6.53);  % perform the Tikhonov SVD -> estimate the source strength

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

G_v_omega(isnan(G_v_omega)) = 0; % NaN points in Green's function can't be be used for reconstruction

v_TSVD = G_v_omega*q_TSVD; % reconstructed velocity with truncated SVD

v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov

%% error evalutation - velocity 

v_ex = velocityFields{mode};

vel_size = numel(v_ex);

v_ex_vector = reshape( v_ex.', [vel_size, 1]);

v_ex_vector(deleteIndexesVirt) = [];
% v_ex_vector(isnan(v_ex_vector))=0; %this may cause the metrics to be very low, probably we must delete the NaN from the vector instead
 
[NCCv , NMSEv] = errorEvaluation(v_ex_vector, v_TIK);

%% error evaluation - pressure

p_TIK = 1i*omega*rho*G_p_omega*q_TIK;

[NCCp , NMSEp] = errorEvaluation(measuredPressure, p_TIK);


%% L curve (Tikhonov) or similar 
% the L curve computed with the reconstructed pressure

rangeTIK = [1e-5,1e2 ]; % range of value for the regularization parameter
rangeTSVD = [1,64 ]; % range of value for the regularization parameter
numParamsTIK = 2e2;
numParamsTSVD = 64;
% L_Curve(G_p_omega, measuredPressureN , range, numberParameters, rho, omega);
[velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexesVirt);

%[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );

%% plot of the reconstruced velocity field vs. exact velocity field
v_TSVD_Fin = addNans(virtualPoints, v_TSVD);
v_TIK_Fin = addNans(virtualPoints, v_TIK);
v_ex_Fin = addNans(virtualPoints, v_ex_vector);

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


