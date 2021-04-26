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

for mode = 1:5

% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Coordinates of the plate surface
platePoints = violinInfos{4};

% Virtual sources points
[virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints, false); % to change: here HologramPOints are not used

%conversion from [mm] to [m]
hologramPoints = hologramPoints*0.001;
virtualPoints = virtualPoints*0.001;
platePoints = platePoints*0.001;


% Green's matrices of the hologram-equivalent sources in a cell array (for each eigenfrequency)
[G_p] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(mode)]);

%% Inverse problem

omega = eigenFreqzRad(mode);

rho = 1.2; % [Kg/m3] air density 

% measurements of pressure vector

measuredPressure = pressureFields{mode};
meshSize = numel(measuredPressure);
measuredPressure = reshape(abs(measuredPressure) , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed

measuredPressureN = whiteNoise(measuredPressure); % add white gaussian noise to the mesurement

G_p_omega = G_p{1}; % take the Green's function matrix of the chosen mode

G_p_omega(isnan(G_p_omega)) = 0; 

q_TSVD = (1/(1i*omega*rho)).*TSVD(G_p_omega, measuredPressureN  , 10); % perform the TSVD -> estimate the source strength

q_TIK= (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega , measuredPressureN  , 5);  % perform the Tikhonov SVD -> estimate the source strength

%% direct problem - green function computation

% Get the normal vectors on the plate surface

gridX = length(unique(platePoints(:, 1)));
gridY = length(unique(platePoints(:, 2)));

X =  violinInfos{1}; 
Y =  violinInfos{2}; 
Z =  violinInfos{3}; 

X =  X*0.001; 
Y =  Y*0.001; 
Z =  Z*0.001; 

Z(isnan(Z)) = 0;

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
%quiver3(X,Y,Z,nx,ny,nz)         % to plot the normal vectors

normalPoints = [reshape(nx, [1024,1]), reshape(ny, [1024,1]), reshape(nz, [1024,1]) ];

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

v_ex_vector = reshape( v_ex, [vel_size, 1]);

v_ex_vector(isnan(v_ex_vector))=0; %this may cause the metrics to be very low, probably we must delete the NaN from the vector instead
 
[NCCv , NMSEv] = errorEvaluation(v_ex_vector, v_TIK);

%% error evaluation - pressure

p_TIK = 1i*omega*rho*G_p_omega*q_TIK;

[NCCp , NMSEp] = errorEvaluation(measuredPressure, p_TIK);

%% L curve (Tikhonov)
% the L curve computed with the reconstructed pressure

range = [1e-10, 1e2]; % range of value for the regularization parameter

numberParameters = 2e2;

% L_Curve(G_p_omega, measuredPressureN , range, numberParameters, rho, omega);

%% plot of the reconstruced velocity field vs. exact velocity field

surfVelRecTSVD = reshape( v_TSVD , [gridY, gridX]).'; 
surfVelRecTIK = reshape( v_TIK , [gridY, gridX]).'; 

surfVelRecTSVD( surfVelRecTSVD == 0) = NaN; % apply a mask on the reconstructed velocity
surfVelRecTIK( surfVelRecTIK == 0) = NaN;

figure(600) 
subplot 131
surf(X, Y, abs(v_ex)./max(abs(v_ex)))
title('Exact velocity')
subplot 132
surf(X, Y, abs(surfVelRecTSVD)./max(abs(surfVelRecTSVD)))
title('TSVD velocity')
subplot 133
surf(X, Y, abs(surfVelRecTIK)./max(abs(surfVelRecTIK)))
title('Tik velocity')

%% plot of the reconstruced pressure field vs. exact pressure field


surfRecP = reshape( p_TIK , [8, 8]); 

figure(6543)
surf(hologramInfos{1},hologramInfos{2},abs(surfRecP))
title('reconstructed pressure')

figure(6544)
surf(hologramInfos{1},hologramInfos{2},abs(pressureFields{mode}))
title('actual pressure')
end

