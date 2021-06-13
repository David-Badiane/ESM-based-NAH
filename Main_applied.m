
% fv1 = stlread('NAH_ESM_mesh_light v15.stl');
% fv = stlread('NAH_ESM_mesh_A.stl');
% points = sortrows(fv.Points);
% 
% nrows = 35; ncols = 35; fileName = ['grid',int2str(nrows),'x',int2str(ncols)];
% [outMatrix] = downsampling_regular(pts, nrows, ncols, fileName);
% [pts] = downsampling_regular(outMatrix, 128, 128)
 
%% INIT
%{
fv1 = stlread('NAH_ESM_mesh_light v15.stl');
fv = stlread('NAH_ESM_mesh_A.stl');
points = sortrows(fv.Points);

nrows = 45; ncols = 15; fileName = ['grid',int2str(nrows),'x',int2str(ncols)];
[outMatrix] = downsampling_regular(pts, nrows, ncols, fileName);
[pts] = downsampling_regular(outMatrix, 128, 128)
%}
clear all
close all
clc
%% NEAR FIELD ACOUSTIC HOLOGRAPHY - ESM METHOD %%

% folders
baseFolder = pwd;
virtualPointsFolder = [baseFolder,'\VPGrids'];
estimationsFolder = [baseFolder, '\Estimations'];
addpath(genpath('functions'));
addpath 'virtualPointsGrids_old';
addpath('violinMeshes');

%% Virtual Points generator

pts = table2array(readtable('grid128x128Fin.csv'));
controller = 1;
zVal = 0; % <-- lattice
nGrids = 1;

%from last file automatically
filesList = ls(virtualPointsFolder);
filesList(1:2,:) = [];
start = [];
for ii = 1:length(filesList(:,1))
    idx = find(filesList(ii,:)== '.');
    start = [start; eval(filesList(ii,4:idx)) + 1];
end
start = max(start);

% all already set and debugged
for ii = 0:nGrids-1
genVirtualPoints(pts,['VP_',int2str(start+ii)], controller, zVal,virtualPointsFolder);
end

%% f = 377

pX = 128;
pY = 128;

nMics = 8;
nMeas = 8;

nPressureMeas = nMics*nMeas;
nViolinPoints = pX*pY;

rho = 1.2; % [Kg/m3] air density 

eigenFreqz = readtable('eigenfrequencies.csv');
eigenFreqz = table2array(eigenFreqz);
eigenFreqz = 2*pi*eigenFreqz; % convert in [rad/s]

nModes = length(eigenFreqz);

violinMesh = table2array(readtable('grid128x128clean')); 
violinMesh = violinMesh.*0.001; % convert in meter
violinMesh(:,1:2) = -1.*violinMesh(:,1:2);

% get normal points for Green's fxs gradient
X = reshape(violinMesh(:,1), pX, pY);
Y = reshape(violinMesh(:,2), pX, pY);
Z = reshape(violinMesh(:,3), pX, pY);

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
normalPoints = [reshape(nx', [nViolinPoints,1]),...
                reshape(ny', [nViolinPoints,1]),...
                reshape(nz', [nViolinPoints,1]) ];

pressureData = table2array(readtable('pressureData.csv'));
measuredPressure = pressureData(:,3:end);
hologramDistance = 0.02;
zHologram = max(Z(:)) + hologramDistance; 
hologramPoints =  [0.001.*pressureData(:,1:2), zHologram*ones(size(pressureData(:,1)))] ; 

velocityData = table2array(readtable('velocityData.csv'));   
velocityData(:,1:2) = 0.001.*velocityData(:,1:2);
nEqSourceGrids = 1;
 

%% COMPUTATION LOOP 

metricsTSVD = [];
metricsTIK = [];

alphaTSVD = [];
alphaTIK = [];

 for ii = 1:length(eigenFreqz)
     
  nccTIKs = [];
  nccTSVDs = [];
  qTSVDs = cell(1,2);
  qTIKs = cell(1,2);
     
     for jj = 1:nEqSourceGrids
        % choose virtual points grid     
        virtualPtsFilename = ['VP_', int2str(jj), '.csv'];
        virtualPoints = table2array(readtable(virtualPtsFilename)) ;
        deleteIndexes = find(isnan(virtualPoints(:,3)));

        % compute Green functions matrix ( hologram 2 virtual points ) 
        [G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , omega );
        G_p_omega = G_p{1}; 

        % compute gradient Green functions matrix 
        % ( virtual points 2 reconstruction plane = violin )
        [G_v] = normalGradient(virtualPoints, violinPoints , omega, normalPoints);
        G_v_omega = G_v{1};


        % APPROACH 1 ) L curve solutions
        % 1) Inverse - individuation of regularization parameter (lambda) 
        [U,s,V] = csvd (G_p_omega);
        figure(2)
        lambda_l = l_curve (U,s,measuredPressure);
        figure(3)
        k_l = l_curve (U,s,measuredPressure,'tsvd');

        % 2) Inverse - calculation of equivalent sources weights
        Lq_TIK =  (1/(1i*omega*rho)).* tikhonov (U,s,V,measuredPressure,lambda_l);
        Lq_TSVD = (1/(1i*omega*rho)).* tsvd (U,s,V,measuredPressure,k_l); % perform the TSVD -> estimate the source strength

        % 3) Direct - solutions ( velocities ) 
        Lv_TSVD = G_v_omega*Lq_TSVD; % reconstructed velocity with truncated SVD
        Lv_TIK = G_v_omega*Lq_TIK; % reconstructed velocity with Tikhonov

        % APPROACH 2) metrics parametrization
        rangeTIK = [0,100]; % range of value for the regularization parameter
        rangeTSVD = [1,32 ]; % range of value for the regularization parameter
        numParamsTIK = 1e2;
        numParamsTSVD = 64;

        % !!! need to revisitate the function and compute the interpolation !!!
        [velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, ... 
        measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD,...
        numParamsTIK, numParamsTSVD,     omega,      rho, deleteIndexes);

        % recalculation with best metrics approach
        q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressure  ,desiredAlpha(4,2)); % perform the TSVD -> estimate the source strength
        q_TIK= (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressure,desiredAlpha(2,2));  % perform the Tikhonov SVD -> estimate the source strength
        v_TSVD = G_v_omega*q_TSVD; % reconstructed velocity with truncated SVD
        v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov

        % !! TO WRITE !! 
        % COMPARE APPROACH 1 & APPROACH 2 by calculating NCC of
        % interpolation
        % write a fx  like : 
        % [ncc_LTIK, ncc_MTIK, ncc_LTSVD, ncc_MTSVD] = function(v_TSVD, v_TIK,v_LTSVD, v_LTIK, velocityGroundtruth)

        % store results
        nccTIKs = [nccTIKs; ncc_LTIK, ncc_MTIK];
        nccTSVDs = [nccTSVDs; ncc_LTSVD, ncc_MTSVD];

        qTSVDs{ii,1} = q_TSVD;
        qTSVDs{ii,2} = Lq_TSVD;
        qTIKs{ii,1} = q_TIK;
        qTIKs{ii,2} = Lq_TIK;
    end 

    % !!! riflettere su architettura dati; trovarne semplice e effettiva !!!
    [TSVDval, TSVDloc] = max(NCCs(:,2), [], 'all');
    [TIKval, TIKloc] = max(NCCs(:,1), [], 'all');
    
    % store best NCC values 
    metricsTSVD = [metricsTSVD; TSVDval];
    metricsTIK = [metricsTIK;   TIKval];
    
    alphaTSVD = [alphaTSVD;   ];
    alphaTIK = [alphaTIK;    ];
    % .... .....
    
    % once individuated the best, let's represent them
    v_TSVD_Fin = addNans(platePoints, v_TSVD);
    v_TIK_Fin = addNans(platePoints, v_TIK);
    v_ex_Fin = addNans(platePoints, v_ex_vector);

    surfVelRecTSVD = reshape( v_TSVD_Fin , [pY, pX]).'; 
    surfVelRecTIK = reshape( v_TIK_Fin , [pY, pX]).'; 
    
    p_TIK = 1i*omega*rho*G_p_omega*q_TIK;
    surfRecP = reshape( p_TIK , [nMics, nMeas]); 

    % velocity
    figure(100) 
    subplot 121
    surf(X, Y, abs(surfVelRecTSVD))
    title('TSVD velocity')
    subplot 122
    surf(X, Y, abs(surfVelRecTIK))
    title('Tik velocity')

    figure(101)
    % pressure
    subplot(121)
    surf(hologramInfos{1},hologramInfos{2},abs(surfRecP))
    title('reconstructed pressure')
    subplot(122)
    surf(hologramInfos{1},hologramInfos{2},abs(pressureFields{mode}))
    title('actual pressure')
    
    % save estiamtion results
    cd(estimationsFolder)
    writeMat2File(data, 'filename.csv', {'' ''}, 2, true);
    cd(baseFolder)
    
    msg = 'press something to go on !!'
    p = input(msg);
    
 end
 
 % plot a figure with the NCC for the various frequencies
 
 % check regularization parameters as frequency varies