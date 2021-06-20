close all
clearvars
clc
VPfilename = 'VPGrids';
addpath(genpath('CSV')); 
addpath(genpath('functions'));
baseFolder = pwd;
virtualPointsFolder = [baseFolder,'\',VPfilename];%VP_Grids'];
addpath(virtualPointsFolder)


%% read the csv

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';
nMics = 8;
nMeas = 8;

% import synthetic data
[violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = ...
    importData(velocityFileName, pressureFileName, nMics, nMeas);
eigenFreqzRad = 2*pi*eigenFreqz; % convert in [rad/s]

%conversion from [mm] to [m]
for ii =1:4
        hologramInfos{ii} = hologramInfos{ii}*0.001;
        violinInfos{ii} = violinInfos{ii}*0.001;
end

%% Setup of global variables 

nModes = 10;
alphas = zeros(8,nModes); % array of regolarizaation parameters and error metrics
rowsNames = { 'NMSE TIK' 'NCC TIK' 'NMSE TSVD' 'NCC TSVD' 'alpha NMSE TIK' 'alpha NCC TIK' 'k NMSE TSVD' 'k NCC TSVD'};
freqzNames = cell(nModes,1);

for ii = 1:nModes
    freqzNames{ii} = ['f', int2str(ii)];
end

rho = 1.2; % [Kg/m3] air density 

% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};

% Coordinates of the plate surface
violinMesh = violinInfos{4};

pX = length(unique(violinMesh(:, 1)));
pY = length(unique(violinMesh(:, 2)));
nNormPoints = pX*pY;                 

X =  violinInfos{1};  Y =  violinInfos{2};  Z =  violinInfos{3}; 
Z(isnan(Z)) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
                                 
                                 
normalPoints = [reshape(nx, [nNormPoints,1]),...
                reshape(ny, [nNormPoints,1]),...
                reshape(nz, [nNormPoints,1]) ];

%% START
nEqSourceGrids = 1;
cd(virtualPointsFolder);
gridTablesNames = {'grid n.', 'zVal', 'lambda_L', 'k_L', 'nmseTSVD_L', 'nccTSVD_L',...
                    'nmseTIK_L','nccTIK_L', 'lambda_nmse_M', 'k_nmse_M', ...
                    'k_ncc_M', 'lambda_ncc_M', 'nmseTSVD_M', 'nccTSVD_M', ...
                    'nmseTIK_M', 'nccTIK_M', };
dataCell = cell(nModes,1);
ZreguFreq = cell(nModes,1);


for mode = 1:nModes
    tStart = tic;
    % Setup of local variables
    omega = eigenFreqzRad(mode); % current eigenfreq mode

    % pressure vector setup
    measuredPressure = pressureFields{mode};
    meshSize = numel(measuredPressure);
    measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed
    measuredPressure = whiteNoise(measuredPressure, 10); % add white gaussian noise to the mesurement

    % velocity vector setup
    v_ex = velocityFields{mode};   
    v_ex_vector = reshape( v_ex.', [numel(v_ex), 1]); 
    v_ex_vector(isnan(v_ex_vector)) = 0;
    deleteIndexes = find(isnan(violinInfos{4}(:,3)));
    v_ex_vector(deleteIndexes,:) = [];
 
    xAx = unique(hologramPoints(:,1));
    yAx = unique(hologramPoints(:,2));
    xStep = abs(xAx(end) - xAx(end-1));
    yStep = abs(yAx(end) - yAx(end-1));
    
    zCenter = -0.15*min([xStep, yStep]);
    nZpoints = 1;
    zSearch = 0;
    transposeGrids = true;
    plotData = true;
    experimentalData = false;
    
    
    [reguData, ZreguDatas] = getBestGrid(nEqSourceGrids, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                           nZpoints, zCenter, zSearch, X, Y , v_ex_vector,...
                           rho, pX, pY,gridTablesNames , transposeGrids, plotData, experimentalData);
    
         
    tempTable = array2table( reguData , 'VariableNames',gridTablesNames);
       
    dataCell{mode} = tempTable;
    ZreguFreq{ii} = ZreguDatas;  
    disp(toc(tStart));

    
    v_ex = velocityFields{mode};   
    v_ex_vector = reshape( v_ex.', [numel(v_ex), 1]); 
    v_ex_vector(isnan(v_ex_vector)) = 0;
    deleteIndexes = find(isnan(violinInfos{4}(:,3)));
    v_ex_vector(deleteIndexes,:) = [];
    for gridN = 2:nEqSourceGrids
        filename = ['VP_', int2str(gridN)];
        [virtualPoints, lattice, deleteIndexes] = getVirtualPoints(violinInfos, hologramPoints, filename, true);
        virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), - 3*zCenter*ones(size(virtualPoints(:,3)))];
        %% Inverse problem (individuation of ES weights)

        % 1) Green's functions matrix
        [G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , [eigenFreqzRad(mode)]);
        G_p_omega = G_p{1}; % take the Green's function matrix of the chosen mode

        % 2) individuation of regularization parameter (lambda)
        [U,s,V] = csvd (G_p_omega);
        lambda_l = l_curve (U,s,measuredPressure);
        k_l = l_curve (U,s,measuredPressure,'tsvd');

        % 3) calculation of equivalent sources weights
        Lq_TIK =  (1/(1i*omega*rho)).* tikhonov (U,s,V,measuredPressure,lambda_l);
        Lq_TSVD = (1/(1i*omega*rho)).* tsvd (U,s,V,measuredPressure,k_l); % perform the TSVD -> estimate the source strength

        %% direct problem - reconstruction

%         % 1) Calculate the gradient of G along the normal vectors
        [G_v] = normalGradient(virtualPoints, violinMesh , eigenFreqzRad(mode), normalPoints);
        G_v_omega = G_v{1};

        % 2) Reconstruction 
        Lv_TSVD = -G_v_omega*Lq_TSVD; % reconstructed velocity with truncated SVD
        Lv_TIK = -G_v_omega*Lq_TIK; % reconstructed velocity with Tikhonov

        %% compute metrics to choose best solution 
        % the L curve computed with the reconstructed pressure

%         rangeTIK = [0,100]; % range of value for the regularization parameter
%         TSVDup = min([length(G_p_omega(1,:)), length(G_p_omega(:,1))]);
%         rangeTSVD = [1,TSVDup];
%         numParamsTIK = 1e2;
%         numParamsTSVD = 64;
%         [velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexes);
%         %[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );

        %% Recalculation of solution with best metrics %%%% CALCULATION OF SOLUTION %%%%%

        %[pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressure , omega , rho , rangeTIK, rangeTSVD , numParamsTIK, numParamsTSVD   );
%          q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressure  ,desiredAlpha(4,2)); % perform the TSVD -> estimate the source strength
%          q_TIK= (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressure,desiredAlpha(2,2));  % perform the Tikhonov SVD -> estimate the source strength
% 
%          v_TSVD = G_v_omega*q_TSVD; % reconstructed velocity with truncated SVD
%          v_TIK = G_v_omega*q_TIK; % reconstructed velocity with Tikhonov

        %% plot of the reconstruced velocity field vs. exact velocity field

%         v_TSVD_Fin = addNans(violinMesh, v_TSVD); 
%         v_TIK_Fin = addNans(violinMesh, v_TIK);
        v_ex_Fin = addNans(violinMesh, v_ex_vector);

        Lv_TSVD_Fin = addNans(violinMesh, Lv_TSVD); 
        Lv_TIK_Fin = addNans(violinMesh, Lv_TIK);

%         surfVelRecTSVD = reshape( v_TSVD_Fin , [pY, pX]).'; 
%         surfVelRecTIK = reshape( v_TIK_Fin , [pY, pX]).'; 
        surfVel = reshape( v_ex_Fin , [pY, pX]).';

        LsurfVelRecTSVD = reshape( Lv_TSVD_Fin , [pY, pX]).'; 
        LsurfVelRecTIK = reshape( Lv_TIK_Fin , [pY, pX]).'; 

%         figure(600) 
%         
%         subplot 311
%         surf(X, Y, abs(surfVel)); view(2);
%         title('Exact velocity')
%         subplot 312
%         surf(X, Y, abs(surfVelRecTSVD));view(2);
%         title('TSVD velocity')
%         subplot 313
%         surf(X, Y, abs(surfVelRecTIK));view(2);
%         title('Tik velocity')
%         sgtitle(['M method', filename]);
%         
        
        figure(601) 
        subplot 311
        surf(X, Y, abs(surfVel)); view(2);
        title('Exact velocity')
        subplot 312
        surf(X, Y, abs(LsurfVelRecTSVD));view(2);
        title('TSVD velocity')
        subplot 313
        surf(X, Y, abs(LsurfVelRecTIK));view(2);
        title('Tik velocity')
        sgtitle(['L curve ', filename]);
        pause(2);
        %% plot of the reconstruced pressure field vs. exact pressure field
%         p_TIK = 1i*omega*rho*G_p_omega*q_TIK;
% 
%         surfRecP = reshape( p_TIK , [nMics, nMeas]); 
% 
% %         figure(602)
% %         subplot(121)
% %         surf(hologramInfos{1},hologramInfos{2},abs(surfRecP))
% %         title('reconstructed pressure')
% %         subplot(122)
% %         surf(hologramInfos{1},hologramInfos{2},abs(pressureFields{mode}))
% %         title('actual pressure')
% 
%         % alphas(:,mode) = desiredAlpha(:);
%         alphasTable = array2table(alphas, 'rowNames', rowsNames,'variableNames', freqzNames);
    end
end
dataStruct = cell2struct(dataCell, freqzNames, 1);

%% SEE Virtual Points grids

cd(baseFolder)
filesList = ls('VPGrids');
filesList(1:2,:) = [];
cd('VPGrids')
figure(150)

for ii = 2

  idx = find(filesList(ii,:)== '.');
  virtualPoints = table2array(readtable(filesList(ii,1:idx+3))); 
%   writeMat2File(virtualPoints, filesList(ii,1:idx+3), {'x' 'y' 'z'}, 3, true);
  plot3(virtualPoints(:,2), virtualPoints(:,1),  virtualPoints(:,3), '.', 'markerSize', 10);
  hold on;
  plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3));
  hold off;
  pause(1);
  
end
