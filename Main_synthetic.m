close all
clearvars
clc
VPfilename = 'VPGrids';
addpath(genpath('CSV')); 
addpath(genpath('functions'));

baseFolder = pwd;
estimationsFolder = [baseFolder, '\Estimations\Synthetic'];
virtualPointsFolder = [baseFolder,'\',VPfilename];%VP_Grids'];
addpath(virtualPointsFolder)
addpath(estimationsFolder)

%% read the synthetic data - csv
nMics = 8;
nMeas = 8;

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';

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

nModes = 20;
rho = 1.2; % [Kg/m3] air density 

freqzNames = cell(nModes,1);
for ii = 1:nModes
    freqzNames{ii} = ['f', int2str(ii)];
end

% Coordinates of the pressure field (hologram)
hologramPoints = hologramInfos{4};
% Coordinates of the plate surface
violinMesh = violinInfos{4};

pX = length(unique(violinMesh(:, 1)));
pY = length(unique(violinMesh(:, 2)));
nViolinPoints = pX*pY;                 
X =  violinInfos{1};  Y =  violinInfos{2};  Z =  violinInfos{3}; 

%% Calculation of the normal vectors wrt the violin mesh

Z(isnan(Z)) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
                                                                  
normalPoints = [reshape(nx, [nViolinPoints,1]),...
                reshape(ny, [nViolinPoints,1]),...
                reshape(nz, [nViolinPoints,1]) ];

%% STRUCT AND TABLE INIT

nEqSourceGrids = 1;
cd(virtualPointsFolder);
gridTablesNames = {'grid n.', 'zVal', 'lambda_L', 'k_L',...
                   'nmseTIK' 'nccTIK' 'normcTIK' 'reTIK'...
                   'nmseTSVD' 'nccTSVD' 'normcTSVD' 'reTSVD'};

estimationCell = cell(nModes,1);

%% COMPUTATION LOOP 

% parameters setup
xAx = unique(hologramPoints(:,1));
yAx = unique(hologramPoints(:,2));
xStep = abs(xAx(end) - xAx(end-1));
yStep = abs(yAx(end) - yAx(end-1));

% The distance should be lattice/2  
% lattice = min(Dxhologram, Dyhologram) = min spacial sampling sampling step
zCenter = -0.5*min([xStep, yStep]);
plotData = true;
experimentalData = false;
 
userControl = input('choose z value // optimize z (!! iterations !!)[1,0] : ');
gridToUse = input('choose grid to use (integer positive): ');
virtualPtsFilename = ['VP_', int2str(gridToUse),'.csv'];

if userControl == 0
    % OPTIMIZATION
    
    maxIter = input('max n of iterations?: ');
    maxFun = input('max n of func evaluations=: ');
    VP_Params = zeros(3,nModes);
    zBound = 0.0075; % minimum value bound on z

    for mode = 1:nModes
        
        tStart = tic;
        % Setup of local variables
        omega = eigenFreqzRad(mode); % current eigenfreq mode
        disp(['f = ', num2str(omega/(2*pi))]);

        % pressure vector setup
        measuredPressure = pressureFields{mode};
        meshSize = numel(measuredPressure);
        measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed
        measuredPressure = whiteNoise(measuredPressure, 20); % add white gaussian noise to the mesurement

        % velocity GroundTruth vector setup
        v_GT = velocityFields{mode};   
        v_GT_vector = reshape( v_GT.', [numel(v_GT), 1]); 
        v_GT_vector(isnan(v_GT_vector)) = 0;
        v_GT_vector(isnan(violinInfos{4}(:,3)),:) = [];

        fun = @(x) applyESM( x, zBound, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                                X, Y , v_GT_vector, virtualPtsFilename,...
                                gridTablesNames, plotData, experimentalData );

        options = optimset('fminsearch');
        options = optimset(options, 'TolFun',1e-8,'TolX',1e-8, 'MaxFunEvals',1,'MaxIter', 1,...
        'DiffMinChange', 1, 'DiffMaxChange', 200); 

        % minimization
        [zpar,fval, exitflag, output] = fminsearch(fun, [zCenter 1 1].', options);
        zpar(1) = -(abs(zpar(1)) + abs(zBound));
        VP_Params(ii,:)= zpar;
    end
    
    %SAVE RESULTS 
    cd([estimationsFolder, '\controlParams'])
    filesList = ls(estimationsFolder);
    filesList(1:2,:) = [];
    numFile = length(filesList) +1;
    disp(['writing File VP_Params_', num2str(numFile),'.csv']);
    
    writeMat2File(VP_Params(1:3,:), ['VP_Params_', int2str(numFile),'.csv'], {'z [m]' 'scaleX' 'scaleY'} , 3, true)
    cd(baseFolder)
    
    elseif userControl == 1
    % CHOOSE z - LOOK PREVIOUS OPTIMUM
    msg = 'want to see previous optimized results ? [0,1] : ';
    paramsSet = input(msg);
    
    if paramsSet == 0
        msg = ['set distance virtual-violinMesh (theoric = ',num2str(zCenter),') ---> : '];
        zCenter = -abs(input(msg)); 
    else
       numFile = input('which params file to read [integer]: ');
       VP_Params = readmatrix(['VP_Params_',int2str(numFile),'.csv']);     
    end
    
    for ii = 1:nModes 
        tStart = tic;
        omega = eigenFreqzRad(ii); % current eigenfreq mode
        disp(['f = ', num2str(omega/(2*pi))]);

        % pressure vector setup
        measuredPressure = pressureFields{ii};
        meshSize = numel(measuredPressure);
        measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed
        measuredPressure = whiteNoise(measuredPressure, 20); % add white gaussian noise to the mesurement

        % velocity GroundTruth vector setup
        v_GT = velocityFields{ii};   
        v_GT_vector = reshape( v_GT.', [numel(v_GT), 1]); 
        v_GT_vector(isnan(v_GT_vector)) = 0;
        v_GT_vector(isnan(violinInfos{4}(:,3)),:) = [];
        
        if paramsSet == 1
            [LossFx, ESM_metrics , ESM_results ] = applyESM(VP_Params(ii,:), 0, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                               X, Y, v_GT_vector, virtualPtsFilename,...
                               gridTablesNames, plotData, experimentalData );   
        else
        [LossFx, ESM_metrics , ESM_results ] = applyESM([zCenter,1,1], 0, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                               X, Y, v_GT_vector, virtualPtsFilename,...
                               gridTablesNames, plotData, experimentalData );
        end
        
        estimationCell{ii} = ESM_metrics; 
        disp(toc(tStart)) 
    end
    
    cd([estimationsFolder,'\metrics'])
    estimationStruct = cell2struct(estimationCell, structNames, 1);
    
    filesList = ls(estimationsFolder);
    filesList(1:2,:) = [];
    numFile = length(filesList) +1;
    disp(['writing File estimationStruct', num2str(numFile),'.mat']);   
    save(['estimationStruct_', num2str(numFile),'.mat'], estimationStruct);
    cd(baseFolder)
end

%% Additional - See Virtual Points grids

cd(baseFolder)
filesList = ls(virtualPointsFolder);
filesList(1:2,:) = [];
cd(virtualPointsFolder)
figure(150)

for ii = 1:length(filesList)
  idx = find(filesList(ii,:)== '.');
  virtualPoints = table2array(readtable(filesList(ii,1:idx+3))); 
  plot3(virtualPoints(:,2), virtualPoints(:,1),  virtualPoints(:,3), '.', 'markerSize', 10);
  hold on;
  plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3));
  hold off;
  pause(1);
end

cd(baseFolder)
