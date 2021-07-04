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

plotImportImgs = false;

% import synthetic data
[violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = ...
    importData(velocityFileName, pressureFileName, nMics, nMeas, plotImportImgs);
fCut = 800;
eigenFreqz = eigenFreqz(eigenFreqz<fCut);
eigenFreqzRad = 2*pi*eigenFreqz; % convert in [rad/s]

%conversion from [mm] to [m]
for ii =1:4
        hologramInfos{ii} = hologramInfos{ii}*0.001;
        violinInfos{ii} = violinInfos{ii}*0.001;
end

%% Setup of global variables 

nModes = length(eigenFreqz);
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

for ii = 1:nModes
    structNames{ii} = ['f', int2str(ii)];
end 

%% COMPUTATION LOOP 

% parameters setup
xAx = unique(hologramPoints(:,1));
yAx = unique(hologramPoints(:,2));
xStep = abs(xAx(end) - xAx(end-1));
yStep = abs(yAx(end) - yAx(end-1));

% The distance should be lattice/2  
% lattice = min(Dxhologram, Dyhologram) = min spacial sampling sampling step
zCenter = -0.5*min([xStep, yStep]);

experimentalData = false;
 
userControl = input('choose z value // optimize z (!! iterations !!)[1,0] : ');
fileList = cellstr(ls(virtualPointsFolder));
disp(fileList(3:end));
gridToUse = input('choose grid to use (integer positive): ');
virtualPtsFilename = ['VP_', int2str(gridToUse),'.csv'];
plotData = input('plot images of the algorithm?  [true/false]');
addNoise = input('Noise on the pressure? [true/false]: ');
if addNoise
SNR = input('Specify Signal To Noise Ratio (SNR): ');
end
%% OPTIMIZATION

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
        if addNoise
        measuredPressure = whiteNoise(measuredPressure, SNR); % add white gaussian noise to the mesurement
        end
        % velocity GroundTruth vector setup
        v_GT = velocityFields{mode};   
        v_GT_vector = reshape( v_GT.', [numel(v_GT), 1]); 
        v_GT_vector(isnan(v_GT_vector)) = 0;
        v_GT_vector(isnan(violinInfos{4}(:,3)),:) = [];

        fun = @(x) applyESM( x, zBound, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                                X, Y , v_GT_vector, virtualPtsFilename,...
                                gridTablesNames, plotData, experimentalData );

        options = optimset('fminsearch');
        options = optimset(options, 'TolFun',1e-4,'TolX',1e-6, 'MaxFunEvals',maxFun,'MaxIter', maxIter); 

        % minimization
        [zpar,fval, exitflag, output] = fminsearch(fun, [zCenter 1 1].', options);
        zpar(1) = -(abs(zpar(1)) + abs(zBound));
        VP_Params(:,mode)= zpar;
     end
    
    %SAVE RESULTS 
    cd([estimationsFolder, '\controlParams'])
    filesList = ls(estimationsFolder);
    filesList(1:2,:) = [];
    numFile = length(filesList(:,1)) +1;
    disp(['writing File VP_Params_', num2str(numFile),'.csv']);
    writeMat2File(VP_Params.', ['VP_Params_', int2str(numFile),'.csv'], {'z [m]' 'scaleX' 'scaleY'} , 3, true)
    cd(baseFolder)
end

%% ESM CHOOSING Z  or LOOK AT PREVIOUS OPTIMUM
if userControl == 1
    % CHOOSE z - LOOK PREVIOUS OPTIMUM
    msg = 'want to see previous optimized results ? [0,1] : ';
    paramsSet = input(msg);
    
    if paramsSet == 0
        msg = ['set distance virtual-violinMesh (theoric = ',num2str(zCenter),') ---> : '];
        zCenter = -abs(input(msg)); 
    else
       cd([estimationsFolder,'\controlParams']);
       fileList = ls([estimationsFolder,'\controlParams']);
       disp(fileList(3:end,:));
       numFile = input('which params file to read [integer]: ');
       VP_Params = readmatrix(['VP_Params_',int2str(numFile),'.csv']);     
       cd(baseFolder);
    end
%%
    for ii = 1:nModes 
        tStart = tic;
        omega = eigenFreqzRad(ii); % current eigenfreq mode
        disp(['f = ', num2str(omega/(2*pi))]);

        % pressure vector setup
        measuredPressure = pressureFields{ii};
        meshSize = numel(measuredPressure);
        measuredPressure = reshape(measuredPressure , [meshSize,1]); % convert the measurement matrix into an array... the magnitude of pressure is needed
        if addNoise
            measuredPressure = whiteNoise(measuredPressure,SNR); % add white gaussian noise to the mesurement
        end
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
    %%
    cd([estimationsFolder,'\metrics'])
    estimationStruct = cell2struct(estimationCell, structNames, 1);
    
    filesList = ls(estimationsFolder);
    filesList(1:2,:) = [];
    numFile = length(filesList(:,1)) +1;
    disp(['writing File estimationStruct', num2str(numFile),'.mat']);   
    save(['estimationStruct_', num2str(numFile),'.mat'], 'estimationStruct');
    cd(baseFolder)
    
    
% plot of metrics
tikNCCsPlot = zeros(11,1);
tsvdNCCsPlot = zeros(11,1);
for ii = 1:11
    tikNCCsPlot(ii) = table2array(estimationCell{ii}(1,5));
    tsvdNCCsPlot(ii) = table2array(estimationCell{ii}(1,9));
end
figure(221)
plot(tikNCCsPlot, 'b-o')
hold on
plot(tsvdNCCsPlot, 'r-o')
xticks(1:11)
xticklabels(eigenFreqz)
xtickangle(90)
ylim([0.4 1])
grid on
xlabel('f  [Hz]');
ylabel('NCC');
title('NCC - synthetic')
legend('Tikhonov', 'TSVD')

tikNMSEPlot = zeros(11,1);
tsvdNMSEPlot = zeros(11,1);
for ii = 1:11
    tikNMSEPlot(ii) = table2array(estimationCell{ii}(1,4));
    tsvdNMSEPlot(ii) = table2array(estimationCell{ii}(1,8));
end
figure(222)
plot(tikNMSEPlot, 'b-o')
hold on
plot(tsvdNMSEPlot, 'r-o')
xticks(1:11)
xticklabels(eigenFreqz)
xtickangle(90)

grid on
xlabel('f  [Hz]');
ylabel('NMSE  [dB]')
title('NMSE - synthetic')
legend('Tikhonov', 'TSVD')
end

%% Additional - See Virtual Points grids

cd(baseFolder)
filesList = ls(virtualPointsFolder);
filesList(1:2,:) = [];
cd(virtualPointsFolder)
figure(150)
if length(filesList(:,1)) == 1
    filesAxis = 1;
else
    filesAxis = 1:length(filesList(:,1));
end
for ii = filesAxis
  idx = find(filesList(ii,:)== '.');
  virtualPoints = table2array(readtable(filesList(ii,1:idx+3))); 
  plot3(virtualPoints(:,2), virtualPoints(:,1),  virtualPoints(:,3), '.', 'markerSize', 10);
  hold on;
  plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3),'x');
  hold off;
  pause(1);
end

cd(baseFolder)


