clear all
close all
clc
%% NEAR FIELD ACOUSTIC HOLOGRAPHY - ESM METHOD %%

% folders
baseFolder = pwd;
virtualPointsFolder = [baseFolder,'\VPGrids'];
estimationsFolder = [baseFolder, '\Estimations\Experimental'];
addpath(genpath('functions'));
addpath(genpath('Data'));
addpath('violinMeshes');
addpath(virtualPointsFolder)
addpath(estimationsFolder)

%% global variables
plotHoloAndVel = false;

nMics = 8;
nMeas = 8;

nPressureMeas = nMics*nMeas;

rho = 1.2; % [Kg/m3] air density 

eigenFreqz = readtable('eigenfrequencies.csv');
eigenFreqz = table2array(eigenFreqz);
eigenFreqz = 2*pi*eigenFreqz; % convert in [rad/s]

nModes = length(eigenFreqz);

violinMesh = table2array(readtable('grid65x25.csv')); 

pX = length(unique(violinMesh(:, 1)));
pY = length(unique(violinMesh(:, 2)));
nViolinPoints = pX*pY;

%% See violin surface
X = reshape(violinMesh(:,1), [pY, pX]).';
Y = reshape(violinMesh(:,2), [pY, pX]).';
Z = reshape(violinMesh(:,3), [pY, pX]).';

if plotHoloAndVel
    figure(101)
    surf(X,Y,Z);
    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('z  [m]');
    zlim([0,0.1]); title('Violin Mesh');
end
%% Fetch and see hologram measured data
pressureData = readtable('pressure_Data.csv');
fNames = pressureData.Properties.VariableNames;
pressureData = table2array(pressureData);

hologramDistance = 0.03;
zHologram = max(Z(:)) + hologramDistance; 
hologramPoints =  [pressureData(:,1:2), zHologram*ones(size(pressureData(:,1)))] ; 
hologramMeshX = reshape( pressureData(:,1) , [nMeas, nMics]).'; 
hologramMeshY = reshape( pressureData(:,2) , [nMeas, nMics]).'; 
if plotHoloAndVel
    for ii = 1:(length(pressureData(1,:)) -2)
        figure(102)
        zPress = reshape(pressureData(:,2+ii), [nMeas, nMics]).';
        surf(hologramMeshX, hologramMeshY, abs(zPress)); view(2);
        xlabel('X   [m]');
        ylabel('Y   [m]');
        titleStr = strrep(strrep(fNames{2+ii},'__', ' '),'_','.');
        title([titleStr(1), titleStr(3:end), ' Hz']);
        hold off;
        pause(0.2);
    end
end


%% Fetch and see velocity groundtruth Data

load('xyDataVlnMeasurements.mat');
velocityFilename = 'velocity_Data';
velocityData = readtable([velocityFilename,'.csv']);
names = velocityData.Properties.VariableNames;
velocityData = table2array(velocityData);
if plotHoloAndVel
    for ii = 1:length(velocityData(1,3:end))
        [XX,YY,surfV] = getVelocityGroundtruth(velocityData(:,ii+2), velocityFilename, 103);
        sgtitle(names{ii+2}) ;
        pause(0.2);
    end
end

%% Calculate normal points wrt the surface for Green's fxs gradient

zNan = find(isnan(Z));
Z(zNan) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
Z(zNan) = nan;
normalPoints = [reshape(nx.', [nViolinPoints,1]),...
                reshape(ny.', [nViolinPoints,1]),...
                reshape(nz.', [nViolinPoints,1]) ];
normalPoints = sortrows(normalPoints);



%% STRUCT AND TABLE INIT
nEqSourceGrids = 1;
gridTablesNames = {'grid n.', 'zVal', 'lambda_L', 'k_L',...
                   'nmseTIK' 'nccTIK' 'normcTIK' 'reTIK'...
                   'nmseTSVD' 'nccTSVD' 'normcTSVD' 'reTSVD'};
for ii = 1:nModes
    structNames{ii} = ['f', int2str(ii)];
end 

estimationCell = cell(nModes,1);

%% COMPUTATION LOOP 

% parameters setup
xAx = unique(hologramPoints(:,1));
yAx = unique(hologramPoints(:,2));
xStep = abs(xAx - circshift(xAx,+1));
yStep = abs(yAx - circshift(yAx,+1));

% The distance should be lattice/2  
% lattice = min(Dxhologram, Dyhologram) = min spacial sampling sampling step
zCenter = -0.5*min([min(xStep), min(yStep)]);
experimentalData = true;

userControl = input('choose z value // optimize z (!! iterations !!)[1,0] : ');
fileList = cellstr(ls(virtualPointsFolder));
disp(fileList(3:end));

gridToUse = input('choose grid to use (integer positive): ');
virtualPtsFilename = ['VP_', int2str(gridToUse),'.csv'];
plotData = input('plot images of the algorithm?  [true/false]');

%% OPTIMIZATION
if userControl == 0
    maxIter = input('max n of iterations?: ');
    maxFun = input('max n of func evaluations=: ');
    VP_Params = zeros(3,nModes);
    zBound = 0.0075;  % minimum value bound on z

    for ii = 1:nModes 
        tStart = tic;
        omega = eigenFreqz(ii);
        disp(['f = ', num2str(omega/(2*pi))]);
        measuredPressure = pressureData(:,ii+2);
        v_GT_vector = velocityData(:,ii+2);
        fun = @(x) applyESM(x, zBound, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                                xData, yData , v_GT_vector, virtualPtsFilename,...
                                gridTablesNames, plotData, experimentalData );

        options = optimset('fminsearch');
        options = optimset(options, 'TolFun',1e-8,'TolX',1e-8, 'MaxFunEvals',maxFun,'MaxIter', maxIter,...
        'DiffMinChange', 1, 'DiffMaxChange', 200); 

        % minimization
        [zpar,fval, exitflag, output] = fminsearch(fun, [zCenter, 1], options);
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
end

%% ESM CHOOSING Z  or LOOK AT PREVIOUS OPTIMUM
if userControl == 1
    
    msg = 'want to see previous optimized results ? [0,1] : ';
    paramsSet = input(msg);
    
    if paramsSet == 0
        msg = ['set distance virtual-violinMesh (theoric = ',num2str(zCenter),') ---> : '];
        zCenter = input(msg); 
    else
       numFile = input('which params file to read [integer]: ');
       VP_Params = readmatrix(['VP_Params_',int2str(numFile),'.csv']);     
    end
    
    for ii = 1:nModes 
        tStart = tic;
        omega = eigenFreqz(ii);
        disp(['f = ', num2str(omega/(2*pi))]);
        measuredPressure = pressureData(:,ii+2);
        v_GT_vector = velocityData(:,ii+2);
        
        if paramsSet == 1
            [LossFx, ESM_metrics , ESM_results ] = applyESM(VP_Params(ii,:), 0, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                               xData, yData , v_GT_vector, virtualPtsFilename,...
                               gridTablesNames, plotData, experimentalData );   
        else
        [LossFx, ESM_metrics , ESM_results ] = applyESM([zCenter,1,1], 0, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                               xData, yData , v_GT_vector, virtualPtsFilename,...
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
    save(['estimationStruct_', num2str(numFile),'.mat'], 'estimationStruct');
    cd(baseFolder)
    
    %% plot of metrics
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
xticklabels(eigenFreqz./(2*pi))
xtickangle(90)
ylim([0.4 1])
grid on
xlabel('N modes');
ylabel('NCC');
title('NCC - experimental')
legend('Tikhonov', 'TSVD')

tikERPlot = zeros(11,1);
tsvdERPlot = zeros(11,1);
for ii = 1:11
    tikERPlot(ii) = table2array(estimationCell{ii}(1,7));
    tsvdERPlot(ii) = table2array(estimationCell{ii}(1,11));
end
figure(222)
plot(tikERPlot, 'b-o')
hold on
plot(tsvdERPlot, 'r-o')
xticks(1:11)
xticklabels(eigenFreqz./(2*pi))
xtickangle(90)

grid on
xlabel('N modes')
ylabel(' RE %');
title('RE % - experimental')
legend('Tikhonov', 'TSVD')

end

