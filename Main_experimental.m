
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
addpath(genpath('Data'));
addpath('violinMeshes');
addpath(virtualPointsFolder)

%% global variables

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

% get normal points for Green's fxs gradient
X = reshape(violinMesh(:,1), [pY, pX]).';
Y = reshape(violinMesh(:,2), [pY, pX]).';
Z = reshape(violinMesh(:,3), [pY, pX]).';
zNan = find(isnan(Z));
Z(zNan) = 0; % for boundaries normal vector

[nx , ny, nz] = surfnorm(X,Y,Z); % returns the x, y, and z components of the three-dimensional surface normals 
                                 % for each point of the surface.
                                 % surfnorm(X',Y',Z') to invert the vector
                                 % direction
Z(zNan) = nan;
normalPoints = [reshape(nx, [nViolinPoints,1]),...
                reshape(ny, [nViolinPoints,1]),...
                reshape(nz, [nViolinPoints,1]) ];
normalPoints = sortrows(normalPoints);

pressureData = table2array(readtable('pressure_Data.csv'));

hologramDistance = 0.02;
zHologram = max(Z(:)) + hologramDistance; 
hologramPoints =  [pressureData(:,1:2), zHologram*ones(size(pressureData(:,1)))] ; 
hologramMeshX = reshape( pressureData(:,1) , [nMeas, nMics]).'; 
hologramMeshY = reshape( pressureData(:,2) , [nMeas, nMics]).'; 

for ii = 1:(length(pressureData(1,:)) -2)
    figure(300)
    zPress = reshape(pressureData(:,2+ii), [nMeas, nMics]);
    surf(hologramMeshX, hologramMeshY, abs(zPress)); view(2);
    xlabel('X   [m]');
    ylabel('Y   [m]');
    title(['f_{', int2str(ii), '}']);
end

velocityData = table2array(readtable('velocityData.csv'));   
load('xyDataVlnMeasurements.mat');

XX = xData.';
YY = yData.';
% for ii = 1:10
%    if ii <= 7
%    XX(ii, isnan(XX(ii,:))) = mean(XX(ii, ~isnan(XX(ii,:))));
%    end
%    YY( isnan(YY(:,ii)), ii) = mean(YY( ~isnan(YY(:, ii)), ii));
% end
% 
% for ii = 1:length(velocityData(1,:)) -2
%     figure(301)
%     F1 = scatteredInterpolant(velocityData(:,1),velocityData(:,2), abs(velocityData(:,2+ii)));
%     figure
%     vq1 = F1(X,Y);
%     
%     subplot 121
%     surf(X,Y,vq1)
%     xlabel('X   [m]');
%     ylabel('Y   [m]');
%     title(['f_{', int2str(ii), '}']);
%     
%     vRows = length(vq1(:,1)); vCols = length(vq1(1,:));
%     zVel = reshape(vq1, [ vCols*vRows, 1 ]);
%     zVel(isnan(reshape(Z,[pX*pY,1]))) = nan;
%     
%     vq = reshape(zVel, [ vRows, vCols]);
%     subplot 122
%     
%     surf(X,Y,vq)
%     xlabel('X   [m]');
%     ylabel('Y   [m]');
%     title(['f_{', int2str(ii), '}']);
%     
% end


%% STRUCT AND TABLE INIT
nEqSourceGrids = 1;
nZpoints = 1;
gridTablesNames = {'grid n.', 'zVal', 'lambda_L', 'k_L', 'nmseTSVD_L', 'nccTSVD_L',...
                    'nmseTIK_L','nccTIK_L', 'lambda_nmse_M', 'k_nmse_M', ...
                    'k_ncc_M', 'lambda_ncc_M', 'nmseTSVD_M', 'nccTSVD_M', ...
                    'nmseTIK_M', 'nccTIK_M', };
structNames = {'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11'};
               
dataCell = cell(length(eigenFreqz),1);
ZreguFreq = cell(length(eigenFreqz),1);


%% COMPUTATION LOOP 
bestZ = table2array(readtable('bestZ.csv'));
 for ii = 1:length(eigenFreqz) 
  tStart = tic;
  omega = eigenFreqz(ii);
  disp(['f = ', num2str(omega/(2*pi))]);
  measuredPressure = pressureData(:,ii);
  v_ex_vector = velocityData(:,  ii);
  
    xAx = unique(hologramPoints(:,1));
    yAx = unique(hologramPoints(:,2));
    xStep = abs(xAx - circshift(xAx,+1));
    yStep = abs(yAx - circshift(yAx,+1));
    
  zCenter = -0.5*min([min(xStep), min(yStep)]);
  zSearch = 0;
  transposeGrids = false;
  plotData = true;
  experimentalData = true;

  [reguData, ZreguDatas] = getBestGrid(nEqSourceGrids, measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
                           nZpoints, zCenter, zSearch, xData, yData , v_ex_vector,...
                           rho, pX, pY,gridTablesNames, transposeGrids, plotData, experimentalData);
  
%   fun = @(x) getBestGridSimple( measuredPressure, hologramPoints, normalPoints, violinMesh , omega,...
%                            x, xData, yData , v_ex_vector,...
%                            rho, pX, pY,gridTablesNames(2:8), transposeGrids, plotData, experimentalData );
%                        
%   options = optimset('fminsearch');
%   options = optimset(options, 'TolFun',1e-8,'TolX',1e-8, 'MaxFunEvals',2e2,'MaxIter', 5e3,...
%     'DiffMinChange', 1, 'DiffMaxChange', 200); 

%   % minimization
%   [zpar,fval, exitflag, output] = fminsearch(fun, zCenter, options);
%                        
  tempTable = array2table( reguData , 'VariableNames',gridTablesNames);
       
    dataCell{ii} = tempTable;
    dataStruct = cell2struct(dataCell, structNames, 1);
    ZreguFreq{ii} = ZreguDatas;
    disp(toc(tStart))
%   bestZ(2,ii)= zpar;
 end
%  writeMat2File(bestZ, ['bestZ.csv'], {'f'} , 1, false);
 
 %% plot metrics
 
 % plot figures for M method
 nccTik_M_max = zeros(1,length(eigenFreqz));
 nccTik_M_max(1) = max(dataStruct.f1.nccTIK_M);
 nccTik_M_max(2) = max(dataStruct.f2.nccTIK_M);
 nccTik_M_max(3) = max(dataStruct.f3.nccTIK_M);
 nccTik_M_max(4) = max(dataStruct.f4.nccTIK_M);
 nccTik_M_max(5) = max(dataStruct.f5.nccTIK_M);
 nccTik_M_max(6) = max(dataStruct.f6.nccTIK_M);
 nccTik_M_max(7) = max(dataStruct.f7.nccTIK_M);
 nccTik_M_max(8) = max(dataStruct.f8.nccTIK_M);
 nccTik_M_max(9) = max(dataStruct.f9.nccTIK_M);
 nccTik_M_max(10) = max(dataStruct.f10.nccTIK_M);
 nccTik_M_max(11) = max(dataStruct.f11.nccTIK_M);
 
  nccTsvd_M_max = zeros(1,length(eigenFreqz));
 nccTsvd_M_max(1) = max(dataStruct.f1.nccTSVD_M);
 nccTsvd_M_max(2) = max(dataStruct.f2.nccTSVD_M);
 nccTsvd_M_max(3) = max(dataStruct.f3.nccTSVD_M);
 nccTsvd_M_max(4) = max(dataStruct.f4.nccTSVD_M);
 nccTsvd_M_max(5) = max(dataStruct.f5.nccTSVD_M);
 nccTsvd_M_max(6) = max(dataStruct.f6.nccTSVD_M);
 nccTsvd_M_max(7) = max(dataStruct.f7.nccTSVD_M);
 nccTsvd_M_max(8) = max(dataStruct.f8.nccTSVD_M);
 nccTsvd_M_max(9) = max(dataStruct.f9.nccTSVD_M);
 nccTsvd_M_max(10) = max(dataStruct.f10.nccTSVD_M);
 nccTsvd_M_max(11) = max(dataStruct.f11.nccTSVD_M);
 
 % plot a figure with the min NCCC for the various frequencies
 nmseTik_M_min = zeros(1,length(eigenFreqz));
 nmseTik_M_min(1) = min(dataStruct.f1.nmseTIK_M);
 nmseTik_M_min(2) = min(dataStruct.f2.nmseTIK_M);
 nmseTik_M_min(3) = min(dataStruct.f3.nmseTIK_M);
 nmseTik_M_min(4) = min(dataStruct.f4.nmseTIK_M);
 nmseTik_M_min(5) = min(dataStruct.f5.nmseTIK_M);
 nmseTik_M_min(6) = min(dataStruct.f6.nmseTIK_M);
 nmseTik_M_min(7) = min(dataStruct.f7.nmseTIK_M);
 nmseTik_M_min(8) = min(dataStruct.f8.nmseTIK_M);
 nmseTik_M_min(9) = min(dataStruct.f9.nmseTIK_M);
 nmseTik_M_min(10) = min(dataStruct.f10.nmseTIK_M);
 nmseTik_M_min(11) = min(dataStruct.f11.nmseTIK_M);
 
 nmseTsvd_M_min = zeros(1,length(eigenFreqz));
 nmseTsvd_M_min(1) = min(dataStruct.f1.nmseTSVD_M);
 nmseTsvd_M_min(2) = min(dataStruct.f2.nmseTSVD_M);
 nmseTsvd_M_min(3) = min(dataStruct.f3.nmseTSVD_M);
 nmseTsvd_M_min(4) = min(dataStruct.f4.nmseTSVD_M);
 nmseTsvd_M_min(5) = min(dataStruct.f5.nmseTSVD_M);
 nmseTsvd_M_min(6) = min(dataStruct.f6.nmseTSVD_M);
 nmseTsvd_M_min(7) = min(dataStruct.f7.nmseTSVD_M);
 nmseTsvd_M_min(8) = min(dataStruct.f8.nmseTSVD_M);
 nmseTsvd_M_min(9) = min(dataStruct.f9.nmseTSVD_M);
 nmseTsvd_M_min(10) = min(dataStruct.f10.nmseTSVD_M);
 nmseTsvd_M_min(11) = min(dataStruct.f11.nmseTSVD_M);

figure(111)
subplot 211
plot(1:11, nmseTik_M_min, '-o')
hold on
plot(1:11, nmseTsvd_M_min, '-o')
hold off
grid on
xticks([1:10])
xticklabels({'111 ','186','285', '375', '498', '684', '730', '841', '984', '1196', '1374'})
legend('TIK','TSVD')
title('min NMSE method M')
xlabel('f [Hz]')

subplot 212
plot(1:11, nccTik_M_max, '-o')
hold on
plot(1:11, nccTsvd_M_max, '-o')
hold off
grid on
xticks([1:10])
xticklabels({'111 ','186','285', '375', '498', '684', '730', '841', '984', '1196', '1374'})
legend('TIK','TSVD')
title('max NCC method M')
xlabel('f [Hz]')

 % check regularization parameters as frequency varies
 
 %% plot the best result as example
 % according to the M method
 
 bestLambda = 31;
 bestMfreq = eigenFreqz(2);
 bestMgrid = 10;
 measuredPressure = pressureData(:, 2 + 2);
 v_ex_vector = velocityData(:, 2 + 2);
 
 virtualPtsFilename = ['VP_', int2str(bestMgrid), '.csv'];
 virtualPoints = table2array(readtable(virtualPtsFilename)) ;
 virtualPoints = 0.001.*virtualPoints;
 [G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , omega );
 G_p_omega = G_p{1};
 [G_v] = normalGradient(virtualPoints, violinMesh , omega, normalPoints);
 G_v_omega = G_v{1};
 [U,s,V] = csvd (G_p_omega);
 q_TIK =  (1/(1i*omega*rho)).* tikhonov(U,s,V,measuredPressure,bestLambda);
 v_TIK = G_v_omega*q_TIK;
 v_TIK_Fin = addNans(violinMesh, v_TIK);
 surfVelRecTIK = reshape( v_TIK_Fin , [pX, pY]); 
 
 figure(500)
 surf(X, Y, abs(surfVelRecTIK))