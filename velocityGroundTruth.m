close all
clear all
clc
%% Ground truth computation

addpath('functions')
addpath('Data\csvData')
% compute FRF of the velocity


baseFolder = pwd;
accFolder = [baseFolder, '\acq_11_06'];
matDataFolder = [baseFolder, '\Data\matData'];

Fs = 48000; % sampling frequency
duration = 2; % [s]
signalLength = duration*Fs;
numberPoints = 27; % number of measurement points in the violin plate
numberAcquisitions = 6; % number of measurements for each point

estimatorMatrix = zeros(3989,numberPoints); % store the H1 estimator
cleanEstimatorMatrix = zeros(3989,numberPoints); % store the H1 estimator denoised
MobilityMatrix = zeros(3989,numberPoints); % store the mobility transfer function (fft(Y)/fft(X))
cleanMobilityMatrix = zeros(3989,numberPoints);% store the Mobility denoised

% import raw data 

measurementPts = [-3 4; -2 4; 0 4; 2 4; 3 4; %coordinates of the measurement points
                  -2 2;  0 2; 2 2;
                  -1 1;  1 1;
                  -2 0;  0 0; 2 0;
                  -1 -1; 1 -1;
                  -2 -2; 0 -2; 2 -2;
                  -2 -4; 0 -4; 2 -4;
                  -3 -5; 3 -5;
                  -2 -6; 0 -6; 2 -6;
                   0 -8];

% temportary file to compute the frequency responses
forceTemp = zeros(signalLength, numberAcquisitions);
accelerationTemp = zeros(signalLength, numberAcquisitions);
FRFTemp = zeros(signalLength/2, numberAcquisitions);


for jj = 1:numberPoints
    
    cd(accFolder)
    fileName = ['x=', int2str(measurementPts(jj,1)),' y=',int2str(measurementPts(jj,2)),'_acq_'];
    
    for ii = 1:numberAcquisitions
        
        tempFileName = [fileName, int2str(ii-1),'.txt'];
        rawTemp = importdata([accFolder,'\', tempFileName]);
        rawTemp = strrep(rawTemp, ',', '.');
        rawMes = zeros(length(rawTemp), 3);
        
        for k = 1:length(rawTemp) % convert txt file into arrays
            rawMes(k,:) = str2num(rawTemp{k});
        end
        
        forceTemp(:,ii) = rawMes(:,1);
        accelerationTemp(:,ii) = rawMes(:,2);
        
        
        % compute the FFT for each acquisition
            Y1 = fft(accelerationTemp(:,ii));
            X = fft(forceTemp(:,ii));           
            Y1 = abs(Y1/signalLength); Y1_cut = Y1(1:signalLength/2); Y1_cut(2:end-1) = 2*Y1_cut(2:end-1);
            X = abs(X/signalLength); X_cut = X(1:signalLength/2); X_cut(2:end-1) = 2*X_cut(2:end-1);
            
            freq = Fs*(0:(signalLength/2)-1)/signalLength;          
            freq(1) = 1;
            freq = freq';
            FRFTemp(:,ii) =  (1./(1i*2*pi*freq)).*(Y1_cut./X_cut);
        
    end
    
    [pxy1, f] = cpsd(accelerationTemp, forceTemp,[],[],signalLength, Fs);
    [pxx, f] = cpsd(forceTemp,forceTemp,[],[],signalLength, Fs);
    
        
    cutIdxs = find(f <2000 & f>5 );
    f = f(cutIdxs); pxx = pxx(cutIdxs,:);pxy1 = pxy1(cutIdxs,:); 
    
    pxy1 = 1./(1i*2*pi*f).*pxy1;
    
    H1 = sum(pxy1,2)./sum(pxx,2);
    
    estimatorMatrix(:,jj) = H1;
    
    MobilityMatrix(:,jj) = mean(FRFTemp(cutIdxs,:),2);
    
end

figure(808)
subplot 211
semilogy(f, (abs(estimatorMatrix)))
title('H1 without SVD')
subplot 212
semilogy(f, (abs(MobilityMatrix)))
title('Y/X')

%% apply SVD to H1 estimator and Mobility function

% SVD Parameters
M = 10;
thresholdIdx = 1;

for ii = 1:numberPoints
    H1Temp = estimatorMatrix(:,ii);
    [H1_clean,singularVals] = SVD(H1Temp.', f, M, thresholdIdx, false);
    cleanEstimatorMatrix(:,ii) = H1_clean;
end

figure(809)
semilogy(f, (abs(cleanEstimatorMatrix)))
title('H1 with SVD')

for ii = 1:numberPoints
    H1Temp = MobilityMatrix(:,ii);
    [H1_clean,singularVals] = SVD(H1Temp.', f, M, thresholdIdx, false);
    cleanMobilityMatrix(:,ii) = H1_clean;
end

figure(890)
semilogy(f, (abs(cleanMobilityMatrix)))
title('Mobility with SVD')

cd(matDataFolder)

save('velocityH1.mat','estimatorMatrix');
save('velocityH1cleaned.mat','cleanEstimatorMatrix');
save('velocityMobility.mat','MobilityMatrix');
save('velocityMobilitycleaned.mat','cleanMobilityMatrix');

cd(baseFolder)
%% find resonance frequencies from velocity

[peakPositions] = peaks(cleanMobilityMatrix, f);
fpeakPositions = f(peakPositions);
idxPks = find(fpeakPositions > 1400);
fpeakPositions(idxPks) = [];

figure(899)
semilogy(f, (abs(cleanMobilityMatrix)))
hold on
for i = 1:length(fpeakPositions)
    
    xline(fpeakPositions(i))
    
end
hold off
title('Mobility with SVD')


%% Fill the velcoity matrix from H1 with pressure peaks

eigenFreqPress = readmatrix('eigenFreqPress.csv');
eigenFreqPressIdx = eigenFreqPress*2;

velocitiesMatrix = zeros(numberPoints, length(eigenFreqPressIdx));

for ii = 1:numberPoints
    for jj = 1:length(eigenFreqPress)
        velocitiesMatrix(ii,jj)= abs(cleanEstimatorMatrix(eigenFreqPressIdx(jj),ii));
    end
end

%% store velocities in the grid points

xyCoord = readmatrix('velocityPoints.csv');
xyCoord(:,3) = [];
xyCoord = xyCoord*10;

for ii = 1:length(eigenFreqPress)
    xyCoord = [xyCoord velocitiesMatrix(:,ii)];
end

xyCoordSorted = sortrows(xyCoord);

label ={'x' 'y'};
    freqLabel = round(eigenFreqPress);

for ii = 1: length(eigenFreqPress)
    label{ii+2} = ['f',num2str(ii),' = ', num2str(freqLabel(ii))];
end


%% export the velocity grid

% auto approach 
geomData = readtable('geometryVelocity.xlsx');

centerData = table2array(geomData(1:5,1:3));
XX = table2array(geomData(12:21,7:13));
YY = table2array(geomData(12:21,16:22));

xData = table2array(geomData(1:10,7:13)).*sign(XX)*10;
yData = table2array(geomData(1:10,16:22)).*sign(YY)*10;
[iMask, jMask] = find(isnan(xData));

zData =nan*ones(length(xData(:,1)),length(xData(1,:)), length(eigenFreqPress));
orderedPoints = zeros(numberPoints, 3);
orderedVelMatrices = zeros(size(velocitiesMatrix));

for ii = 1: numberPoints
    [xx, yy] = find(XX == measurementPts(ii,1) & YY == measurementPts(ii,2));
    orderedPoints(ii,1) = xData(xx,yy);
    orderedPoints(ii,2) = yData(xx,yy);
    disp(orderedPoints(ii,:));
    disp([xx,yy]);
    zData(xx,yy,:) = velocitiesMatrix(ii,:);
%   orderedPoints(ii,3) = velocities(ii); 
%   zData(xx,yy) = velocities(ii);
end  
save('zDataVInMeasurements', 'zData');

writeMat2File(orderedPoints, 'velocityPoints.csv', {'x' 'y' 'z'}, 3, true );
% x = reshape(xData, numberPoints,1);
% y = reshape(yData, numberPoints,1);
% z = reshape(zData, numberPoints, 1);
% 
% figure()
% subplot 121
% plot3(x, y, z, '.');
% subplot 122
% surf(xData, yData, zData);

