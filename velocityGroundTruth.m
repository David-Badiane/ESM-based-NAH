close all
clear all
clc
%% Ground truth computation

addpath('functions')
addpath(genpath('Data'))
addpath('violinMeshes');
baseFolder = pwd;
addpath(genpath('Exp_Measurements'))

% compute FRF of the velocity


baseFolder = pwd;
matDataFolder = [baseFolder, '\Data\matData\Velocity'];
csvDataFolder = [baseFolder, '\Data\csvData'];

expFolder = [baseFolder, '\Exp_Measurements'];
accFolder = [expFolder, '\ACC_acq_11_06'];

velocityFilename = 'velocity_Data.csv';

%% READ DATA
Fs = 48000; % sampling frequency
duration = 2; % [s]
t = [0:1/Fs:duration - 1/Fs].';
signalLength = duration*Fs;
numberPoints = 27; % number of measurement points in the violin plate
numberAcquisitions = 6; % number of measurements for each point

%coordinates of the measurement points
measurementPts = [-3 4; -2 4; 0 4; 2 4; 3 4; 
                  -2 2;  0 2; 2 2;
                  -1 1;  1 1;
                  -2 0;  0 0; 2 0;
                  -1 -1; 1 -1;
                  -2 -2; 0 -2; 2 -2;
                  -2 -4; 0 -4; 2 -4;
                  -3 -5; 3 -5;
                  -2 -6; 0 -6; 2 -6;
                   0 -8];

% preallocation
estimatorMatrix = []; % store the H1 estimator
mobilityMatrix = []; % store the mobility transfer function (fft(Y)/fft(X))        
  
highCut = 2000;
lowCut = 5;
alpha = 10; % exponential filter decay constant
offsetFilter = 0.1 * Fs; % 100 ms of offset

% temportary file to compute the frequency responses
forceTemp = zeros(signalLength, numberAcquisitions);
accelerationTemp = zeros(signalLength, numberAcquisitions);
FRFTemp = zeros(signalLength/2, numberAcquisitions);

for jj = 1:numberPoints
    cd(accFolder)
    fileName = ['x=', int2str(measurementPts(jj,1)),' y=',int2str(measurementPts(jj,2)),'_acq_'];
    for ii = 1:numberAcquisitions
        tempAcqFilename = [fileName, int2str(ii-1),'.csv'];
        expFilter = circshift(exp(-alpha*t),offsetFilter);
        expFilter(1:offsetFilter) = ones(1,offsetFilter);
        
        rawAccMeasurements = readmatrix(tempAcqFilename);
        forceTemp(:,ii) = expFilter.*rawAccMeasurements(:,1);
        accelerationTemp(:,ii) = expFilter.*rawAccMeasurements(:,2);

        % compute the FFT for each acquisition
        Y1 = FFT(accelerationTemp(:,ii));
        X = FFT(forceTemp(:,ii));           
        freq = (Fs*(1:(signalLength/2))/signalLength).';      
        FRFTemp(:,ii) =  (1./(1i*2*pi*freq)).*(Y1./X);        
    end
    
    [pxy1, f] = cpsd(accelerationTemp, forceTemp,[],[],signalLength, Fs);
    [pxx , f] = cpsd(forceTemp,forceTemp,[],[],signalLength, Fs);
   
    cutIdxs = find(f <highCut & f>lowCut );
    
    fAxis = f(cutIdxs); pxx = pxx(cutIdxs,:);pxy1 = pxy1(cutIdxs,:); 
        
    H1 =  1./(1i*2*pi*fAxis).*sum(pxy1,2)./sum(pxx,2);
%     figure(114)
%     semilogy(f, abs(H1));
%     xlim([0,2000]);

    estimatorMatrix = [estimatorMatrix H1];
    mobilityMatrix =[ mobilityMatrix   mean(FRFTemp(cutIdxs,:),2)]; 
end

figure(808)
subplot 211
semilogy(fAxis, (abs(estimatorMatrix)))
title('H1 without SVD')
subplot 212
semilogy(fAxis, (abs(mobilityMatrix)))
title('Y/X')


%% apply SVD to H1 estimator and Mobility function

% SVD Parameters
singValsNum = 7;
usedSingVals = 1;

cleanEstimatorMatrix = []; % store the H1 estimator denoised
cleanMobilityMatrix = [];% store the Mobility denoise       

% to choose correct parameters for SVD
[H1_clean,singularVals] = SVD(H1.', fAxis, singValsNum, usedSingVals, true);

for ii = 1:numberPoints
    H1Temp = estimatorMatrix(:,ii);
    [H1_clean,singularVals] = SVD(H1Temp.', fAxis, singValsNum, usedSingVals, false);
    cleanEstimatorMatrix = [cleanEstimatorMatrix H1_clean.'];
    
    mobilityTemp = mobilityMatrix(:,ii);
    [mobility_clean,singularVals] = SVD(mobilityTemp.', fAxis,singValsNum, usedSingVals, false);
    cleanMobilityMatrix = [ cleanMobilityMatrix mobility_clean.'];
end

figure(809)
semilogy(fAxis, (abs(cleanEstimatorMatrix)))
title('H1 with SVD')

figure(890)
semilogy(fAxis, (abs(cleanMobilityMatrix)))
title('Mobility with SVD')

cd(matDataFolder)
save('velocityH1.mat','estimatorMatrix');
save('velocityH1cleaned.mat','cleanEstimatorMatrix');
save('velocityMobility.mat','mobilityMatrix');
save('velocityMobilitycleaned.mat','cleanMobilityMatrix');
cd(baseFolder)

%% find resonance frequencies from velocity
% min peak prominence and min peak width B - high setting
highPeaksParams = [120,10];
% min peak prominence and min peak width A - low setting
lowPeaksParams = [100,5];
% fCut low for peaks to ignore
ignorePeaksLow = 100;
% fCut high for peaks to ignore
ignorePeaksHigh = 1400;
% frequency threshold for peak finders A - B configuration
fThreshold = 250;

[peakPositions, fPeaks] = peaks(cleanEstimatorMatrix, fAxis, fThreshold,...
    ignorePeaksLow, ignorePeaksHigh, highPeaksParams, lowPeaksParams);


figure(899)
semilogy(fAxis, (abs(cleanMobilityMatrix)))
hold on  
xline(fPeaks) 
hold off
title('Mobility with SVD')


%% Fill the velocity matrix 
nPeaks = length(peakPositions);
velocitiesMatrix = zeros(numberPoints, nPeaks);

for ii = 1:numberPoints
    for jj = 1:nPeaks
        velocitiesMatrix(ii,jj)= cleanEstimatorMatrix(peakPositions(jj),ii);
    end
end

velocitiesMatrix2 = zeros(numberPoints, nPeaks);
for ii = 1:numberPoints
        [Hv,f0, fLocs, csis, Q, shapes] = EMASimple(cleanEstimatorMatrix(:,ii), fAxis, 1e-4, 8, false);
        
        map = [];
        tracked = [];
        for jj = 1:nPeaks
            diff = abs(f0 - fPeaks(jj));
            [minVal, minLoc] = min(diff);
            map = [map minLoc];
            
            if minVal >= 0.05*fPeaks(jj) & ii > 1
                tracked = [tracked 0];
            else
                tracked = [tracked 1];
            end
        end
        
        disp([ 'fPeaks = ' num2str(fPeaks.')]);
        disp([ 'f0     = ' num2str(f0(map).')]);
        disp(['tracked = ' num2str(tracked)]);
        
        modeShapes = Hv(fLocs(map));
        modeShapes(~tracked) = Hv(peakPositions(~tracked));
        velocitiesMatrix(ii,:)= modeShapes;
end

%% store velocities in the grid points

velocityData = readmatrix('velocityPoints.csv');
velocityData(:,3) = [];
velocityData = velocityData;

velocityData = [velocityData velocitiesMatrix];


xyCoordVelSorted = sortrows(velocityData);

realVel = real(velocityData);
imagVel = [real(velocityData(:,1:2)), imag(velocityData(:,3:end))];
velocityData = sortrows(realVel) + 1i*(sortrows(imagVel));
velocityData(:,1:2) = real(velocityData(:,1:2));

velocityLabel ={'x' 'y'};
freqLabel = round(fPeaks);

for ii = 1: length(fPeaks)
    velocityLabel{ii+2} = ['f_{',num2str(ii),'} = ', num2str(freqLabel(ii)), ' Hz'];
end

velocityTable = array2table(velocityData, 'variableNames', velocityLabel);

cd(csvDataFolder)
writeMat2File(velocityData, velocityFilename, velocityLabel, length(velocityLabel), true);
cd(baseFolder)
%% Take a look on the velocities

velocityData = readmatrix(velocityFilename);
for ii = 1:nPeaks
    
    [X,Y,surfV] = getVelocityGroundtruth(velocityData(:,ii+2), velocityFilename, ii);
    title(['f_{', int2str(ii),'} = ', num2str(fPeaks(ii)), 'Hz'])
    

end

%     figure(ii+100)
%     plot3(velocityData(:,1), velocityData(:,2), abs(velocityData(:,ii+2)),'.')


%% Once you have the velocity, link the velocity 
  % on the measured grid of the violin - needed to interpolate on the
  % ESM resulting surface 

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
end 
cd(matDataFolder)
save('zDataVInMeasurements', 'zData');
cd(baseFolder)
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