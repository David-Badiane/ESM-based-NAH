close all
clear all
clc
%% Ground truth computation

addpath('functions')
% compute FRF of the velocity


baseFolder = pwd;
accFolder = [baseFolder, '\acq_11_06'];

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

% temportary file to compute the frequency responces
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

save('velocityH1.mat','estimatorMatrix');
save('velocityH1cleaned.mat','cleanEstimatorMatrix');
save('velocityMobility.mat','MobilityMatrix');
save('velocityMobilitycleaned.mat','cleanMobilityMatrix');

figure(890)
semilogy(f, (abs(cleanMobilityMatrix)))
title('Mobility with SVD')

%% find resonance frequencies

referenceFreqValues = [283, 375, ... % approximate position of peaks in Hz
                       498,682, 729, 839, ...
                       983, 1197, 1374, 1630 ];
referenceFreqValues = referenceFreqValues -4;
                  
intervalResonance = 10; % Hz to look around the reference to find resonances

% peakPositions = zeros(numberPoints,length(referenceFreqValues)); % store the resonance frequencies
peakValues = zeros(numberPoints,length(referenceFreqValues));
% peakPositions = zeros(100);
% peakValues = zeros(100);

% frequency resolution 
fResolution = ((f(end) - f(1))/(length(f)-1));

% EMA simple parameters
MinPeakProminence = 2e-4;
MinPeakWidth = 20;

for k = 1:numberPoints
      checkMobility = cleanMobilityMatrix(:, k);
    
%         minSample = (referenceFreqValues(j)-intervalResonance)/fResolution;
%         maxSample = (referenceFreqValues(j)+intervalResonance)/fResolution;
%         [Hv,f0, fLocs] = EMASimple(checkMobility, f, MinPeakProminence, MinPeakWidth,  true);   
%         peakPositions = [peakPositions; fLocs];
%         peakValues = [peakValues; Hv(fLocs)];

%         [Hv,f0, fLocs] = EMASimple(checkMobility, f, MinPeakProminence, MinPeakWidth,  false); 
%         peakPositions(k,1:length(fLocs)) = fLocs.';
%         peakValues(k, 1:length(fLocs))= checkMobility(fLocs).';
           
          for jj = 1:length(referenceFreqValues)
              peakValues(k,jj) = checkMobility(round(referenceFreqValues(jj)/fResolution));
          end

end


figure(899)
semilogy(f, (abs(cleanMobilityMatrix)))
hold on
for hh = 1:numberPoints
%     iddx = find(peakPositions(hh,:) ~= 0);
    
    stem(f(round(referenceFreqValues/fResolution)), abs(peakValues(hh,:)))
end
hold off
title('Mobility with SVD')

%% store velocities in the grid points


%% THIS IS A BACKUP CODE for THE VIOLIN GRID

% take one single resonance frequency anound the interval 
resInterval = find(f>360&f<390); % indx

% take the velocity (H1 value) at that resonance
velocities = zeros(1,18);
peaks = zeros(1,18);

for ii = 1:length(velocities)
    checkEstimator = cleanEstimatorMatrix(:,ii);
    [pks, locs] = findpeaks(abs(checkEstimator(resInterval)), f(resInterval),'MinPeakProminence',2e-4);
    velocities(ii) = pks;
    peaks(ii) = locs;
    %debug
    figure()
    semilogy(f, abs(checkEstimator))
    hold on
    xline(locs)
    hold off
end


%% export the velocity grid

% auto approach 
geomData = readtable('geometryVelocity.xlsx');

centerData = table2array(geomData(1:5,1:3));
XX = table2array(geomData(12:21,7:13));
YY = table2array(geomData(12:21,16:22));

xData = table2array(geomData(1:10,7:13)).*sign(XX);
yData = table2array(geomData(1:10,16:22)).*sign(YY);


zData = zeros(size(xData));
orderedPoints = zeros(numberPoints, 3);

for ii = 1: numberPoints
    [xx, yy] = find(XX == measurementPts(ii,1) & YY == measurementPts(ii,2));
    orderedPoints(ii,1) = xData(xx,yy);
    orderedPoints(ii,2) = yData(xx,yy);
    disp(orderedPoints(ii,:));
    disp([xx,yy]);
%     orderedPoints(ii,3) = velocities(ii); 
%     zData(xx,yy) = velocities(ii);
end  

writeMat2File(orderedPoints, 'velocityPoints.csv', {'x' 'y' 'z'}, 3, true )
% x = reshape(xData, numberPoints,1);
% y = reshape(yData, numberPoints,1);
% z = reshape(zData, numberPoints, 1);
% 
% figure()
% subplot 121
% plot3(x, y, z, '.');
% subplot 122
% surf(xData, yData, zData);

%% write approach
xCoord = 10*[-5, -4.7, -4, -4.1, -4.8, -5.1, 0, 0, 0, 0, 0, 0, 4.9 ,5, 4.2, 4, 4.8, 5];
xCoord = xCoord';

yCoord = 10*[8, 4.6, 0, -4.1, -8, -12.1, 8, 5, 0, -4, -7.9, -11.9, 7.9, 4.9, 0 , -4.1, -8.1, -12.1];
yCoord = yCoord';

velocities = [velocities(3), velocities(4), velocities(1),  velocities(2),  velocities(5), ...
       velocities(6) ,   velocities(11),   velocities(9),   velocities(7),   velocities(13), ...
       velocities(15),   velocities(17),    velocities(12),   velocities(10),   velocities(8), ...
       velocities(14),   velocities(16),   velocities(18)]';
 
   
figure()
plot3(xCoord, yCoord, velocities, '.');

XX = reshape(xCoord, [3,6])'
YY = reshape(yCoord, [3,6])'
ZZ = reshape(velocities, [3,6])'
figure(999)
surf(XX,YY,ZZ)


  