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

estimatorMatrix = zeros(3989,18); % store the H1 estimator
cleanEstimatorMatrix = zeros(3989,18); % store the H1 estimator denoised
frfMatrix = zeros(3989,18); % store the standard FRF (fft(Y)/fft(X))

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
                   0 8];

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
    
    frfMatrix(:,jj) = mean(FRFTemp(cutIdxs,:),2);
    
end

figure(808)
subplot 211
semilogy(f, (abs(estimatorMatrix)))
title('H1 without SVD')
subplot 212
semilogy(f, (abs(frfMatrix)))
title('Y/X')

%% apply SVD to H1 estimator

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

%% find resonance frequencies
% impossible with our data :(
% peakMat = zeros(numberPoints, 28);
% 
% for k = 1:numberPoints
%     checkEstimator = cleanEstimatorMatrix(:, k);
%     checkEstimator = lowpass(checkEstimator, 1000 , 2*length(checkEstimator));  
% %     [Hv,f0, fLocs, csis, Q] = EMASimple(checkEstimator, f, 0.0001, 10,  true);
%     [pks, locs] = findpeaks(abs(checkEstimator), f,'MinPeakProminence',2e-4);
%     peakMat(k, 1:length(locs)) = locs';
% end

%% velocity field

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
msrPoints = [-2  0; -2 -2; -2  4; -2  2;
             -2 -4; -2 -6;  0  0;  2  0;
              0  2;  2  2;  0  4;  2  4;
              0 -2;  2 -2;  0 -2;  2 -4;  
              0 -6;  2  6];
          
nMsrPoints = length(msrPoints(:,1));

centerData = table2array(geomData(1:6,1:3));
XX = table2array(geomData(11:16,7:9));
YY = table2array(geomData(11:16,11:13));

xData = table2array(geomData(1:6,7:9)).*sign(XX);
yData = table2array(geomData(1:3,11:16)).'.*sign(YY);

zData = zeros(size(xData));

for ii = 1: nMsrPoints
    [xx, yy] = find(XX == msrPoints(ii,1) & YY == msrPoints(ii,2));
    zData(xx,yy) = velocities(ii);
end  

x = reshape(xData, nMsrPoints,1);
y = reshape(yData, nMsrPoints,1);
z = reshape(zData, nMsrPoints, 1);

figure()
subplot 121
plot3(x, y, z, '.');
subplot 122
surf(xData, yData, zData);

% write approach
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


  