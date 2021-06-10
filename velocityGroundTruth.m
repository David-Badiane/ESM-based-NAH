close all
clear all
clc
%% Ground truth computation

% compute FRF of the velocity

Fs = 48000; % sampling frequency
duration = 2; % [s]
signalLength = duration*Fs;
numberPoints = 18; % number of measurement points in the violin plate
numberAcquisitions = 6; % number of measurements for each point

estimatorMatrix = zeros(3989,18); % store the H1 estimator
cleanEstimatorMatrix = zeros(3989,18); % store the H1 estimator denoised
frfMatrix = zeros(3989,18); % store the standard FRF (fft(Y)/fft(X))

% import raw data 

% list of file names that contains acquisitions
listPath = {'-2_0_-2_-2', '-2_4_-2_2', '-2_-4_-2_-6', '0_0_2_0', '0_2_2_2', '0_4_2_4', '0_-2_2_-2', '0_-4_2_-4', '0_-6_2_-6'};
firstIndexes = [1,3,5,7,9,11,13,15,17];
secondIndexes = [2,4,6,8,10,12,14,16,18];

forceTemp = zeros(signalLength, numberAcquisitions);
firstTemp = zeros(signalLength, numberAcquisitions);
secondTemp = zeros(signalLength, numberAcquisitions);
firstFRFTemp = zeros(signalLength/2, numberAcquisitions);
secondFRFTemp = zeros(signalLength/2, numberAcquisitions);


for jj = 1:numberPoints/2
    currentPath = listPath{jj};
    addpath(strcat('accelerations/', currentPath))
    
    for ii = 1:numberAcquisitions
        
        rawTemp = importdata(strcat( currentPath, '_frf_acc', int2str(ii-1) , '.txt'));
        rawTemp = strrep(rawTemp, ',', '.');
        rawMes = zeros(length(rawTemp), 3);
        
        for k = 1:length(rawTemp) % convert txt file into 
            rawMes(k,:) = str2num(rawTemp{k});
        end
        
        forceTemp(:,ii) = rawMes(:,1);
        firstTemp(:,ii) = rawMes(:,2);
        secondTemp(:,ii) = rawMes(:,3);
        
        % compute the FFT for each acquisition
            Y1 = fft(firstTemp(:,ii));
            Y2 = fft(secondTemp(:,ii));
            X = fft(forceTemp(:,ii));
            
            Y1 = abs(Y1/signalLength); Y1_cut = Y1(1:signalLength/2); Y1_cut(2:end-1) = 2*Y1_cut(2:end-1);
            Y2 = abs(Y2/signalLength); Y2_cut = Y2(1:signalLength/2); Y2_cut(2:end-1) = 2*Y2_cut(2:end-1);
            X = abs(X/signalLength); X_cut = X(1:signalLength/2); X_cut(2:end-1) = 2*X_cut(2:end-1);
            
            freq = Fs*(0:(signalLength/2)-1)/signalLength;          
            freq(1) = 1;
            freq = freq';
            firstFRFTemp(:,ii) =  (1./(1i*2*pi*freq)).*(Y1_cut./X_cut);
            secondFRFTemp(:,ii) =  (1./(1i*2*pi*freq)).*(Y2_cut./X_cut);
        
    end
    
    [pxy1, f] = cpsd(firstTemp, forceTemp,[],[],signalLength, Fs);
    [pxy2, f] = cpsd(secondTemp, forceTemp,[],[],signalLength, Fs);
    [pxx, f] = cpsd(forceTemp,forceTemp,[],[],signalLength, Fs);
    
        
    cutIdxs = find(f <2000 & f>5 );
    f = f(cutIdxs); pxx = pxx(cutIdxs,:);pxy1 = pxy1(cutIdxs,:); pxy2 = pxy2(cutIdxs,:);
    
    pxy1 = 1./(1i*2*pi*f).*pxy1;
    pxy2 = 1./(1i*2*pi*f).*pxy2;
    
    H1 = sum(pxy1,2)./sum(pxx,2);
    H2 = sum(pxy2,2)./sum(pxx,2);
    
    estimatorMatrix(:,firstIndexes(jj)) = H1;
    estimatorMatrix(:,secondIndexes(jj)) = H2;
    
    frfMatrix(:,firstIndexes(jj)) = mean(firstFRFTemp(cutIdxs,:),2);
    frfMatrix(:,secondIndexes(jj)) = mean(secondFRFTemp(cutIdxs,:),2);
    
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

for ii = 1:length(velocities)
    checkEstimator = cleanEstimatorMatrix(:,ii);
    [pks, locs] = findpeaks(abs(checkEstimator(resInterval)), f(resInterval),'MinPeakProminence',2e-4);
    velocities(ii) = pks;
    %debug
%     figure()
%     semilogy(f, abs(checkEstimator))
%     hold on
%     xline(locs)
%     hold off
end

