% close all
% clear all
% clc
%% Preprocessing for NAH

% this script is used to compute the pressure data in the hologram points
% from measurement data

clear all;
close all;
addpath(genpath('Exp_Measurements')); % contains the hologram measurements
addpath(genpath('Data'))
addpath(genpath('functions'))

% store in a cell the set of measurements of the  8 mic array 
% +1 reference microphone over 8.
nMics = 8 + 1;
nMeasurements = 8;          % number of measurements, moving the array along the vertical axis
nTakes = 7;                 % number of takes for each measurement
mics_Data = cell(nMeasurements,nMics); % 8 vertical position and 9 channels

for ii = 1:nMeasurements % for each vertical position
    for jj = 1: nMics
        micFilename = ['y0',int2str(ii),'-00',int2str(jj),'.wav'];
        mics_Data{ii,jj} = audioread(micFilename);
    end
end

% store in a cell the set of measurements of the hammer and another
% reference microphone
% (force, acceleration and sound pressure)
impulseData = cell(nMeasurements,nTakes);
% cosa sono gli indici ? 
for ii = 1:nMeasurements
    for jj = 1:nTakes
        forceFilename = ['28_05_backplate_y', int2str(ii), '_', int2str(jj-1)];
        impulseData{ii,jj} = table2array(readtable([forceFilename,'.csv']));        
    end
end

% save the pressure and force cells into .mat files
save('mics_Data.mat', 'mics_Data', '-v7.3'); % must be saved with older version because the file is too big
save('impulseData.mat','impulseData');

%% import mat files
load('mics_Data.mat'); %  microphone
load('impulseData.mat'); % acquisition system

forces = cell(nMeasurements,1);

for ii = 1:length(impulseData(:,1))
    tempImpulses = [];
    for jj = 1:length(impulseData(1, :))
        tempImpulses = [tempImpulses; impulseData{ii,jj}(:,1).'];
    end
    forces{ii} = tempImpulses;
end

save('forces.mat','forces');

%% cut the pressure measurement signals e

Fs = 48000; % sampling frequency
duration = 2;
signalLength = duration*Fs; % signal length (the same as the force signal)
micSignals = cell(nMeasurements,nMics); % store the pressure signals
t = linspace(0, duration , signalLength); % time axis to plot 

for ii = 1:nMeasurements % for each vertical configuration    
    ReferenceMic = mics_Data{ii,nMics}; % take the reference microphone signal
    % find the highest peaks % che è min peak distance ? 
    [pks , loc] = findpeaks(ReferenceMic,'MinPeakDistance',2*Fs,'MinPeakHeight',0.1);  
    indexesCut = zeros(1,nTakes); % initiate the vector of indexes for cutting the signals
    % find the sample index relative to the peaks 
    for jj = 1:nTakes 
        indexesCut(jj) = findchangepts(ReferenceMic((loc(jj) - 0.01*Fs):(loc(jj) + 0.01*Fs)),'Statistic','std');
        indexesCut(jj) = indexesCut(jj) + loc(jj) - 0.01*Fs;
    end

    % init the temporary vector of cut signals
    cutSignalTemp = zeros(nTakes,signalLength); 
    for jj = 1:nMics % for each microphone signal
        signalToCut = mics_Data{ii,jj}; 
        for kk = 1:nTakes
            cutSignalTemp(kk,:) = signalToCut(indexesCut(kk):indexesCut(kk) + signalLength -1); % create 7 signal  
        end
        micSignals{ii,jj} = cutSignalTemp; % store them 
    end
end

%% apply the temporal exponential filter

alpha = 30; % exponential factor

micSignalsFiltered = cell(nMeasurements,nMics);
forceSignalsFiltered = cell(nMeasurements,1);
offsetFilter = 9600;

for ii = 1:nMeasurements % to the pressure signals
    for jj = 1:nMics
        filteredSignalMatrix = zeros(nTakes, signalLength);
        for kk = 1:nTakes           
             filteredSignalTemp = micSignals{ii,jj}(kk,:);
             expFilter = circshift(exp(-alpha*t),offsetFilter);
             expFilter(1:offsetFilter) = ones(1,offsetFilter);
             filteredSignalTemp = filteredSignalTemp.*expFilter;
             filteredSignalMatrix(kk,:) = filteredSignalTemp; 
        end
        micSignalsFiltered{ii,jj} = filteredSignalMatrix;
    end
end
figure(219)
semilogy(t, filteredSignalTemp/max(filteredSignalTemp), t ,expFilter);
xlim([0, 1]);
pause(0.01)

for ii = 1:nMeasurements % to the force signals
    filteredSignalTemp = forces{ii};
    filteredSignalTemp = filteredSignalTemp.*(exp(-alpha*t));
    forceSignalsFiltered{ii} = filteredSignalTemp;
end

% debug plot
figure(102)
subplot 211
a = randi([1 nMeasurements],1,1);
b = randi([1 nMics],1,1);
c = randi([1 nTakes],1,1);
plot(t, micSignalsFiltered{a,b}(c, :))
hold on
plot(t, exp(-alpha*t))
hold off
title(['y ',num2str(a),', mic: ', num2str(b),', meas: ', num2str(c)])

subplot 212
plot(t, forceSignalsFiltered{a}(c,:))
hold on
plot(t, exp(-alpha*t))
hold off
title(['force: y ',num2str(a),', meas: ', num2str(c)])

%% compute the H1 estimator
nMics_acq = 8;
nMeasurements = 8;
fCut = 2000;
H1 = cell(nMeasurements,1);
% apply a low pass filter to 3000 Hz
               

for ii = 1:nMeasurements
    H1Temp = [];
    for jj = 1:nMics_acq
        % pressureTemp = micSignalsFiltered{ii,jj};
        pressureTemp = micSignalsFiltered{ii,jj};
        forceTemp = forceSignalsFiltered{ii};
        
        [pxy, f] = cpsd(pressureTemp.', forceTemp.',[],[],signalLength, Fs);
        [pxx, f] = cpsd(forceTemp.',forceTemp.',[],[],signalLength, Fs);
        
        cutIdxs = find(f <fCut);
        f = f(cutIdxs); pxx = pxx(cutIdxs,:); pxy = pxy(cutIdxs,:);
        
        H1_est = sum(pxy,2)./sum(pxx,2);
%         figure(333)
%         semilogy(f, abs(H1_est));
%         xlim([0,500]);
%         ylim([0, 1e-1]);
        H1Temp = [H1Temp, H1_est];
    end
    
    H1{ii} = H1Temp( cutIdxs,:);
end

for ii = 1:nMeasurements % plot of H1 for each mic and vertical postion
    
    figure(ii)
    subplot 211
    semilogy(f, abs(H1{ii}), 'lineWidth', 1.2)
    xlim([0 1500])
    legend('mic1','mic2','mic3','mic4','mic5','mic6','mic7','mic8');
    title(['quota ', num2str(ii)]);
    
    colors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F','#000'};
   colororder(colors);
   
end

%% applying the SVD to H1 to reduce the noise

H1_cleaned = [];
singValsNum = 17;
usedSingVals = 1;
% to decide the threshold and n of singVals 
Hcheck = H1{2}(:,4);
[H1_est,singularVals] = SVD(Hcheck, f,singValsNum, usedSingVals, true);

for ii = 1:nMeasurements
    H1_temp = [];
    for jj = 1:nMics_acq
        [H1_est,singularVals] = SVD(H1{ii}(:,jj), f, singValsNum, usedSingVals, false);
        H1_temp = [H1_temp, H1_est];
    end
    H1_cleaned = [H1_cleaned, H1_temp];
end

for ii = 1:nMeasurements % plot the reduced-noise H1 estimator
    
    figure(ii)
    subplot(2,1,2)
    semilogy(f, abs(H1_cleaned(:,((ii-1)*nMics+1):ii*nMics)).', 'lineWidth', 1.2)
    xlim([0 1500])
    legend('mic1','mic2','mic3','mic4','mic5','mic6','mic7','mic8');
    title(['quota ', num2str(ii)]);
    
    colors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F','#000'};
    colororder(colors);
end


% save the H1 estimator (both with and without SVD)

% save('H1.mat','H1');
% 
% save('H1_cleaned.mat','H1_cleaned');

%% import data

% load('f.mat');
% load('H1.mat');
% load('H1_cleaned.mat');
% load('forces.mat');

%% Extract eigenfrequencies from pressure

% min peak prominence and min peak width B - high setting
highPeaksParams = [120,20];
% min peak prominence and min peak width A - low setting
lowPeaksParams = [15,3];
% fCut low for peaks to ignore
ignorePeaksLow = 100;
% fCut high for peaks to ignore
ignorePeaksHigh = 1410;
% frequency threshold for peak finders A - B configuration
fThreshold = 258;

[peakPositions, fPeaks] = peaks(H1_cleaned, f, fThreshold,...
    ignorePeaksLow, ignorePeaksHigh, highPeaksParams, lowPeaksParams);

%% PressureField
nPeaks = length(peakPositions);
% to adjust the parameters before
[Hv,f0, fLocs, csis, Q, pressureVals] = EMASimple(H1_cleaned(:,ii), f, 1e-4, 16, false);
pressureMatrix = zeros(nMeasurements*nMics, nPeaks);
for ii = 1:nMics*nMeasurements
        [Hv,f0, fLocs, csis, Q, ~] = EMASimple(H1_cleaned(:,ii), f, 1e-4, 16, false);
        
        map = [];
        tracked = [];
        for jj = 1:nPeaks
            diff = abs(f0 - fPeaks(jj));
            [minVal, minLoc] = min(diff);
            map = [map minLoc];
            
            if minVal >= 0.1*fPeaks(jj)
                tracked = [tracked 0];
            else
                tracked = [tracked 1];
            end
        end
        
        disp([ 'fPeaks = ' num2str(fPeaks.')]);
        disp([ 'f0     = ' num2str(f0(map).')]);
        disp(['tracked = ' num2str(tracked)]);
        
        pressureVals = Hv(fLocs(map));
        pressureVals(~tracked) = abs(Hv(peakPositions(~tracked)));
        pressureMatrix(ii,:)= pressureVals;
end

%% Pressure field 

highCut = 2000;
lowCut = 50;
cutIdxs = find(f <highCut & f>lowCut );
fAxis = f(cutIdxs);
yElements = 8;
xElements = 8;

figure(474)

for kk = 1:xElements
    H1cleanplot = H1_cleaned(cutIdxs,(kk-1)*nMics +1:kk*nMics);    
    semilogy(fAxis, abs(H1cleanplot).');
    title('H1 w/ SVD')
    hold on
end
xline(fPeaks)
hold off

% % take the velocity (H1 value) at that resonance
% pressure = zeros(nMics,nMeasurements, length(fPeaks));
% 
% 
% for ii = 1:xElements*yElements
%     for kk = 1:length(fPeaks)
%          pressure(ii,jj,kk) = H1_cleaned(peakPositions(kk),ii);
%     end
% end


xHologram = flip([176, 126, 74, 23, -27, -76 -126, -178]);
yHologram = linspace(0, 8*51.5, 8);
yHologram = yHologram - 211.3280;

[X,Y] = meshgrid(xHologram, yHologram);

pressureData = sortrows([X(:) Y(:)],2)
pressureData = [pressureData pressureMatrix];
% collect pressures in a matrix
for ii = 1:length(fPeaks)
    figure(55)
    surf(X, Y, reshape(abs(pressureData(:,2+ii)), nMics, nMeasurements).');
    pause(0.5);
end

realPress = real(pressureData);
imagPress = [real(pressureData(:,1:2)), imag(pressureData(:,3:end))];
pressureData = sortrows(realPress) + 1i*(sortrows(imagPress));
pressureData(:,1:2) = 0.001* real(pressureData(:,1:2));  

% [xSorted idxSorted] = sort(xyCoord(:,1));
% 
% xyCoord(:,2:end) = xyCoord(idxSorted,2:end);
% xyCoordSorted = flip(xyCoord);

pressureNames ={'x' 'y'};
freqLabel = round(fPeaks,1);

for ii = 1: length(fPeaks)
    pressureNames{ii+2} = ['f_{',num2str(ii),'} = ', num2str(freqLabel(ii))];
end

writeMat2File(pressureData, ['pressure_Data2.csv'], pressureNames , length(pressureNames), true);
