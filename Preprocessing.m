close all
clear all
clc
%% Preprocessing for NAH

% this script is used to compute the pressure data in the hologram points
% from measurement data

% import raw files of pressure measurements

addpath('measurements'); % contains the hologram measurements
addpath('measurements/forceMeasurements'); % contains the acquisition system measurements

% store in a cell the set of measurements of the microphone array and
% reference microphone.
arrayData = cell(8,9); % 8 vertical position and 9 channels

for i = 1:8
    for j = 1:9
        arrayData{i,j} = audioread(strcat('y0',int2str(i),'-00',int2str(j),'.wav'));
    end
end

% store in a cell the set of measurements of the hammer and another
% reference microphone
% (force, acceleration and sound pressure)
impulseData = cell(8,7);

for i = 1:8
    for j = 1:7
        forceData = importdata(strcat('28_05_backplate_y', int2str(i), '_', int2str(j-1), '.txt'));
        forceData = strrep(forceData, ',', '.');
        forceMes = zeros(length(forceData), 3);
        for k = 1:length(forceData)
            forceMes(k,:) = str2num(forceData{k});
        end
        impulseData{i,j} = forceMes;
    end
end


% save the pressure and force cells into .mat files
save('arrayData.mat', 'arrayData', '-v7.3'); % must be saved with older version because the file is too big
save('impulseData.mat','impulseData');

%% import mat files
load('arrayData.mat'); %  microphone
load('impulseData.mat'); % acquisition system

forces = cell(8,1);
for ii = 1:length(impulseData(:,1))
    tempImpulses = [];
    for jj = 1:length(impulseData(1, :))
        tempImpulses = [tempImpulses; impulseData{ii,jj}(:,1).'];
    end
    forces{ii} = tempImpulses;
end

save('forces.mat','forces');

%% cut the pressure measurement signals 

Fs = 48000; % sampling frequency
duration = 2;
signalLength = duration*Fs; % signal length (the same as the force signal)
micSignals = cell(8,9); % store the pressure signals

for aa = 1:8 % for each vertical configuration
    
    ReferenceMic = arrayData{aa,9}; % take the reference microphone signal
  
    [pks , loc] = findpeaks(ReferenceMic,'MinPeakDistance',2*Fs,'MinPeakHeight',0.1);  % find the highest peaks
    
    indexesCut = zeros(1,7); % initiate the vector of indexes for cutting the signals
    
    for j = 1:7 % find the initial starting samples of the impulse 
        indexesCut(j) = findchangepts(ReferenceMic((loc(j) - 0.01*Fs):(loc(j) + 0.01*Fs)),'Statistic','std');
        indexesCut(j) = indexesCut(j) + loc(j) - 0.01*Fs;
    end
    
    cutSignalTemp = zeros(7,signalLength); % init the temporary vector of cut signals
    
    for bb = 1:9 % for each microphone signal
        signalToCut = arrayData{aa,bb}; 
        for cc = 1:7
            cutSignalTemp(cc,:) = signalToCut(indexesCut(cc):indexesCut(cc) + signalLength -1); % create 7 signal
            % debug
            figure(101)
            t = linspace(0, signalLength/Fs , signalLength); % time axis to plot 
            subplot(2,4,cc)
            plot(t, cutSignalTemp(cc,:))         
        end
        micSignals{aa,bb} = cutSignalTemp; % store them 
    end
    
end

%% apply the temporal exponential filter

alpha = 10; % exponential factor

micSignalsFiltered = cell(8,9);
forceSignalsFiltered = cell(8,1);

for i = 1:8 % to the pressure signals
    for j = 1:9
        filteredSignalVector = zeros(7, signalLength);
        for k = 1:7           
             filteredSignalTemp = micSignals{i,j}(k,:);
             filteredSignalTemp = filteredSignalTemp.*(exp(-alpha*t));
             filteredSignalTemp =  lowpass(filteredSignalTemp,3000,Fs); % apply a low pass
                                                                        % filter to 3000 Hz
             filteredSignalVector(k,:) = filteredSignalTemp;
             
        end
        micSignalsFiltered{i,j} = filteredSignalVector;
    end
end

for ii = 1:8 % to the force signals
    filteredSignalTemp = forces{ii};
    filteredSignalTemp = filteredSignalTemp.*(exp(-alpha*t));
    forceSignalsFiltered{ii} = filteredSignalTemp;
end

% debug plot
figure(102)
subplot 211
a = randi([1 8],1,1);
b = randi([1 9],1,1);
c = randi([1 7],1,1);
plot(t, micSignalsFiltered{a,b}(c, :))
hold on
plot(t, exp(-alpha*t))
hold off
title(strcat('y ',num2str(a),', mic: ', num2str(b),', meas: ', num2str(c)))

subplot 212
plot(t, forceSignalsFiltered{a}(c,:))
hold on
plot(t, exp(-alpha*t))
hold off
title(strcat('force: y ',num2str(a),', meas: ', num2str(c)))

%% compute the H1 estimator

nMics = 8;
nHeights = 8;

H1 = cell(nHeights,1);

for ii = 1:nHeights
    H1Temp = [];
    for jj = 1:nMics
        pressureTemp = micSignalsFiltered{ii,jj};
        forceTemp = forceSignalsFiltered{ii};
        
        [pxy, f] = cpsd(pressureTemp.', forceTemp.',[],[],signalLength, Fs);
        [pxx, f] = cpsd(forceTemp.',forceTemp.',[],[],signalLength, Fs);
        
        cutIdxs = find(f <2000);
        f = f(cutIdxs); pxx = pxx(cutIdxs,:); pxy = pxy(cutIdxs,:);
        
        H1_est = sum(pxy,2)./sum(pxx,2);
        H1Temp = [H1Temp, H1_est];
    end
    
    H1{ii} = H1Temp( cutIdxs,:);
end

for ii = 1:nHeights % plot of H1 for each mic and vertical postion
    
    figure(ii)
    subplot(2,1,1)
    semilogy(f, abs(H1{ii}), 'lineWidth', 1.2)
    xlim([0 1500])
    legend('mic1','mic2','mic3','mic4','mic5','mic6','mic7','mic8');
    title(['quota ', num2str(ii)]);
    
    colors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F','#000'};
   colororder(colors);
end

%% applying the SVD to H1 to reduce the noise

H1_cleaned = cell(8,1);
M = 5;
thresholdPerc = 5;

for ii = 1:nHeights
    H1_temp = [];
    for jj = 1:nMics
        [H1_est,threshold,singularVals] = SVD(H1{ii}(:,jj), f, M, thresholdPerc, false);
        H1_temp = [H1_temp, H1_est];
    end
    H1_cleaned{ii} = H1_temp;
end

for ii = 1:8 % plot the reduced-noise H1 estimator
    
    figure(ii)
    subplot(2,1,2)
    semilogy(f, abs(H1_cleaned{ii}), 'lineWidth', 1.2)
    xlim([0 1500])
    legend('mic1','mic2','mic3','mic4','mic5','mic6','mic7','mic8');
    title(['quota ', num2str(ii)]);
    
    colors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F','#000'};
    colororder(colors);
end

% take a look at the SVD in a signal
Hcheck = H1{2}(:,4);
[H1_est,threshold,singularVals] = SVD(Hcheck, f, M, thresholdPerc, true);

% save the H1 estimator (both with and without SVD)

save('H1.mat','H1');

save('H1_cleaned.mat','H1_cleaned');

%% import data

load('f.mat');
load('H1.mat');
load('H1_cleaned.mat');
load('forces.mat');

figure(474)

for k = 1:8
    
    H1plot = H1{k};
    H1cleanplot = H1_cleaned{k};
    
    subplot 211
    semilogy(f, abs(H1plot));
    title('H1 w/o SVD')
    
    subplot 212
    semilogy(f, abs(H1cleanplot));
    title('H1 w/ SVD')
    
    hold on
end

hold off
%% Pressure field 

addpath('functions')
Fs = 48000;
cutIdxs = find(f <2000 & f>5 );
f = f(cutIdxs);

yElemends = 8;
xElements = 8;

% take one single resonance frequency anound the interval 
resInterval = find(f>360&f<390); % indx

% take the velocity (H1 value) at that resonance
pressure = zeros(8,8);

checkEstimator = H1_cleaned{4};
checkEstimator = checkEstimator(:,4);
[pks, locs] = findpeaks(abs(checkEstimator(resInterval)), f(resInterval),'MinPeakHeight',1e-3);

resFreq = locs;
resFreq = 377;
resIdx = find(f == resFreq);

for jj = 1:yElemends
    for ii = 1:xElements
        pressure(ii,jj) = H1_cleaned{jj}(resIdx,ii);
    end
end

xHologram = [176, 126, 74, 23, -27, -76 -126, -178];
yHologram = linspace(0, 8*51.5, 8);
yHologram = yHologram -   211.3280;

[X,Y] = meshgrid(xHologram, yHologram);
figure(908)
surf(X,Y,abs(pressure))
xlabel('x'); ylabel('y')

