close all
clear all
clc
%% Preprocessing for NAH
%{
% import raw files

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

% debug
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
%}
%% measured FRF of the acceleration

Fs = 48000;
duration = 2;
signalLength = duration*Fs;

numberAcquisitions = 6;

% SVD Parameters
M = 10;
thresholdPerc = 30;

% EstimatorCell = cell(6,3);
estimatorMatrix = zeros(3989,18);
cleanEstimatorMatrix = zeros(3989,18);
frfMatrix = zeros(3989,18);

listPath = {'-2_0_-2_-2', '-2_4_-2_2', '-2_-4_-2_-6', '0_0_2_0', '0_2_2_2', '0_4_2_4', '0_-2_2_-2', '0_-4_2_-4', '0_-6_2_-6'};
firstIndexes = [1,3,5,7,9,11,13,15,17];
secondIndexes = [2,4,6,8,10,12,14,16,18];

forceTemp = zeros(signalLength, numberAcquisitions);
firstTemp = zeros(signalLength, numberAcquisitions);
secondTemp = zeros(signalLength, numberAcquisitions);
firstFRFTemp = zeros(signalLength, numberAcquisitions);
secondFRFTemp = zeros(signalLength, numberAcquisitions);

for jj = 1:9
    currentPath = listPath{jj};
    addpath(strcat('accelerations/', currentPath))
    
    for ii = 1:numberAcquisitions
        
        rawTemp = importdata(strcat( currentPath, '_frf_acc', int2str(ii-1) , '.txt'));
        rawTemp = strrep(rawTemp, ',', '.');
        rawMes = zeros(length(rawTemp), 3);
        
        for k = 1:length(rawTemp)
            rawMes(k,:) = str2num(rawTemp{k});
        end
        
        forceTemp(:,ii) = rawMes(:,1);
        firstTemp(:,ii) = rawMes(:,2);
        secondTemp(:,ii) = rawMes(:,3);
        
        firstFRFTemp(:,ii) = fft(rawMes(:,2))./ftt(rawMes(:,1));
        secondFRFTemp(:,ii) = fft(rawMes(:,3))./ftt(rawMes(:,1));
        %%%%%% GUARDARE DOCUMEMTAZIONE DI MATçLAB PER LA FFT %%%%%%%
        
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
    
    
    
    [H1_clean,threshold,singularVals] = SVD(H1.', f, M, thresholdPerc, false);
    [H2_clean,threshold,singularVals] = SVD(H2.', f, M, thresholdPerc, false);
    
    % EstimatorCell{3,2} = H1_clean;
    % EstimatorCell{3,3} = H2_clean;
    
    estimatorMatrix(:,firstIndexes(jj)) = H1;
    estimatorMatrix(:,secondIndexes(jj)) = H2;
    
    cleanEstimatorMatrix(:,firstIndexes(jj)) = H1_clean;
    cleanEstimatorMatrix(:,firstIndexes(jj)) = H2_clean;
    
    
end

figure(808)
subplot 211
plot(f, 10*log10(abs(estimatorMatrix)))
title('H1 without SVD')
subplot 212
plot(f, 10*log10(abs(cleanEstimatorMatrix)))
title('H1 with SVD')

%% Pressure amplitude computation

% import frequency bins,  H1 estimator and forces
load('f.mat');
load('H1.mat');
load('H1_cleaned.mat');
load('forces.mat');

%%
% TEMPORARY CODE BELOW (WAITING FOR THE FRF)
Hcheck = H1{2}(:,4);
Fs = 48000;
thresholdPerc = 5;
addpath('functions')

[fAmps, fLocs] = findpeaks(abs(Hcheck),'minPeakProminence', 1e-2, 'MinPeakWidth',3);
f0 = f(fLocs);

 figure()
 plot(f,20*log10(abs(Hcheck)));
hold on
%  for ii = 1:length(fLocs)
%      stem(f0, abs(Hcheck(fLocs)));
%  end
 
%{
 fUse = fLocs(3);
 pressureMatrix = zeros(8,8);
 
 for ii = 1:8
    
    forceSig = mean(forces{ii});
    F = fft(forceSig);
    L = length(F);
    P2 = abs(F/L);
    P1 = P2(:,1:L/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    fAxis = Fs*(0:(L/2))/L;
    [F_clean,threshold,singularVals] = SVD(P1(1:3000), f(1:3000), 100, thresholdPerc, false);
    pressureMatrix(ii,:) = H1{ii}(fUse,:);
    
 end
 

%}
 %% velocity trial


addpath('accelerations/0_0_2_0')
velocityTrial = importdata('0_0_2_0_frf_acc0.txt');
velocityTrial = strrep(velocityTrial, ',', '.');
velocityMes = zeros(length(velocityTrial), 3);
for k = 1:length(velocityTrial)
    velocityMes(k,:) = str2num(velocityTrial{k});
end
accelerationZero = velocityMes(:,2);
hammerZero = velocityMes(:,1);
Fs = 48000;
w = 1:48000;
w = 2*pi*w;
w = w';
V = abs(1./(1i.*w).*fft(accelerationZero, Fs));
F = abs(fft(hammerZero, Fs));
FRF = 10*log10(V./F);
plot(w./(2*pi), FRF)
xlim([0 2000])
hold off