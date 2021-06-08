% Get velocity signals
% READ ACCELERATION 
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = ;             % Length of signal
t = (0:L-1)*T; 
f = Fs*(0:(L/2))/L;
w = 2*pi*f;

velocity = 1i*w*fft(acceleration);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

M = 10;
thresholdPerc = 30;
sampleName
counter 
[vSVD,threshold,singularVals] = SVD(velocity, f, M, 30, smpleName, jj);


[Hv,f0, fLocs, csis, Q] = EMASimple(vSVD, f,minPeakVal, minPeakWidth)
% Get pressure signals

fv = stlread('NAH_ESM_mesh_fin.stl');
points = sortrows(fv.vertices);
[outMatrix] = downsampling_regular(points, 128, 128)
[pts] = downsampling_regular(outMatrix, 128, 128)


