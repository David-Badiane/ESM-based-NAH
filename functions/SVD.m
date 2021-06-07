function [HvSVD, threshold,singularVals] = SVD(frf, freq, M, thresholdPerc, plotData)
% SVD Summary of this function goes here
%
% frf  (array)  = signal on which perform SVD,
% freq (array)  = frequency axis relative to the signal
% M    (int)    = max number of singular values computed
% thresholdPerc (double) = set threshold for useful singular values in percentage of
%                          the maximum value of the singularValues 
%                          (ex. thresholdPerc = 30 --> 30%)
% samplesNames (cell array) = cell array containing the names of the samples
% counter (int)             = value indicating the sample number

% a) Generate Hankel matrix and preallocation of variables
H = hankel(frf);
N = length(freq); % total number of bins
L = N - M + 1;   % nrows
H = H(1:L,1:M);  % the piece of the Hankel matrix we want
index = min([M,L]);
axis = 1:index;

% b) Singular values decomposition
[U,S,V] = svd(H);
singularVals = diag(S);

% c) Set threshold in percentage of maximum value between all singular values
maxSingVal = max(singularVals);
threshold = thresholdPerc/100*maxSingVal;
th = ones(size(axis))*threshold;
singVals = singularVals(singularVals>= threshold);
nRec = length(singVals);
Vp = V';

% d) Reconstruction of the cleaned Hankel matrix
H = U(:,1:nRec)*S(1:nRec,1:nRec)*Vp(1:nRec,:);

% e) Plot singular values and threshold
if plotData
    figure()
    stem(axis,singularVals);
    hold on;
    plot(axis,th);
    xlabel('n° singular value');
    ylabel('singular values');
    title(['SVD']);
end

% f) Retrieve the spectrum
HvSVD = zeros(size(frf));
for ii = 1:N
    k = min([M,ii]);
    l = max([1,ii-L+1]);
    if (k>l)
        sumOver = l:k;
    else
        sumOver = k:l;
    end
    sum = 0;
    for jj = sumOver
        sum = sum + H(ii-jj+1,jj);
    end
    HvSVD(ii) = 1/(k-l+1)*sum; 
end

if plotData
    figure()
    semilogy(freq, abs(frf(1:length(freq))));
    hold on
    semilogy(freq,abs(HvSVD(1:length(freq))));
    xlabel('f');
    ylabel('|H|');
    title(['Noisy and cleaned frfs ']);
end
end

