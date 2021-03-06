function [HvSVD, singularVals] = SVD(frf, freq, M, nUsedSingVals, plotData)
% SVD - Singular Values Decomposition - reduces noise
% computes SVD of an FRF thourgh the hankel matrix, uses nUsedSingVals to
% reconstruct

%   INPUTS
%   frf          (1Darray)    = FRF on which perform SVD;
%   freq         (1Darray)    = frequency axis relative to the FRF;
%   M            (int)        = max number of singular values computed;
%   nUsedSingVals (double)    = number of singular values used to reconstruct 
%   plotData     (boolean     = true if you want to generate images;

%   OUTPUTS
%   HvSVD        (2Darray)    = H1 estimator after SVD;
%   singularVals (array)      = sigualar values;

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

singVals = singularVals(1:nUsedSingVals);
nRec = length(singVals);
Vp = V';

% d) Reconstruction of the cleaned Hankel matrix
H = U(:,1:nRec)*S(1:nRec,1:nRec)*Vp(1:nRec,:);

% e) Plot singular values and threshold
if plotData
    figure(17)
    stem(axis,singularVals);
    hold on
    stem(1:length(singVals), singVals)
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

%g) show results
if plotData
    figure(18)
    semilogy(freq, abs(frf(1:length(freq))));
    hold on
    semilogy(freq,abs(HvSVD(1:length(freq))));
    xlabel('f');
    ylabel('|H|');
    title(['Noisy and cleaned frfs ']);
    hold off
end
end

