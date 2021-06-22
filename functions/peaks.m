function [peaksLoc, fpeakPositions] = peaks(matrix, f, fThreshold, ignorePeaksLow, ignorePeaksHigh,...
    highPeaksParams, lowPeaksParams)
%PEAKS find peaks of the fRF

%   INPUTS
%   matrix           (2Darray) = matrix of the H1 estimator;
%   f                (array)   = rfequency axis;
%   fThreshold       (double)  = frequency thershold for peak fining;
%   ignorePeaksLow   (double)  = parameter for ignore the low peaks for findpeaks function;
%   ignorePeaksHigh  (double)  = parameter for ignore the high peaks for findpeaks function;
%   highPeaksParams  (array)   = parameter for findpeaks function;
%   lowPeaksParams   (array)   = parameter for findpeaks function;

%   OUTPUTS
%   peaksLoc         (array)   = peaks location values;
%   fpeakPositions   (array)   = peaks location indices;

absMat = abs(matrix);
absMat = absMat ./ mean(absMat);
absMat = sum(absMat, 2).^2;

[~,low2HighIdx] = min(abs(f - (fThreshold*ones(size(f))) ));
[~, peaksLocHigh] = findpeaks(absMat(low2HighIdx+1:end), 'MinPeakProminence', highPeaksParams(1), 'minPeakWidth', highPeaksParams(2));
[~, peaksLocLow] = findpeaks(absMat(1:low2HighIdx), 'MinPeakProminence', lowPeaksParams(1), 'minPeakWidth', lowPeaksParams(2));

peaksLocHigh =low2HighIdx + peaksLocHigh;
peaksLoc = [ peaksLocLow; peaksLocHigh];
peaksLoc = sortrows(peaksLoc);

figure(901)
semilogy(f, absMat)
hold on
xline(f(peaksLoc))
hold off

fpeakPositions = f(peaksLoc);
idxPks = find(fpeakPositions < ignorePeaksLow |  fpeakPositions > ignorePeaksHigh);
if ~isempty(idxPks)
fpeakPositions(idxPks) = [];
peaksLoc(idxPks) = [];
end
end

