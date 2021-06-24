function [peaksLoc, fpeakPositions] = peaks(matrix, f, fThreshold, ignorePeaksLow, ignorePeaksHigh,...
    highPeaksParams, lowPeaksParams)
%PEAKS finds peaks of the FRF by analysing their cumulative sum
%   INPUTS
%   matrix           (2Darray) = matrix of the H1 estimator;
%   f                (1Darray)   = frequency axis;
%   fThreshold       (double)  = frequency threshold for peak fining (we use two find peaks);
%   ignorePeaksLow   (double)  = low threshold - ignore peaks at lower frequency than it;
%   ignorePeaksHigh  (double)  = high threshold - ignore peaks at higher frequency than it;
%   highPeaksParams  (1Darray) = [2x1] = minPeakProminence, minPeak width for highFreq findpeaks;
%   lowPeaksParams   (1Darray) = [2x1] = minPeakProminence, minPeak width for lowFreq findpeaks;
%   OUTPUTS
%   peaksLoc         (array)   = peaks location values;
%   fpeakPositions   (array)   = peaks location indices;

absMat = abs(matrix);
absMat = absMat ./ mean(absMat);
absMat = sum(absMat, 2).^2;

[~,low2HighIdx] = min(abs(f - (fThreshold*ones(size(f))) ));
[~, peaksLocHigh] = findpeaks(absMat(low2HighIdx+1:end), 'MinPeakProminence', highPeaksParams(1), 'minPeakWidth', highPeaksParams(2));
[~, peaksLocLow]  = findpeaks(absMat(1:low2HighIdx), 'MinPeakProminence', lowPeaksParams(1), 'minPeakWidth', lowPeaksParams(2));

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

