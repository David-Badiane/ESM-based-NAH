function [peaksLoc, fpeakPositions] = peaks(matrix, f, fThreshold, ignorePeaksLow, ignorePeaksHigh,...
    highPeaksParams, lowPeaksParams)
%PEAKS find peaks of the FRF
%   Detailed explanation goes here

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

