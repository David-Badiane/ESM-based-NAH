function [peaksLoc] = peaks(matrix, f)
%PEAKS find peaks of the FRF
%   Detailed explanation goes here

thresh = 200*2;

temp = abs(matrix);
temp = temp ./ mean(temp);
temp = sum(temp.^2, 2);

figure(901)
semilogy(f, temp)

[~, peaksLocHigh] = findpeaks(temp(thresh:end, 1), 'MinPeakProminence', 10);
[~, peaksLocLow] = findpeaks(temp(1:thresh, 1), 'MinPeakProminence', 1);

peaksLocHigh = peaksLocHigh+thresh;

peaksLoc = cat(1, peaksLocLow, peaksLocHigh);
peaksLoc = sort(peaksLoc, 1);

end

