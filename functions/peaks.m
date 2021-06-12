function [peaks, peaksLoc] = peaks(cleanMobilityMatrix, plotData)
%PEAKS find peaks of the FRF
%   Detailed explanation goes here

temp = abs(cleanMobilityMatrix);
temp = temp ./ mean(temp);
temp = sum(temp, 2);
[peaks, peaksLoc] = findpeaks(temp);

end

