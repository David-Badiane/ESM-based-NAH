function [peaksLoc] = peaks(cleanMobilityMatrix)
%PEAKS find peaks of the FRF
%   Detailed explanation goes here

temp = abs(cleanMobilityMatrix);
temp = temp ./ mean(temp);
temp = sum(temp, 2);
[~, peaksLoc] = findpeaks(temp);

end

