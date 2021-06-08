function createHologram(pressureMatrix, yCoordinate)
%% this function creates the hologram of pressure points and store it in a .csv

% hologram geometry. measures in mm
 xHologram = [-178, 128, -76, -27, 23, 74, 126, 176]; % from measurements with error of +-2 mm CHECK THIS??
 yHologram = linspace(0, 7*51.5, 8); % from measurements with error of +-2 mm
 yHologram = yHologram - yCoordinate; % shift the zero to the center of the violin. The distance is 360.5 vs. the measured one in 362
 zHologram = 20*ones(8,1);
 
 numberFrequencies = length(pressureMatrix);
 % use writeMat2File !!

end

