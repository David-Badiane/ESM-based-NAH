function createHologram(pressureMatrix, yCoordinate)
%% this function creates the hologram of pressure points and store it in a .csv

% hologram geometry. measures in mm
 xHologram = linspace(-178, 176, 8); % from measurements with error of +-2 mm
 centerCoordinateY = 133; % measured distance of the center to the bottom of the plate
 yHologram = linspace(0, 7*51.5, 8); % from measurements with error of +-2 mm
 yHologram = yHologram - centerCoordinateY; % shift the zero to the center of the violin
 zHologram = 20*ones(8,1);

end

