function [violinInfos, velocityFields, hologramInfos, pressureFields] = importData(velocityFileName, pressureFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION READS THE CSV OF PRESSURE AND NORMAL  %%%%%%%%%%%%
%%% VELOCITY, EXTRACTS THEIR DATA AND STORES THEM       %%%%%%%%%%%%
%%% IN CELLS                                            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reading csv
csvVel = readtable(velocityFileName);
csvPress = readtable(pressureFileName);
violinMesh = table2array(csvVel(:,1:3));
hologramMesh = table2array(csvPress(:,1:3));

%% Reconconstruct the plate geometry
gridX = length(unique(violinMesh(:,1)));
gridY = length(unique(violinMesh(:, 2)));

% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> 
% ex. try  reshape(1:10,[5,2])  what we actuallly want is  reshape(1:10,[2,5])'

X =  reshape(violinMesh(:,1), [gridY, gridX]).'; 
Y =  reshape(violinMesh(:,2), [gridY, gridX]).'; 
Z =  reshape(violinMesh(:,3), [gridY, gridX]).';

figure(1)
surf(X,Y,Z);
title('Violin surface');
zlim([0,100]);
violinInfos = {X,Y,Z}; 

%% Velocity Fields
numFreqBins = length(table2array(csvVel(1,:)))-3;
velsMatrix = table2array( csvVel(:,4:numFreqBins+3) );

velocityFields = cell(numFreqBins,1);
for ii = 1:numFreqBins
    velocityFields{ii} = reshape(velsMatrix(:,ii), [gridY, gridX]).';   
end

figure(2)
surf(X,Y,abs(velocityFields{1}));
title('Velocity field surface - f1');

%% Reconstruct the hologram geometry

gridX = length(unique(hologramMesh(:,1)));
gridY = length(unique(hologramMesh(:,2)));

% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> 
% ex. try  reshape(1:10,[5,2])  what we actuallly want is  reshape(1:10,[2,5])'

X =  reshape(hologramMesh(:,1), [gridY, gridX]).'; 
Y =  reshape(hologramMesh(:,2), [gridY, gridX]).'; 
Z =  reshape(hologramMesh(:,3), [gridY, gridX]).';

figure(3)
surf(X,Y,Z);
hologramInfos = {X,Y,Z}; 
title('Hologram surface');

%% Pressure Fields
numFreqBins = length(table2array(csvPress(1,:)))-3;
pressMatrix = table2array( csvPress(:,4:numFreqBins+3) );

pressureFields = cell(numFreqBins,1);
for ii = 1:numFreqBins
    pressureFields{ii} = reshape(pressMatrix(:,ii), [gridY, gridX]).';   
end

figure(4) 
surf(X,Y,abs(pressureFields{1}));
title('Pressure field surface - f1');
end

