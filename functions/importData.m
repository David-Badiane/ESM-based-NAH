function [violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = importData(velocityFileName, pressureFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION READS THE CSV OF PRESSURE AND NORMAL VELOCITY,             %%%%%%%%%%%%
%%% EXTRACTS THEIR DATA AND STORES THEM IN CELL ARRAYS                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                   %%%%%%%%%%%%
%%% velocityFileName = filename of the velocity file (string)                %%%%%%%%%%%%
%%% pressureFileName = filename of the pressure file (string)                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                  %%%%%%%%%%%%
%%% violinInfos = violin {X mesh, Y mesh, Z mesh, unrolled mesh}(cell array) %%%%%%%%%%%%
%%% veleocityFields = velocity fields matrices (cell array)                  %%%%%%%%%%%%
%%% hologramInfos = hologram meshes, same structure of violinMesh(cell array)%%%%%%%%%%%%
%%% pressureFields = pressure fields matrices (cell array)                   %%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reading csv
csvVel = readtable(velocityFileName);
csvPress = readtable(pressureFileName);
violinMesh = table2array(csvVel(:,1:3));
hologramMesh = table2array(csvPress(:,1:3));

%% Import eigenfrequencies 

    eigenFreqz = char(csvPress.Properties.VariableNames(4:end));
    eigenFreqzdouble = zeros(size(eigenFreqz(:,1)));
    eigenFreqz = eigenFreqz(:,3:end);
    indexes = find(eigenFreqz == '_');
    eigenFreqz(indexes) = '.';
    
    for ii = 1:length(eigenFreqz(:,1)) 
        eigenFreqzdouble(ii, :) = str2double(eigenFreqz(ii, :));
    end
    
    eigenFreqz = eigenFreqzdouble;


%% Reconconstruct the plate geometry

gridX = length(unique(violinMesh(:, 1)));
gridY = length(unique(violinMesh(:, 2)));

% clean the violin border to accident zeros
violinMeshZ = violinMesh(:, 3);
violinMeshZ( violinMeshZ <= 0) = NaN;
violinMesh(:,3) = violinMeshZ;

% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> 
% ex. try  reshape(1:10,[5,2])  what we actuallly want is  reshape(1:10,[2,5])'

X =  reshape(violinMesh(:,1), [gridY, gridX]).'; 
Y =  reshape(violinMesh(:,2), [gridY, gridX]).'; 
Z =  reshape(violinMesh(:,3), [gridY, gridX]).';

%{
figure(1)
surf(X,Y,Z);
title('Violin surface');
%zlim([0,100]);
%}

violinInfos = {X,Y,Z, violinMesh}; 
%% Velocity Fields

numFreqBins = length(table2array(csvVel(1,:)))-3;
velsMatrix = table2array( csvVel(:,4:numFreqBins+3) );

velsMatrix( velsMatrix == 0) = NaN; % clean the violin border to accident zeros

velocityFields = cell(numFreqBins,1);

for ii = 1:numFreqBins
    velocityFields{ii} = reshape(velsMatrix(:,ii), [gridY, gridX]).';  
    velocityFields{ii}(27,2) = NaN; % cleaning undesired velocities!
    velocityFields{ii}(27,15) = NaN;
end

%{
figure(2)
surf(X,Y,abs(velocityFields{1}));
title('Velocity field surface - f1');
%}

%% Reconstruct the hologram geometry

gridX = length(unique(hologramMesh(:,1))); % useless?
gridY = length(unique(hologramMesh(:,2))); % useless?

% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> 
% ex. try  reshape(1:10,[5,2])  what we actuallly want is  reshape(1:10,[2,5])'

X =  reshape(hologramMesh(:,1), [gridY, gridX]).'; % useless?
Y =  reshape(hologramMesh(:,2), [gridY, gridX]).'; % useless?
Z =  reshape(hologramMesh(:,3), [gridY, gridX]).';  % rename?

%{
figure(3)
surf(X,Y,Z);
title('Hologram surface');
%}

hologramInfos = {X,Y,Z}; % rename Z ad Zh?
%% Pressure Fields 

numFreqBins = length(table2array(csvPress(1,:)))-3; % useless?
pressMatrix = table2array( csvPress(:,4:numFreqBins+3) );
pressureFields = cell(numFreqBins,1);

for ii = 1:numFreqBins
    pressureFields{ii} = reshape(pressMatrix(:,ii), [gridY, gridX]).';   
end

for ii = 1:length(pressureFields)
    pressureFields{ii} = downsampling(pressureFields{ii}, 8, 8);
end

for ii = 1:length(hologramInfos) % useful to plot the reconstructed pressure.
    hologramInfos{ii} = downsampling(hologramInfos{ii}, 8, 8);
    
end

% Show first frequency pressure field
%{
figure(500)
surf(hologramInfos{1}, hologramInfos{2}, abs(pressureFields{1}));
title('actual pressure')
%}

hologramMesh =  [reshape(hologramInfos{1}', [64, 1]),...
                 reshape(hologramInfos{2}', [64, 1]),...
                 reshape(hologramInfos{3}, [64, 1])];
             
hologramInfos = {hologramInfos{1},hologramInfos{2},hologramInfos{3},hologramMesh};


end

