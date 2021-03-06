function [violinInfos, velocityFields, hologramInfos, pressureFields, eigenFreqz] = importData(velocityFileName, pressureFileName, resampleX, resampleY, plotData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPORTDATA this function read the csv pressure and normal velocity       %%%%%%%%%%%%
%%% extract their data and sotre the in a cell arrays                        %%%%%%%%%%%%
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

pX = length(unique(violinMesh(:, 1)));
pY = length(unique(violinMesh(:, 2)));

% clean the violin border to accident zeros
violinMeshZ = violinMesh(:, 3);
violinMeshZ( violinMeshZ <= 0) = NaN; % TO DO asl if they are an error
violinMesh(:,3) = violinMeshZ;

% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> 
% ex. try  reshape(1:10,[5,2])  what we actuallly want is  reshape(1:10,[2,5])'

X =  reshape(violinMesh(:,1), [pY, pX]).'; 
Y =  reshape(violinMesh(:,2), [pY, pX]).'; 
Z =  reshape(violinMesh(:,3), [pY, pX]).';

if plotData
figure(1)
surf(X,Y,Z);
title('Violin surface');
zlim([0,100]);
end

violinInfos = {X,Y,Z, violinMesh}; 
%% Velocity Fields

numFreqBins = length(table2array(csvVel(1,:)))-3;
velsMatrix = table2array( csvVel(:,4:numFreqBins+3) );
velocityFields = cell(numFreqBins,1);

xAxis = unique(violinMesh(:,1));
yAxis = unique(violinMesh(:,2));

for ii = 1:numFreqBins
    idxWrongPts = intersect ( find( violinMesh(:,2) == yAxis(3) | violinMesh(:,2) == yAxis(end-2)),...
        find(violinMesh(:,1) == xAxis(end-18))); 
    idxGoodPts = intersect ( find( violinMesh(:,2) == yAxis(3) | violinMesh(:,2) == yAxis(end-2)),...
        find(violinMesh(:,1) == xAxis(end-17)));

    vel = velsMatrix(:,ii);
    vel(isnan(violinMesh(:,3))) = nan;
    vel(idxWrongPts) = vel(idxGoodPts);
    velocityFields{ii} = reshape(vel, [pY, pX]).';
    if plotData
        figure(2)
        surf(X,Y,abs(velocityFields{ii}));
        title('Velocity field surface - f1');

        figure(3)
        plot3(violinMesh(:,1), violinMesh(:,2), abs(vel),'.');
        pause(0.5);
    end
end

%% Reconstruct the hologram geometry

pX = length(unique(hologramMesh(:,1))); % useless?
pY = length(unique(hologramMesh(:,2))); % useless?

X =  reshape(hologramMesh(:,1), [pY, pX]).'; % useless?
Y =  reshape(hologramMesh(:,2), [pY, pX]).'; % useless?
Z =  reshape(hologramMesh(:,3), [pY, pX]).';  % rename?

if plotData
    figure(3)
    surf(X,Y,Z);
    title('Hologram surface');
end

hologramInfos = {X,Y,Z}; % rename Z ad Zh?
%% Pressure Fields 

numFreqBins = length(table2array(csvPress(1,:)))-3; % useless?
pressMatrix = table2array( csvPress(:,4:numFreqBins+3) );
pressureFields = cell(numFreqBins,1);

for ii = 1:numFreqBins
    pressureFields{ii} = reshape(pressMatrix(:,ii), [pY, pX]).';   
end

for ii = 1:length(hologramInfos) % useful to plot the reconstructed pressure.
    hologramInfos{ii} = downsampling(hologramInfos{ii}, resampleX, resampleY);    
end

hologramMesh =  [reshape(hologramInfos{1}, [resampleX*resampleY, 1]),...
                 reshape(hologramInfos{2}, [resampleX*resampleY, 1]),...
                 reshape(hologramInfos{3}, [resampleX*resampleY, 1])];
          

for ii = 1:length(pressureFields)
    pressureFields{ii} = downsampling(pressureFields{ii}, resampleX, resampleY);
    p = reshape(pressureFields{ii}, [resampleX*resampleY,1]);
    press = pressMatrix(:,ii);
    if plotData
        figure(5)
        plot3(hologramMesh(:,1), hologramMesh(:,2), abs(p),'.');
        hold on
        oldHMesh = table2array(csvPress(:,1:3));
        plot3(oldHMesh(:,1), oldHMesh(:,2), abs(press), 'o');
        hold off

        figure(6)
        surf(hologramInfos{1}, hologramInfos{2}, abs(pressureFields{ii}))
    end
end

   
hologramInfos = {hologramInfos{1},hologramInfos{2},hologramInfos{3},hologramMesh};
if plotData
    figure(500)
    surf(hologramInfos{1}, hologramInfos{2}, abs(pressureFields{1}));
    title('actual pressure')
end
end

