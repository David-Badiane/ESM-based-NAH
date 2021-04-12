close all
clearvars
clc
addpath 'CSV\pressure' 
addpath 'CSV\velocity' 
addpath 'functions'

%% read the csv

velocityFileName = 'vel_1.csv'; 
pressureFileName = 'acpr_1.csv';
[violinInfos, velocityFields, hologramInfos, pressureFields] = importData(velocityFileName, pressureFileName);

pressFieldsDown = cell(size(pressureFields));

% downsampling the simulated data in order to obtain a 8x8 grid

for ii = 1:length(pressureFields)
    pressFieldsDown{ii} = downsampling(pressureFields{ii}, 8, 8);
end

for ii = 1:length(hologramInfos)
    hologramInfosDown{ii} = downsampling(hologramInfos{ii}, 8, 8);
end

figure(5)
surf(hologramInfosDown{1}, hologramInfosDown{2}, abs(pressFieldsDown{1}));

%% Green's functions matrix



% 
% 
% 
% csvVel = readtable("vel_1.csv");
% meshcord = table2array(csvVel(:,1:3));
% %meshcord(:, :, meshcord == NaN) = 0;
% 
% %% extract x,y axis values without repetitions 
% 
% x = unique(meshcord(:,1));
% y = unique(meshcord(:, 2));
% gridX = length(x);
% gridY = length(y);
% 
% %% reconconstruct the plate geometry
% % reshape into size (16*64) then transpose. because reshape orders the
% % vector by columns --> ex. try  reshape(1:10,[5,2]) we want
% % reshape(1:10,[2,5])'   
% z = meshcord(:,3);
% platemesh = reshape(z, [gridY, gridX]).';   
% figure 
% mesh(platemesh);
% 
% %% reconstruct the field
% csvvelmatrix = table2array(csvVel(:,1:49));
% 
% selfreq = 3 + 1;
% meshfreq = [csvvelmatrix(:,1), csvvelmatrix(:,2), csvvelmatrix(:,selfreq)];
% 
% platemesh = reshape(z, [gridY, gridX]).';   
% 
% velocityfield = zeros(length(x), length(y));
% for i = 1:length(x)
%     for j = 1:length(y)
%         b = meshfreq((16*(i-1)) + j, :);
%         velocityfield(i, j) = abs(b(3));
%     end
% end        
% 
% figure 
% mesh(velocityfield);
% 
% figure
% %surf(velocityfield);
% 
% %% Pressure field
% csvpress = readtable("acpr_1.csv")
% 
% % Reconstruct the field
% csvpressmatrix = zeros(size(csvpress));
% 
% % Insert double in the zero matrix
% for i = 1 : 1024
%     for j = 1 : 3
%         csvpressmatrix(i, j) = table2array(csvpress(i, j));
%     end
% end
% 
% % Convert cell to double
% for i = 1 : 1024
%     for j = 4 : 20 % TO DO, there is a data error at index 21
%         csvpressmatrix(i, j) =table2array(csvpress(i, j));
%     end
% end
% 
% % Field
% selfreq = 3 + 1;
% meshpress = [csvpressmatrix(:,1), csvpressmatrix(:,2), csvpressmatrix(:,selfreq)];
% 
% pressfield = zeros(length(x), length(y));
% for i = 1:length(x)
%     for j = 1:length(y)
%         b = meshpress((16*(i-1)) + j, :);
%         pressfield(i, j) = abs(b(3));
%     end
% end        
% 
% figure 
% surf(pressfield);
% 

