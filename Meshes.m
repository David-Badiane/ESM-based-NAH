close all
clearvars
clc
addpath 'CSV\pressure' 
addpath 'CSV\velocity' 
addpath 'functions'


%% read the csv
csvVel = readtable("vel_1.csv");
meshcord = table2array(csvVel(:,1:3));
%meshcord(:, :, meshcord == NaN) = 0;

%% extract x,y axis values without repetitions 

x = unique(meshcord(:,1));
y = unique(meshcord(:, 2));
gridX = length(x);
gridY = length(y);

%% reconconstruct the plate geometry
% reshape into size (16*64) then transpose. because reshape orders the
% vector by columns --> ex. try  reshape(1:10,[5,2]) we want
% reshape(1:10,[2,5])'   
z = meshcord(:,3);
platemesh = reshape(z, [gridY, gridX]).';   
figure 
mesh(platemesh);

%% reconstruct the field
csvvelmatrix = zeros(size(csvvel));

% insert double in the zero matrix
csvvelmatrix(:,1:3) = table2array(csvVel(:,1:3));


% convert cell to double
for i = 1 : 1024
    for j = 4 : 49
        csvvelmatrix(i, j) = str2double(cell2mat(table2array(csvvel(i, j))));
    end
end

%% velocity field 
selfreq = 3 + 1;
meshfreq = [csvvelmatrix(:,1), csvvelmatrix(:,2), csvvelmatrix(:,selfreq)];

velocityfield = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        b = meshfreq((16*(i-1)) + j, :);
        velocityfield(i, j) = abs(b(3));
    end
end        

figure 
mesh(velocityfield);

figure
%surf(velocityfield);

%% pressure field

csvpress = readtable("acpr_1.csv");

% reconstruct the field
csvpressmatrix = zeros(size(csvpress));

% insert double in the zero matrix
for i = 1 : 1024
    for j = 1 : 3
        csvpressmatrix(i, j) = table2array(csvpress(i, j));
    end
end

% convert cell to double
for i = 1 : 1024
    for j = 4 : 20 % TO DO, there is a data error at index 21
        csvpressmatrix(i, j) = str2double(cell2mat(table2array(csvpress(i, j))));
    end
end

% field
selfreq = 3 + 17;
meshpress = [csvpressmatrix(:,1), csvpressmatrix(:,2), csvpressmatrix(:,selfreq)];

pressfield = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        b = meshpress((16*(i-1)) + j, :);
        pressfield(i, j) = abs(b(3));
    end
end        

figure 
mesh(pressfield);

% downsampling

pressfielddown = downsampling(pressfield, 8, 8);

figure
mesh(pressfielddown);

figure
surf(pressfielddown);
