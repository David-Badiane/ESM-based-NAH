close all
clearvars
clc
        
%% read the csv
csvvel = readtable("vel_1.csv");
meshcord = [csvvel(:,1), csvvel(:,2), csvvel(:,3)];
meshcord = table2array(meshcord);
%meshcord(:, :, meshcord == NaN) = 0;

%% take the axes indices
x = meshcord(:, 1);
x = x(1:16:(64*16));

y = meshcord(:, 2);
y = y(1:16);

%% reconconstruct the plate geometry
platemesh = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        a = meshcord((16*(i-1)) + j, :);
        platemesh(i, j) = a(3);
    end
end        

figure 
mesh(platemesh);

%% reconstruct the field

csvvelmatrix = zeros(size(csvvel));

% insert double in the zero matrix
for i = 1 : 1024
    for j = 1 : 3
        csvvelmatrix(i, j) = table2array(csvvel(i, j));
    end
end

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

%figure
%surf(velocityfield);

%% pressure field

csvpress = readtable("C:\Users\Andrea\Downloads\violin_nah_example\violin_nah_example\pressure\acpr_1.csv");

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

%% Green's matrix
% acoustic hologram at point r in matrix form -> p(w) = i*w*G_p(w)*q(w)
% G_p(w) Green's matrix mxn (M points on the hologram, N point on the surface) 
%   calculated from virtual points a and r
%   g_w(r_m, a_n) = (1/4/pi)*exp(-1i* w/c * norm((r - a), 2)) / norm((r - a), 2)       
% q(w) equivalent source weights

c = 343; % sound speed m/s

presscord = meshcord;
presscord(:, 3) = csvpressmatrix(1, 3);

virtualcord = presscord;
virtualcord(:, 3) = csvvelmatrix(:, 3);
virtualcord(isnan(virtualcord)) = 0;

lattice = min(abs(presscord(:, 3) - virtualcord(:, 3)));

% points on the pressure field
r_m = presscord;
% virtual points
a_n = virtualcord - [0, 0, lattice];

% extract the eigenfrequencies
eigfreq = csvpress.Properties.VariableNames;
eigfreq = char(eigfreq(1, 4:end));
eigfreq = eigfreq(:,3:end);

eigfreqdouble = zeros(size(eigfreq(:,1)));
for i = 1:length(eigfreq(:,1)) % for cycle because strrep needs vectors
    eigfreq(i, :) = strrep(eigfreq(i, :),'_','.');
    eigfreqdouble(i, :) = str2double(eigfreq(i, :));
end

eigfreq = eigfreqdouble;

G_w = zeros(length(r_m(:, 1)), length(a_n(:, 1)), length(eigfreq(:, 1))); 
for i = 1 : length(r_m(:, 1))
    for j = 1 : length(a_n(:, 1))
        for k = 1 : length(eigfreq(:, 1))
        G_w(i, j, k) = (1/4/pi) * exp(-1i * eigfreq(k, 1)./c *... 
            0.001*norm((r_m(i, :) - a_n(j, :)), 2))...
            / 0.001*norm((r_m(i, :) - a_n(j, :)), 2);
        end
    end
end
