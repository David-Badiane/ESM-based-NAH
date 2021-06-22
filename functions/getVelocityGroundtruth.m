function [X,Y,surfV] = getVelocityGroundtruth(v_ex_vector, velocityFilename, figureNum)
% GETVELOCITYGROUNDTRUTH this function calculate the array velocity into
% surface

%   INPUTS
%   v_ex_vector      (array)   = vector of the velocities groung truth ;
%   velocityFilename (cell)    = names of the file to choose;
%   figureNum        (double)  = number of the figure to plot;

%   OUTPUTS
%   X                (array)   = x coordinates of the mesh;
%   Y                (array)   = y coordinates of the mesh;
%   surfV            (2Darray) = z value of the velocity;

violinMesh = table2array(readtable('grid75x35.csv')); 
velocityData = readmatrix([velocityFilename, '.csv']); 

pX = length(unique(violinMesh(:, 1)));
pY = length(unique(violinMesh(:, 2)));

X = reshape(violinMesh(:,1), [pY, pX]).';
Y = reshape(violinMesh(:,2), [pY, pX]).';
Z = reshape(violinMesh(:,3), [pY, pX]).';


F1 = scatteredInterpolant(velocityData(:,1),velocityData(:,2), abs(v_ex_vector));

vq1 = F1(X,Y);
vRows = length(vq1(:,1)); 
vCols = length(vq1(1,:));
zVel = reshape(vq1, [ vCols*vRows, 1 ]);
zVel(isnan(reshape(Z,[pX*pY,1]))) = nan;

surfV = reshape(zVel, [ vRows, vCols]);

figure(figureNum)
subplot 121
surf(X,Y,vq1); view(2);
xlabel('X   [m]');
ylabel('Y   [m]');


subplot 122
surf(X,Y,surfV); view(2);
xlabel('X   [m]');
ylabel('Y   [m]');

end

