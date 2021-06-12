function [zCordsinter] = interpGrid(meshPoints, xData, yData, plotData)
%INTERPGRID
%   Detailed explanation goes here

xCords = reshape(meshPoints(:,1), [128,128]);
yCords = reshape(meshPoints(:,2), [128,128]);
zCords = reshape(meshPoints(:,3), [128,128]);

zCordsinter = griddata(xCords,yCords,zCords,xData,yData);

if plotData
figure()
plot3(xData,yData,zCordsinter, '.', 'markerSize', 20);
hold on 
plot3(meshPoints(:,1),meshPoints(:,2),meshPoints(:,3), '.');

end

