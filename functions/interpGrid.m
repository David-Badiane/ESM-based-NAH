function [zCordsinter] = interpGrid(meshPoints, xData, yData, pX, pY, plotData)
%INTERPGRID
%   Detailed explanation goes here

xCords = reshape(meshPoints(:,1), [pX,pY]);
yCords = reshape(meshPoints(:,2), [pX,pY]);
zCords = reshape(meshPoints(:,3), [pX,pY]);

xData = 0.001.*xData;
yData = 0.001.*yData;

zCordsinter = griddata(xCords,yCords,zCords,xData,yData);

if plotData
figure()
plot3(xData,yData,zCordsinter, '.', 'markerSize', 20);
hold on 
plot3(meshPoints(:,1),meshPoints(:,2),meshPoints(:,3), '.');

end

