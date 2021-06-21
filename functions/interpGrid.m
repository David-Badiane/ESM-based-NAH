
function [zCordsinter] = interpGrid(meshPoints, xData, yData, pX, pY, plotData)
%INTERPGRID
%   Detailed explanation goes here

xCords = reshape(meshPoints(:,1), [pY,pX]).';
yCords = reshape(meshPoints(:,2), [pY,pX]).';
zCords = reshape(meshPoints(:,3), [pY,pX]).';


x = reshape(xData, [length(xData(1,:))*length(xData(:,1)),1]);
y = reshape(yData, [length(yData(1,:))*length(yData(:,1)),1]);

zCordsinter = griddata(xCords,yCords,zCords,xData.',yData.');
 
if plotData
figure()
plot3(xData,yData,zCordsinter, '.', 'markerSize', 20);
hold on 
plot3(meshPoints(:,1),meshPoints(:,2),meshPoints(:,3), '.');

end

