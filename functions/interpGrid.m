
function [zCordsinter] = interpGrid(meshPoints, xData, yData, pX, pY, plotData)
% INTERPGRID this function interpolate one grid with a couple of
% coordinates

%   INPUTS
%   meshPoints  (2Darray) = points of the grid   
%   xData       (array)   = value of x for the interpolation;
%   yData       (array)   = value of y for the interpolation;
%   pX          (array)   = x coordinates of the mesh;
%   pY          (array)   = y coordinates of the mesh;
%   plotData    (boolean) = choose if the plot have to be shown;

%   OUPUTS
%   zCordInter  (array)   = z value of the coordinates interpolated;

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

