
function [zCordsinter] = interpGrid(meshPoints, xData, yData, pX, pY, plotData)
% INTERPGRID this function interpolates in the surface defined by meshPoints 
%            over the points of coordinate xData and yData

%   INPUTS
%   meshPoints  (2Darray)  = points of the grid on which we interpolate   
%   xData       (2Darray)  = X matrix - target points for the interpolation;
%   yData       (2Darray)  = Y matrix - target points for the interpolation;
%   pX          (double)   = number x coordinates of the mesh;
%   pY          (double)   = number y coordinates of the mesh;
%   plotData    (boolean)  = choose if the plot have to be shown;

%   OUPUTS
%   zCordInter  (2Darray)  = Z matrix of the interpolated coordinates;

xCords = reshape(meshPoints(:,1), [pY,pX]).';
yCords = reshape(meshPoints(:,2), [pY,pX]).';
zCords = reshape(meshPoints(:,3), [pY,pX]).';


x = reshape(xData, [length(xData(1,:))*length(xData(:,1)),1]);
y = reshape(yData, [length(yData(1,:))*length(yData(:,1)),1]);

zCordsinter = griddata(xCords,yCords,zCords,xData.',yData.');
 
if plotData
figure(120)
plot3(xData,yData,zCordsinter, '.', 'markerSize', 20);
hold on 
plot3(meshPoints(:,1),meshPoints(:,2),meshPoints(:,3), '.');
xlabel(' x [m] ');
ylabel('z [m]');

end
end

