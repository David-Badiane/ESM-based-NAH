%function [outputArg1,outputArg2] = adjustBoundary(violinInfos, plotData)
%ADJUSTBOUNDARY Summary of this function goes here
%   Detailed explanation goes here
plotData = true;
violinPoints = violinInfos{4};
index  = find(isnan(violinPoints(:,3)));

violinPoints(index,:) = [];
x = violinPoints(:,1);
y = violinPoints(:,2);
z = violinPoints(:,3);

boundaryIndex = boundary(x,y);
if plotData
    figure()
    plot(x,y,'.')
    hold on;
    plot(x(boundaryIndex),y(boundaryIndex));
end

xBorder = x(boundaryIndex);
yBorder = y(boundaryIndex);
zBorder = z(boundaryIndex);


%end

