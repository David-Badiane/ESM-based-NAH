function [ellipse] = ellipseVirtualPoints(pts, maxRx, maxRy,minRx, minRy,zVal)
%FORMVIRTUALPOINTS Summary of this function goes here
%   Detailed explanation goes here

    x = pts(:,1); y = pts(:,2); z = pts(:,3);
    maxX = max(x); maxY = max(y);
    
    xAxis = unique(x);
    yAxis = unique(y);
     
    ncols = length(xAxis);
    nrows = length(yAxis);
    
    maxRx = round(ncols*maxRx/maxX/2);
    maxRy = round(nrows*maxRy/maxX/2);
    minRx = round(ncols*minRx/maxX/2);
    minRy = round(nrows*minRy/maxX/2);
   
    notNanIdxs = find(~isnan(z));
    z(notNanIdxs) = nan;

    figure(3)
    plot(x(notNanIdxs),y(notNanIdxs),'.');
    hold on;
    

    [X Y] = meshgrid(1:nrows, 1:ncols);
    % Next create the circle in the image.
    maxX = max(x); minX = min(x);
    ratio = abs(minX)./(maxX + abs(minX));
    centerX = round(ratio*nrows);
    centerY = round(ncols/2);
    ellipse = zeros(nrows,ncols);

    
    ellipse(((Y-centerY).^2 ./ maxRy^2 + (X-centerX).^2 ./ maxRx^2 <= 1 ))=1;
    ellipse(((Y-centerY).^2 ./ minRy^2 + (X-centerX).^2 ./ minRx^2 <= 1 ))=0;
    %meshIdxs = zeros(size(X));
    
    ellipseIdxs = find(ellipse == 1);
    z(ellipseIdxs) = zVal;
    plot3(x,y,z,'.')
    hold off;
    ellipse = [x,y,z];
    end
 

