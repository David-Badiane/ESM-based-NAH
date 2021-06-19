function [rectPts] = rectVirtPoints(pts, xRect, yRect, xCenter, yCenter, zVal)
%RECTVIRTPOINTS Summary of this function goes here
%   Detailed explanation goes here
    x = xRect/2*pts(:,1)/max(pts(:,1));
    y = yRect/2*pts(:,2)/max(pts(:,2));
    z = pts(:,3);
    z(:) = nan;
    
    xUp = xCenter+xRect/2;
    xDown = xCenter-xRect/2;
    
    yUp = yCenter+yRect/2;
    yDown = yCenter-yRect/2;
    
    indexes = intersect(find(x>= xDown & x <= xUp), find(y>= yDown & y <= yUp));
    z(indexes) = zVal;
    rectPts = [x,y,z];
    figure(1)
    plot3(pts(:,1),pts(:,2),pts(:,3), '.');
    hold on
    plot3(rectPts(:,1),rectPts(:,2),rectPts(:,3), '.')
end

