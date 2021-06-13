function [rectPts] = rectVirtPoints(pts, xRect, yRect, xCenter, yCenter, zVal)
%RECTVIRTPOINTS Summary of this function goes here
%   Detailed explanation goes here
    x = pts(:,1); y = pts(:,2); z = pts(:,3);
    z(:) = nan;
    
    xUp = xCenter+round(xRect/2);
    xDown = xCenter-round(xRect/2);
    
    yUp = yCenter+round(yRect/2);
    yDown = yCenter-round(yRect/2);
    
    indexes = intersect(find(x>= xDown & x <= xUp), find(y>= yDown & y <= yUp));
    z(indexes) = zVal;
    rectPts = [x,y,z];
    figure(1)
    plot3(x,y,z, '.');
    hold on
    plot3(pts(:,1),pts(:,2),pts(:,3), '.')
end

