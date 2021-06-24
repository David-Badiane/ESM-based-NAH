function [rectPts] = rectVirtPoints(pts, xRect, yRect, xCenter, yCenter, zVal)
%RECTVIRTPOINTS creates a rectangular grid

% INPUTS
% pts  [nPts x 3] (2DArray) = points of the mesh
% xRect           (double)  = edgeX value
% yRect           (double)  = edgeY value
% xCenter         (double)  = x coordinate of the center
% y Center        (double)  = y coordinate of the center
% zVal            (double)  = value of the z coordinate of the grid

% OUTPUTS
% rectPts [nPts x 3](2DArray) = rectangluar grid points
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

