function [intPts] = innerVirtualPoints(pts, xBorder, yBorder, zVal)
%innerVirtualPoints form a grid with the inner points of the geometry
% INPUTS
% pts  [nPts x 3] (2DArray) = points of the mesh
% xBorder         (double)  = how many points interior to the border along X
% yBorder         (double)  = how many points interior to the border along y
% zVal            (double)  = value of the z coordinate of the grid
% OUTPUTS
% intPts [nPts x 3](2DArray) = grid with the inner points of the geometry

close all

    x = pts(:,1);
    y = pts(:,2);
    z = pts(:,3);
    
    xUnique = unique(x);
    yUnique = unique(y);
    
    notNanIdxs = find(~isnan(z));
    z(notNanIdxs) = zVal;
%     idxsX = intersect(find(x>xMin),find( x<xMax));
% %     idxsY = intersect(find(y>yMin),find( y<yMax));
%       z(unique([idxsX; idxsY])) = zVal;
%     
    figure()
    hold on;

        for ii = 1:length(yUnique)
            idxs = find(y==yUnique(ii));
            idxs = intersect(idxs, notNanIdxs);
            xCheck = x(idxs);
            xLen = length(xCheck);
            for jj = 1:xLen
                        xlimDown = round(xBorder*xLen);
                        xlimUp = round(abs((1-xBorder)*xLen));
                        
                        if (jj >= xlimDown && jj <= xlimUp)
                        else
                            z(idxs(jj)) = nan;
                        end
            end
        end
        
        for ii = 1:length(xUnique)
            idxs = find(x==xUnique(ii));
            idxs = intersect(idxs, notNanIdxs);
            yCheck = y(idxs);
            yLen = length(yCheck);
            for jj = 1:yLen
                ylimDown = round(yBorder*yLen);
                ylimUp = round((1-yBorder)*yLen);
                        if (jj >= ylimDown && jj <= ylimUp) 
      
                        else
                            z(idxs(jj)) = nan;
                        end
            end
        end
        
    intPts = [x(:), y(:), z(:)];
    figure()
    plot3(intPts(:,1), intPts(:,2), intPts(:,3),'.', 'markerSize', 10);
    hold on
    plot3(intPts(:,1), intPts(:,2), -0.01*ones(size(z)),'.');
end