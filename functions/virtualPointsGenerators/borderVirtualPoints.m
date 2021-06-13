function [border] = borderVirtualPoints(pts, xBorder, yBorder, zVal)
%FORMpts Summary of this function goes here
%   Detailed explanation goes here
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
                        xlimDown = max(1,round(xBorder*xLen));
                        xlimUp = min(round(abs((1-xBorder)*xLen)),xLen);
                        
                        if (jj > xlimDown && jj < xlimUp) %& mod(jj,cut) == 0
                          z(idxs(jj)) = nan;
                        else
                        end
            end
        end
        
        for ii = 1:length(xUnique)
            idxs = find(x==xUnique(ii));
            idxs = intersect(idxs, notNanIdxs);
            yCheck = y(idxs);
            yLen = length(yCheck);
            for jj = 1:yLen
                ylimDown = max(1,round(yBorder*yLen));
                ylimUp = min(round((1-yBorder)*yLen),yLen);
                        if (jj > ylimDown && jj < ylimUp) 
                            z(idxs(jj)) = nan;
                        else
                        end
            end
        end
        
    border = [x(:), y(:), z(:)];
    figure()
    plot3(border(:,1), border(:,2), border(:,3),'.', 'markerSize', 10);
%     hold on
%     plot3(border(:,1), border(:,2), -0.01*ones(size(z)),'.');
end