function [virtualPoints] = formVirtualPoints(virtualPoints, scale,  xCut,yCut,xBorder, yBorder, active)
%FORMVIRTUALPOINTS Summary of this function goes here
%   Detailed explanation goes here

    x = virtualPoints(:,1);
    y = virtualPoints(:,2);
    z = virtualPoints(:,3);
    
    xUnique = unique(x);
    yUnique = unique(y);
    
    notNanIdxs = find(~isnan(z));
%     xMin = min(x(notNanIdxs));
     xMax = max(x(notNanIdxs));
%     yMin = min(y(notNanIdxs));
    yMax = max(x(notNanIdxs));
    zVal = mean(z(notNanIdxs));
% 
%     
%     idxsX = intersect(find(x>xMin),find( x<xMax));
% %     idxsY = intersect(find(y>yMin),find( y<yMax));
%       z(unique([idxsX; idxsY])) = zVal;
%     
    figure()
    hold on;
    if active
        
        for ii = 1:length(xUnique)
            idxs = find(x==xUnique(ii));
            idxs = intersect(idxs, notNanIdxs);
            yCheck = y(idxs);

            
            for jj = 1:length(yCheck)
                if abs(yCheck(jj))< 0.2*yMax
                       z(idxs(jj)) = nan;
                 end 
                if (jj >= yBorder+1 && jj <= length(yCheck) -yBorder) && abs(yCheck(jj))< 0.15
                    
                    if mod(jj,yCut) == 0 
                        z(idxs(jj)) = nan;
                    end
                else 
                end
                    
                plot3(x(idxs(jj)), y(idxs(jj)),z(idxs(jj)), '.');
            end
            

        end

        for ii = 1:length(yUnique)
            idxs = find(y==yUnique(ii));
            notNanIdxs = find(~isnan(z));
            idxs = intersect(idxs, notNanIdxs);
            xCheck = x(idxs);
            
            for jj = 1:length(xCheck)
%                 if abs(xCheck(jj))< 0.15*xMax
%                     z(idxs(jj)) = nan;
%                 end 
                    if jj >= xBorder+1 && jj <= length(xCheck) - xBorder
                        if mod(jj,xCut) == 0 
      
                        else
                            z(idxs(jj)) = nan;
                        end
                    end
            end
        end
    end
    virtualPoints = [scale*x(:), scale*y(:), z(:)];
end

