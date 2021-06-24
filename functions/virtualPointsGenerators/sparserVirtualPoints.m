function [outPts] = sparserVirtualPoints(pts,controller, xCut, yCut)
%SPARSERVIRTUALPOINTS Summary of this function goes here
% INPUTS
% pts        (2DArray)= points of the mesh
% controller (double) = 0, modular sparsing
% controller          = 1, random sparsing
% controller          = 2, modular random sparsing - modular weighted randomically
% xCut       (double) = for modular sparsing along x, takes 1 point out of xCut
% yCut       (double) = for modular sparsing along x, takes 1 point out of xCut
% OUTPUTS
% outPs      (2DArray)= out matrix with sparser points

    x = pts(:,1); y = pts(:,2); z = pts(:,3);    
    xAxis = unique(x);
    yAxis = unique(y);
    notNanIdxs = find(~isnan(z));


switch(controller)
    case 0 
        for ii = 1:length(xAxis)
            idxs = find(x==xAxis(ii));
            idxs = intersect(idxs, notNanIdxs);
            yCheck = y(idxs);
            
            for jj = 1:length(yCheck)                    
                    if mod(jj,yCut) == 0 
                        z(idxs(jj)) = nan;
                    end
            end
        end
            
        for ii = 1:length(yAxis)
            idxs = find(y==yAxis(ii));
            idxs = intersect(idxs, notNanIdxs);
            xCheck = x(idxs);
            
            for jj = 1:length(xCheck)
                if mod(jj,xCut) == 0 

                else
                    z(idxs(jj)) = nan;
                end
            end
        end
    case 1
        
        deactivateIdxs = find(randi([0 1],1,length(notNanIdxs)) == 1);
        z(notNanIdxs(deactivateIdxs)) = nan;
    case 2
        
        for ii = 1:length(xAxis)
            idxs = find(x==xAxis(ii));
            idxs = intersect(idxs, notNanIdxs);
            yCheck = y(idxs);
            
            for jj = 1:length(yCheck)                    
                    if mod(jj,round(rand(1)*yCut) )== 0 
                        z(idxs(jj)) = nan;
                    end
            end
        end
            
        for ii = 1:length(yAxis)
            idxs = find(y==yAxis(ii));
            idxs = intersect(idxs, notNanIdxs);
            xCheck = x(idxs);
            
            for jj = 1:length(xCheck)
                if mod(jj,round(rand(1)*xCut)) == 0 

                else
                    z(idxs(jj)) = nan;
                end
            end
        end
end

outPts =[x,y,z];
figure()
plot3(x,y,z, '.')
        
end

