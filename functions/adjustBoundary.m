function [finalPoints] = adjustBoundary(violinInfos, plotData)
%ADJUSTBOUNDARY Summary of this function goes here
%   Detailed explanation goes here
plotData = false;
violinInfos{4} = violinInfos{4}*1000;
violinPointsNaN =  violinInfos{4};

violinPoints = violinInfos{4};
index  = find(isnan(violinPoints(:,3)));


violinPoints(index,:) = [];
x = violinPoints(:,1);
y = violinPoints(:,2);
z = violinPoints(:,3);

xUnique = unique(violinInfos{4}(:,1));
yUnique = unique(violinInfos{4}(:,2));



boundaryIndex = boundary(x,y);
pointsToAdd = [];
edge = 10;
polyDegree = 3;
for ii = 1:length(yUnique)
    idxs = find(y==yUnique(ii));
    if idxs~[]
    
    
    xCheck = x(idxs);
    xStep = xUnique(2)-xUnique(1);
    
    points2Check = [];
    for jj = 2:length(xCheck)-1
        if xCheck(jj) -xCheck(jj-1) > 3* xStep
            points2Check = [xCheck(jj), xCheck(jj-1)];
            ind = jj;
        end
    end
    
    if ~isempty(points2Check)
        idxsCLeft = idxs(ind-2:ind-1);
        xCLeft = x(idxsCLeft);
    
        xCL = (max(xCLeft)+xStep);
        p = polyfit(xCLeft,z(idxsCLeft),polyDegree);
        zCL = polyval(p,xCL);
        
        if zCL>z(idxsLeft(1))
            zCL = z(idxsLeft(2)) -((z(idxsLeft(1))-z(idxsLeft(2)))/2);
        end
        
        idxsCRight = idxs(ind:ind+1);
        xCRight = x(idxsCRight);
    
        xCR = (min(xCRight)-xStep);
        p = polyfit(xCRight,z(idxsCRight),polyDegree);
        zCR = polyval(p,xCR);
        if zCR>z(idxsCRight())
            zCR = z(idxsCRight(1)) - ((z(idxsCRight(2))-z(idxsCRight(1)))/2);
        end
    end
    
    idxsLeft = idxs(1:edge);
    xLeft = x(idxsLeft);
    
    xL = (min(xLeft)-xStep);
    p = polyfit(xLeft,z(idxsLeft),polyDegree);
    zL = polyval(p,xL);
    if zL>z(idxsLeft(1))
       zL = z(idxsLeft(1)) -((z(idxsLeft(2))-z(idxsLeft(1)))/2);
    end
    
    idxsRight = idxs(length(idxs)-edge+1:end);
    xRight = x(idxsRight);
    xR = max(xRight)+xStep; 
    p = polyfit(xRight,z(idxsRight),polyDegree);
    zR = polyval(p,xR);

    if zR>z(idxsRight(end))
       zR = z(idxsRight(end)) - ((z(idxsRight(end-1))-z(idxsRight(end)))/2);
    end
    
    pointsToAdd = [pointsToAdd; xL, yUnique(ii), zL; xR, yUnique(ii), zR];
    
    if ~isempty(points2Check)
        pointsToAdd = [pointsToAdd; xCL, yUnique(ii), zCL; xCR, yUnique(ii), zCR];
    end
    
    
    
    
    
    if plotData
        figure()
        plot(x(idxs),z(idxs),'.')
        hold on;
        plot(xL,zL,'*',xR,zR,'*');
        hold on;
        if ~isempty(points2Check)
            plot(xCL,zCL,'*',xCR,zCR,'*');
        end
    end   
    end
end


points2Add = [];
polyDegree = 1;

for ii = 1:length(xUnique)
    idxs = find(x==xUnique(ii));
    if idxs~[]
    edge = 2;
    
    yCheck = y(idxs);
    yStep = yUnique(2)-yUnique(1);
    
    idxsLeft = idxs(1:edge);
    yLeft = y(idxsLeft);
    
    yL = (min(yLeft)-yStep);
    p = polyfit(yLeft,z(idxsLeft),polyDegree);
    zL = polyval(p,yL);
    if zL < 0
        zL =z(idxsLeft(1))/2; 
    end  
%     if zL>z(idxsLeft(1))
%        zL = z(idxsLeft(1)) -((z(idxsLeft(2))-z(idxsLeft(1)))/2);
%     end
    
    idxsRight = idxs(length(idxs)-edge+1:end);
    yRight = y(idxsRight);
    yR = max(yRight)+yStep; 
    p = polyfit(yRight,z(idxsRight),polyDegree);
    zR = polyval(p,yR);
    if zR < 0
        zR = z(idxsRight(end-1))/2; 
    end    
%     if zR>z(idxsRight(end))
%        zR = z(idxsRight(end)) - ((z(idxsRight(end-1))-z(idxsRight(end)))/2);
%     end
      add = true;
      
      for jj = 1:length(pointsToAdd(:,1))
         if pointsToAdd(jj,1) == [xUnique(ii)]
             %if pointsToAdd(jj,2) == yL || pointsToAdd(jj,2) == yR
                add = false;
             %end
         end
      end
      
      notadd = find(round(pointsToAdd(:,1:2)) == round([xUnique(ii),yL]));
      notadd = intersect(notadd, find(round(pointsToAdd(:,1:2)) == round([xUnique(ii),yR])));
      if isempty(notadd)
        points2Add = [points2Add; xUnique(ii), yL, zL; xUnique(ii), yR, zR];
      end
%     if ~isempty(points2Check)
%         pointsToAdd = [pointsToAdd; xCL, yUnique(ii), zCL; xCR, yUnique(ii), zCR];
%     end
    
    
    
    
    if plotData
        figure()
        plot(y(idxs),z(idxs),'.')
        hold on;
        plot(yL,zL,'*',yR,zR,'*');
        hold on;
%         if ~isempty(points2Check)
%             plot(xCL,zCL,'*',xCR,zCR,'*');
%         end
    end   
    end
end


figure; plot3(x,y,z,'.','MarkerSize',10);xlabel('x');
hold on; plot3(pointsToAdd(:,1),pointsToAdd(:,2),pointsToAdd(:,3),'.','MarkerSize',10)
hold on; plot3(points2Add(:,1),points2Add(:,2),points2Add(:,3),'.','MarkerSize',5)

xBorder = x(boundaryIndex);
yBorder = y(boundaryIndex);
zBorder = z(boundaryIndex);

pointsToAdd = [pointsToAdd; points2Add];
notNaN = find(~isnan(violinPointsNaN));
replaceIdxs = zeros(size(pointsToAdd(:,1)));

for ii = 1:length(pointsToAdd)
    replace = find(round(violinPointsNaN(:,1)) == round(pointsToAdd(ii,1)));
    replaceIdxs(ii) = intersect(replace,find(round(violinPointsNaN(:,2)) == round(pointsToAdd(ii,2))));

end

violinPointsNaN = violinInfos{4};
violinPointsNaN(replaceIdxs,3) = pointsToAdd(:,3);

index = find(isnan(violinPointsNaN(:,3)));
violinPointsNaN(index,3) = 0;

finalPoints = sortrows(violinPointsNaN,1)/1000;
figure; plot3(pointsToAdd(:,1),pointsToAdd(:,2),pointsToAdd(:,3),'.','MarkerSize',10)
hold on; plot3(finalPoints(:,1),finalPoints(:,2),finalPoints(:,3),'.','MarkerSize',5)

  
end

