function [outMatrix] = downsampling_regular(inputMatrix, nrows, ncols, fileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION PERFORMS DOWNSAMPLING ON A MATRIX OR ARRAY               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                 %%%%%%%%%%%%
%%% inputMatrix = matrix to downsample (2DArray)                           %%%%%%%%%%%%
%%% nrows = target number of rows (double)                                 %%%%%%%%%%%%
%%% ncols = target number of columns (double)                              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% outMatrix = resampled matrix or array (2DArray)                        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = inputMatrix(:,1);
y = inputMatrix(:,2);
z = inputMatrix(:,3);

nanIdx = find(~isnan(z));
% offsetX = 603.72;
% offsetY = 353;
% 
% x = x + offsetX;
% y = y + offsetY;

maxX = max(x(nanIdx)); minX = min(x(nanIdx));
maxY = max(y(nanIdx)); minY = min(y(nanIdx));

deltaX = (maxX - minX)/ncols;
deltaY = (maxY - minY)/nrows;

pts = [];

% for ii = 1:ncols
%     diffX = abs (x + deltaX/2) - (maxX - (ii-1)*deltaX);
%     idxs = find(diffX < 100);
%     yy = [];
%     zz = [];
%     for jj = 1:nrows
%         [fval, flocs] = min(abs((y(idxs)-deltaY/2) - (maxY - (jj-1)*deltaY)));
%         yy = [yy; y(idxs(flocs))];
%         zz = [zz; z(idxs(flocs))];
%     end
%     [fval, floc] = min(diffX);
%     xx = ones(size(yy))*x(floc);
%     pts = [pts; xx, yy, zz];
% end

xRect = linspace(minX,maxX,ncols);
yRect = linspace(minY,maxY,nrows);
[X,Y] = meshgrid(xRect, yRect);
xRect = reshape(X, length(X(:,1))*length(X(1,:)) ,1);
yRect = reshape(Y, length(Y(:,1))*length(Y(1,:)) ,1);

pts = [];
for ii = 1:nrows*ncols
   [fval, floc] = min(sqrt((x - xRect(ii)).^2 + (y - yRect(ii)).^2));
   if fval > 10* min([deltaX ,deltaY])
       pts = [pts; xRect(ii) yRect(ii) NaN ];
   else
        pts = [pts; xRect(ii) yRect(ii) z(floc)];
   end
   
end


figure()
plot3(pts(:,1), pts(:,2), pts(:,3),'.');
hold on 
plot3(xRect, yRect, zeros(size(xRect)), '.', 'markerSize', 1);




outMatrix = pts;
writeMat2File(outMatrix, [fileName,'.csv'], {'x' 'y' 'z'}, 3, true );
end

