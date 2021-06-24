function [outMatrix] = downsampling_regular(inputMatrix, nrows, ncols, fileName, saveData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOWNSAMPLING this fucntion performs donwnsapling over a rectangular    %%%%%%%%%%%%
%%%              grid on a matrix                                          %%%%%%%%%%%%

%%% INPUTS                                                                 %%%%%%%%%%%%
%%% inputMatrix = matrix to downsample (2DArray)                           %%%%%%%%%%%%
%%% nrows    = target number of rows (double)                              %%%%%%%%%%%%
%%% ncols    = target number of columns (double)                           %%%%%%%%%%%%
%%% fileName = fileName for saving [only the name, not .csv] (double)      %%%%%%%%%%%%
%%% saveData = true if you want to save data on .csv file (boolean)        %%%%%%%%%%%%

%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% outMatrix = resampled matrix or array (2DArray)                        %%%%%%%%%%%%

x = inputMatrix(:,1);
y = inputMatrix(:,2);
z = inputMatrix(:,3);

nanIdx = find(~isnan(abs(z)));
maxX = max(x(nanIdx)); minX = min(x(nanIdx));
maxY = max(y(nanIdx)); minY = min(y(nanIdx));

deltaX = (maxX - minX)/ncols;
deltaY = (maxY - minY)/nrows;


xRect = linspace(minX,maxX,ncols);
yRect = linspace(minY,maxY,nrows);
[X,Y] = meshgrid(xRect, yRect);
xRect = reshape(X, length(X(:,1))*length(X(1,:)) ,1);
yRect = reshape(Y, length(Y(:,1))*length(Y(1,:)) ,1);

pts = [];
for ii = 1:nrows*ncols
   [fval, floc] = min(sqrt((x - xRect(ii)).^2 + (y - yRect(ii)).^2));
   if fval >  2* min([deltaX ,deltaY])
       pts = [pts; xRect(ii) yRect(ii) NaN ];
   else
       pts = [pts; xRect(ii) yRect(ii) z(floc)];
   end
   
end


figure(100)
plot3(pts(:,1), pts(:,2), abs(pts(:,3)),'.');
hold on 
plot3(xRect, yRect, zeros(size(xRect)), '.', 'markerSize', 0.5);
plot3(x,y,abs(z), 'x');
pause(0.1);
hold off;



outMatrix = pts;
if saveData
    writeMat2File(outMatrix, [fileName,'.csv'], {'x' 'y' 'z'}, 3, true );
end
end

