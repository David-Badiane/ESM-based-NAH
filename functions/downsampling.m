function [outMatrix] = downsampling(inputMatrix, nrows, ncols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOWNSAMPLING this fucntion performs donwnsapling on a matrix or array  %%%%%%%%%%%%
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

[X, Y] = meshgrid(1:length(inputMatrix(1, :)), 1:length(inputMatrix(:, 1)));
[xq,yq] = meshgrid(linspace(1,length(inputMatrix(1, :)),ncols),...
    linspace(1,length(inputMatrix(:, 1)),nrows));
outMatrix = interp2(X, Y, inputMatrix, xq, yq);

end

