function [outMatrix] = downsampling(inputMatrix, nrows, ncols)
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

rowsfactor = round(length(inputMatrix(:, 1))/nrows);
colsfactor = round(length(inputMatrix(1, :))/ncols);

outMatrix = inputMatrix(1:(rowsfactor):end, 1:(colsfactor):end);
end

