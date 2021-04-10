function [outMatrix] = downsampling(inputMatrix, nrows, ncols)

rowsfactor = round(length(inputMatrix(:, 1))/nrows);
colsfactor = round(length(inputMatrix(1, :))/ncols);

outMatrix = inputMatrix(1:(rowsfactor):end, 1:(colsfactor):end);
end

