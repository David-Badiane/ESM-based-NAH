function [nanVel] = addNans(points, velocity)
%ADDNANS this function creates an array of size [length(points(:,3),1]
% where the output array has the same nan indexes of points(:,3)
% and its not nan indexes are filled by the value of velocity

%   INPUTS
%   points        (2Darray) = matrix of points whose nans mask is taken;
%   velocity      (1Darray) = signal on which the mask is applied;
%   OUTPUT
%   nanVel        (1Darray) = output signal with NaNs inserted;


notNanIdxs = find(~isnan(points(:,3)));
nanVel = points(:,3);
nanVel(notNanIdxs) = velocity;

%{
figure()
plot3(virtualPoints(:,1),virtualPoints(:,2), abs(nanVel), '.')

%}
end

