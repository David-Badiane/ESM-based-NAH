function [nanVel] = addNans(points, velocity)
%ADDNANS this function creates an array of size [length(points(:,3),1]
% where the output array has the same nan indexes of points(:,3)
% and its not nan indexes are filled by the value of velocity

%   INPUTS
%   virtualPoints (array) = signal on which the mask is applied;
%   velocity      (array) = singal where NaNs index positons are taken;

%   OUTPUT
%   nanVel        (array) = output signal with NaNs inserted;


% nanIndexes = find(isnan(virtualPoints(:,3)));
% nanVelocity = zeros(size(virtualPoints(:,1)));
% nanVelocity(nanIndexes) = nan;
% notNanIndexes = find(~isnan(nanVelocity));
% nanVelocity(notNanIndexes) = velocity;
% figure()
% plot3(virtualPoints(:,1),virtualPoints(:,2),nanVelocity, '.')


notNanIdxs = find(~isnan(points(:,3)));
nanVel = points(:,3);
nanVel(notNanIdxs) = velocity;

%{
figure()
plot3(virtualPoints(:,1),virtualPoints(:,2), abs(nanVel), '.')

%}
end

