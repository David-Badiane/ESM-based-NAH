function [nanVel] = addNans(virtualPoints, velocity)
%ADDNANS this function create a binary mask in an array taking index 
% positions where NaNs are present in another array

%   INPUTS
%   virtualPoints (array) = signal on which the mask is applied;
%   velocity      (array) = singla where NaNs index positons are taken;

%   OUTPUT
%   nanVel        (array) = output signal with NaNs inserted;


% nanIndexes = find(isnan(virtualPoints(:,3)));
% nanVelocity = zeros(size(virtualPoints(:,1)));
% nanVelocity(nanIndexes) = nan;
% notNanIndexes = find(~isnan(nanVelocity));
% nanVelocity(notNanIndexes) = velocity;
% figure()
% plot3(virtualPoints(:,1),virtualPoints(:,2),nanVelocity, '.')


notNanIdxs = find(~isnan(virtualPoints(:,3)));
nanVel = virtualPoints(:,3);
nanVel(notNanIdxs) = velocity;

%{
figure()
plot3(virtualPoints(:,1),virtualPoints(:,2), abs(nanVel), '.')

%}
end

