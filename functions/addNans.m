function [nanVelocity] = addNans(virtualPoints, velocity)
%ADDNANS Summary of this function goes here
%   Detailed explanation goes here
nanIndexes = find(isnan(virtualPoints(:,3)));
nanVelocity = zeros(size(virtualPoints(:,1)));
nanVelocity(nanIndexes) = NaN;
notNanIndexes = find(~isnan(nanVelocity));
nanVelocity(notNanIndexes) = velocity;
end

