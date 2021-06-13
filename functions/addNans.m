function [nanVel] = addNans(virtualPoints, velocity)
%ADDNANS Summary of this function goes here
%   Detailed explanation goes here
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

