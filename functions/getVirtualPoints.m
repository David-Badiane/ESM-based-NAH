function [virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints)
%GETVIRTUALPOINTS Summary of this function goes here
%   Detailed explanation goes here


    % Virtual points (equivalent source)
    virtualPoints = violinInfos{4};
    virtualPoints(isnan(virtualPoints)) = 0;

    % The lattice is the minimum distance btw z positions of hologram and
    % equivalent sources. 
    %Since the array have different sizes, a for cycle is necessary
    diff = zeros(length(virtualPoints(:,3)));
    for ii = 1:length(virtualPoints(:,3))
        diff(ii) = min(abs(hologramPoints(:, 3) - virtualPoints(ii, 3)));
    end
    lattice = min(diff(diff~=0));
    
    virtualPoints = virtualPoints - [0, 0, lattice];
end

