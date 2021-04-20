function [virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION CALCULATES THE VIRTUAL POINTS GRID                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                 %%%%%%%%%%%%
%%% violinInfos = violin meshes and points (2DArray)                       %%%%%%%%%%%%
%%% hologramPoints = hologram points (2DArray)                             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% virtualPoints = virtual sources grid (2DArray)                         %%%%%%%%%%%%
%%% lattice = distance from the object surface (double)                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Detailed explanation goes here


    % Virtual points (equivalent source)
    virtualPoints = violinInfos{4};
    % virtualPoints(isnan(virtualPoints)) = 0; THIS IS AN ERROR!!!!

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

