function [virtualPoints, lattice,deleteIndexes] = getVirtualPoints(violinInfos, hologramPoints, params, plotData)
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
    x = unique(hologramPoints(:,1));
    y = unique(hologramPoints(:,2));
    diffX = min( abs(x - circshift(x,1)));
    diffY = min( abs(y - circshift(y,1)));
    lattice = min([diffX, diffY]);
    minZ = min(virtualPoints(:,3));
%     virtualPoints(isnan(virtualPoints(:,3))) = [];
    idxs = find(~isnan(virtualPoints(:,3)));
    deleteIndexes = find(isnan(virtualPoints(:,3)));
    virtualPoints(idxs,3) = minZ - lattice;
    
% [virtualPoints] = formVirtualPoints(virtualPoints, scale, offSet, xCut,yCut ,deleteX, deleteY,xBorder, yBorder, active)

    [virtualPoints] = formVirtualPoints(virtualPoints, params{1}, params{2},params{3}, params{4}, params{5}, params{6});   
    
    if plotData
        figure()
        surf(violinInfos{1},violinInfos{2},violinInfos{3});
        title('Violin surface');
        zlim([-0.1, 0.1]);
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('z [m]')
        hold on 
        %surf(X,Y,Z);
        plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3),'.', 'markerSize',10 );
    end
    
end

