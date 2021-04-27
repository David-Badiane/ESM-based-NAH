function [virtualPoints, lattice] = getVirtualPoints(violinInfos,hologramPoints, plotData)
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
    diffY = min(abs(y - circshift(y,1)));
    lattice = min([diffX, diffY]);
    lattice = lattice * 0.001;
    
    virtualPoints = virtualPoints - [0, 0, lattice];
    
    X =  reshape(virtualPoints(:,1), [16, 64]).'; 
    Y =  reshape(virtualPoints(:,2), [16, 64]).'; 
    Z =  reshape(virtualPoints(:,3), [16, 64]).';
    
    if plotData
    figure()
    surf(violinInfos{1},violinInfos{2},violinInfos{3});
    title('Violin surface');
    zlim([-100,100]);
    hold on 
    surf(X,Y,Z);
    end
    
end

