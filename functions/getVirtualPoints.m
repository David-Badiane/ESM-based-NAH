function [virtualPoints, lattice,deleteIndexes] = getVirtualPoints(violinInfos,hologramPoints, plotData)
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
    minZ = min(virtualPoints(:,3));
%     virtualPoints(isnan(virtualPoints(:,3))) = [];
    idxs = find(~isnan(virtualPoints(:,3)));
    deleteIndexes = find(isnan(virtualPoints(:,3)));
    virtualPoints(idxs,3) = minZ - lattice;
    
    x = virtualPoints(:,1);
    y = virtualPoints(:,2);
    z = virtualPoints(:,3);

    xUnique = unique(x);
    yUnique = unique(y);

%     for ii = 1:length(xUnique)
%         idxs = find(x==xUnique(ii));
%         yCheck = y(idxs);
%         for jj = 1:length(yCheck)
%             if ii >2 || ii < length(xUnique) -2
%                 if mod(jj,2) ==1
%                     z(idxs(jj)) = nan;
%                 end
%             end
%         end
%     end
    
    for ii = 1:length(yUnique)
        idxs = find(y==yUnique(ii));
        xCheck = x(idxs);
        for jj = 1:length(xCheck)
            if ii >4 && ii < length(yUnique)-3
                if jj > 7 && jj < length(xCheck)-5
                    if mod(jj,3) == 1
                    else
                        z(idxs(jj)) = nan;
                    end
                end
            end
        end
    end
    
    scale = 1.2;
    virtualPoints = [scale *x(:), scale*y(:), z(:)];
    %virtualPoints = virtualPoints - [0, 0, lattice];
    %virtualPoints = sort(downsampling(virtualPoints, 440, 3));

%     X =  reshape(virtualPoints(:,1), [8, 8]).'; 
%     Y =  reshape(virtualPoints(:,2), [8, 8]).'; 
%     Z =  reshape(virtualPoints(:,3), [8, 8]).';
    
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

