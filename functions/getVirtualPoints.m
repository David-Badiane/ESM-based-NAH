function [virtualPoints, lattice,deleteIndexes] = getVirtualPoints(violinInfos, hologramPoints, filename, plotData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GETVIRTUALPOINTS this function calculates the virtual points grid      %%%%%%%%%%%%
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
    violinMesh = violinInfos{4};
    % virtualPoints(isnan(virtualPoints)) = 0; THIS IS AN ERROR!!!!

    % The lattice is the minimum distance btw z positions of hologram and
    % equivalent sources. 
    x = unique(hologramPoints(:,1));
    y = unique(hologramPoints(:,2));
    diffX = min( abs(x - circshift(x,1)));
    diffY = min( abs(y - circshift(y,1)));
    lattice = min([diffX, diffY]);
    minZ = min(violinMesh(:,3));
%     virtualPoints(isnan(virtualPoints(:,3))) = [];
    idxs = find(~isnan(violinMesh(:,3)));
    deleteIndexes = find(isnan(violinMesh(:,3)));
%     virtualPoints(idxs,3) = minZ - lattice;
    

    %     [virtualPoints] = formVirtualPoints(virtualPoints, params{1}, params{2},params{3}, params{4}, params{5}, params{6});   
    disp(filename);
    virtualPoints = table2array(readtable([filename, '.csv'])); %/1000; 
   % virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), virtualPoints(:,3)]
    if plotData
        figure(5)
        plot3(violinInfos{4}(:,1),violinInfos{4}(:,2), minZ*violinInfos{4}(:,3),'.');
         title('Virtual Points ');
%         zlim([-0.1, 0.1]);
%         xlabel('x [m]')
%         ylabel('y [m]')
%         zlabel('z [m]')
        hold on 
        plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3),'x', 'markerSize',10 );
        view(2);
        hold off;
    end
    
end

