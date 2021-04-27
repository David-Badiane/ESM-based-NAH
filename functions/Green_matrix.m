function [G, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints, omega)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S MATRIX                              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                 %%%%%%%%%%%%
%%% hologramPoints = x,y,z coordinates of the measurement points (2DArray) %%%%%%%%%%%%
%%% virtualPoints  = x,y,z coordinates of virtual sources (2DArray)        %%%%%%%%%%%%
%%% omega = eigenfrequencies in radians per seconds (1DArray)              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% G = Green's matrices (cell array)  for each omega                      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deleteIndexesVirt = find(isnan(virtualPoints(:,3)));
virtualPoints(deleteIndexesVirt,:) = [];

c = 343; % [m/s] speed of sound

% preallocate cell array
G = cell(length(omega),1);

% preallocate temporary matrix
G_w = zeros(length(hologramPoints(:, 1)), length(virtualPoints(:, 1)));  %Green's function init


for kk = 1:length(omega)
    k = omega(kk)/c;   % wave number

    for ii = 1:length(hologramPoints(:, 1)) 
        for jj = 1:length(virtualPoints(:, 1))      
            dist = norm((hologramPoints(ii, :) - virtualPoints(jj, :)), 2);  % L2 norm of (xr,yr,zr)-(xs,ys,zs)
            G_w(ii , jj) = (1/4/pi) * exp(-1i*k*dist)/dist;    % Green Matrix cell
        end   
    end
    
G{kk} = G_w;  % fill cell array
end

end

