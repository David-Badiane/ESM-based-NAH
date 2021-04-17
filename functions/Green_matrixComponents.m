function [G] = Green_matrixComponents(hologramPoints , virtualPoints , omega)
                                                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S FUNCTIONS MATRIX      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: s  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the hologram plane measurement points     %%%%%%%%%%%%
%%%        r  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the equivalent surface points             %%%%%%%%%%%%
%%%        omega = eigenfrequency array                 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: G_w  = cell array containing Green's       %%%%%%%%%%%%%
%%%                function matrices                   %%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 343; % [m/s] speed of sound

% preallocate cell array
Gx = cell(length(omega),1);
Gy = cell(length(omega),1);
Gz = cell(length(omega),1);
G = {Gx,Gy,Gz};
% preallocate temporary matrix
G_w = zeros(length(hologramPoints(:, 1)), length(virtualPoints(:, 1)));  %Green's function init

%conversion from [mm] to [m]
hologramPoints = hologramPoints*0.001;
virtualPoints = virtualPoints*0.001;


for component = 1:length(G)
    for kk = 1:length(omega)
        k = omega(kk)/c;   % wave number

        for ii = 1:length(hologramPoints(:, 1)) 
            for jj = 1:length(virtualPoints(:, 1))      
                dist = hologramPoints(ii, component) - virtualPoints(jj, component);            % L2 norm of (xr,yr,zr)-(xs,ys,zs)
                G_w(ii, jj) = (1/4/pi) * exp(-1i*k*dist)/dist;    % Green Matrix cell
            end   
        end

    G{component}{kk} = G_w;  % fill cell array
    end
end
end

