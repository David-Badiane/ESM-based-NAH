function [G] = Green_matrix(r , s , omega)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S FUNCTIONS MATRIX      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: r  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the hologram plane measurement points     %%%%%%%%%%%%
%%%        or the surface coordinate points             %%%%%%%%%%%%
%%%        s  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the equivalent surface points             %%%%%%%%%%%%
%%%        omega = eigenfrequency array                 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: G_w  = cell array containing Green's       %%%%%%%%%%%%%
%%%                function matrices                   %%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 343; % [m/s] speed of sound

% preallocate cell array
G = cell(length(omega),1);

% preallocate temporary matrix
G_w = zeros(length(r(:, 1)), length(s(:, 1)));  %Green's function init

%conversion from [mm] to [m]
r = r*0.001;
s = s*0.001;


for kk = 1:length(omega)
    k = omega(kk)/c;   % wave number

    for ii = 1:length(r(:, 1)) 
        for jj = 1:length(s(:, 1))      
            dist = norm((r(ii, :) - s(jj, :)), 2);            % L2 norm of (xr,yr,zr)-(xs,ys,zs)
            G_w(ii, jj) = (1/4/pi) * exp(-1i*k*dist)/dist;    % Green Matrix cell
        end

    G{kk} = G_w;  % fill cell array
    end

end

end

