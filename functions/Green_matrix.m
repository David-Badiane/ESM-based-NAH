function [G] = Green_matrix(r , s , omega)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S FUNCTIONS MATRIX      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: r  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the hologram plane measurement points     %%%%%%%%%%%%
%%%        s  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the surface or equivalent surface points  %%%%%%%%%%%%
%%%        omega = eigenfrequency of the plate          %%%%%%%%%%%%
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

<<<<<<< HEAD
for kk = 1:length(omega)
    k = omega(kk)/c;   % wave number

    for ii = 1:length(r(:, 1)) 
        for jj = 1:length(s(:, 1))      
            dist = norm((r(ii, :) - s(jj, :)), 2);            % L2 norm of (xr,yr,zr)-(xs,ys,zs)
            G_w(ii, jj) = (1/4/pi) * exp(-1i*k*dist)/dist;    % Green Matrix cell
        end
=======
for i = 1 : length(r(:, 1)) 
    for j = 1 : length(s(:, 1))      
        G_w(i, j, k) = (1/4/pi) * exp(-1i * omega/c *... 
            0.001*norm((r(i, :) - s(j, :)), 2))... %conversion from [mm] to [m]
            / 0.001*norm((r(i, :) - s(j, :)), 2);
>>>>>>> 382eeb2c66a0de3861df7ef8836dc1d15a4544bf
    end

    G{kk} = G_w;  % fill cell array
end

end

