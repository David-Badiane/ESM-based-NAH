function [G_w] = Green_matrix(r , s , omega)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S FUNCTIONS MATRIX      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: r is a matrix containing the coordinates     %%%%%%%%%%%%
%%%          of the hologram plane measurement points   %%%%%%%%%%%%
%%%        s is a matrix containing the coordinates     %%%%%%%%%%%%
%%%          of the surface or equivalent surface points  %%%%%%%%%%
%%%        omega is the eigenfrequency of the plate     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: G_w is the Green's function matrix          %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_w = zeros(length(r(:, 1)), length(s(:, 1)));  %Green's function init

c = 343; % [m/s] speed of sound

for i = 1 : length(r_m(:, 1)) 
    for j = 1 : length(a_n(:, 1))      
        G_w(i, j, k) = (1/4/pi) * exp(-1i * omega/c *... 
            0.001*norm((r_m(i, :) - a_n(j, :)), 2))... %conversion from [mm] to [m]
            / 0.001*norm((r_m(i, :) - a_n(j, :)), 2);
    end
end


end

