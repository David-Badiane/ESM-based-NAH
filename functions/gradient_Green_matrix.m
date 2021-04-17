function [G] = gradient_Green_matrix(s , r , G_in, omega, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE GREEN'S FUNCTIONS             %%%%%%%
%%% DERIVATIVE ALONG THE NORMAL SURFACE'S DIRECTION          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: s  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the surface coordinate points             %%%%%%%%%%%%
%%%        r  =   matrix containing the coordinates     %%%%%%%%%%%%
%%%        of the equivalent surface points             %%%%%%%%%%%%
%%%        omega = eigenfrequency array                 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT: G_v  = cell array containing Green's       %%%%%%%%%%%%%
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
            der_x = (G_in(
        end   
    end
    
G{kk} = G_w;  % fill cell array
end


end

