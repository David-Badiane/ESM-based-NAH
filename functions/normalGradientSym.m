function [G_v] = normalGradientSym(virtualPoints , platePoints , omega, normalPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE NORMAL GRADIENT OF THE GREEN'S MATRIX       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                 %%%%%%%%%%%%
%%% virtualPoints = x,y,z coordinates of the virtual point (2DArray)       %%%%%%%%%%%%
%%% platePoints  = x,y,z coordinates of surface points (2DArray)           %%%%%%%%%%%%
%%% omega = eigenfrequencies in radians per seconds (1DArray)              %%%%%%%%%%%%
%%% normalPoints = x,y,z components of the normal vector (2DArray)         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% G_v = normal gradient of Green's matrix for each omega (cell array)    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343; % [m/s] speed of sound

% preallocate cell array

deleteIndexesVirt = find(isnan(virtualPoints(:,3)));
virtualPoints(deleteIndexesVirt,:) = [];

deleteIndexesPlate = find(isnan(platePoints(:,3)));
platePoints(deleteIndexesPlate,:) = [];

normalPoints(deleteIndexesPlate, :) = [];

G_v = cell(length(omega),1);
% preallocate temporary matrix
G_w = zeros(length(virtualPoints(:, 1)), length(platePoints(:, 1)));  %Green's function init

syms x y z
k = omega/c;
dist = sqrt( x^2 + y^2 + z^2 );
green = 1/4/pi* (exp(-1i*k*dist))/(dist);

greenX = diff(green, x,1);
greenY = diff(green, y,1);
greenZ = diff(green, z,1);

for kk = 1:length(omega)
    for ii = 1:length(virtualPoints(:, 1))      % for each virtual point
        for jj = 1:length(platePoints(:, 1))   % for each surface point 
            % retrieve G_v cell with L2 norm
            distX =  virtualPoints(ii, 1) - platePoints(jj, 1);        
            distY =  virtualPoints(ii, 2) - platePoints(jj, 2); 
            distZ =  virtualPoints(ii, 3) - platePoints(jj, 3); 
            
            vars = [distX, distY, distZ];
            Gx = eval(subs(greenX,[x,y,z],vars))* normalPoints(jj,1);
            Gy = eval(subs(greenY,[x,y,z],vars))* normalPoints(jj,2);
            Gz = eval(subs(greenZ,[x,y,z],vars))* normalPoints(jj,3);          
            
            G_w(ii,jj) = sqrt(Gx^2 + Gy^2 + Gz^2 );
        end   
    end
    G_v{kk} = G_w;
end

end

