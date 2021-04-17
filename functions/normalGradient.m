function [G_v] = normalGradient(hologramPoints , virtualPoints , omega, normalPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FUNCTION COMPUTES THE NORMAL GRADIENT OF THE GREEN'S MATRIX       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS                                                                 %%%%%%%%%%%%
%%% hologramPoints = x,y,z coordinates of the measurement points (2DArray) %%%%%%%%%%%%
%%% virtualPoints  = x,y,z coordinates of virtual sources (2DArray)        %%%%%%%%%%%%
%%% omega = eigenfrequencies in radians per seconds (1DArray)              %%%%%%%%%%%%
%%% normalPoints = x,y,z components of the normal vector (2DArray)         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS                                                                %%%%%%%%%%%%
%%% G_v = normal gradient of Green's matrix for each omega (cell array)    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343; % [m/s] speed of sound

% preallocate cell array
Gx = cell(length(omega),1);
Gy = cell(length(omega),1);
Gz = cell(length(omega),1);
normalGradientComponents = {Gx,Gy,Gz};
G_v = cell(length(omega),1);
% preallocate temporary matrix
G_w = zeros(length(hologramPoints(:, 1)), length(virtualPoints(:, 1)));  %Green's function init

%conversion from [mm] to [m]
hologramPoints = hologramPoints*0.001;
virtualPoints = virtualPoints*0.001;


for component = 1:length(normalGradientComponents)   % for each component (x,y,z)
    for kk = 1:length(omega)            % for each eigenfrequency
        k = omega(kk)/c;   % wave number

        for ii = 1:length(hologramPoints(:, 1))      % for each hologram point
            for jj = 1:length(virtualPoints(:, 1))   % for each virtual  point  
                % calculation of the distance vector
                distX = hologramPoints(ii, 1) - virtualPoints(jj, 1);        
                distY = hologramPoints(ii, 2) - virtualPoints(jj, 2); 
                distZ = hologramPoints(ii, 3) - virtualPoints(jj, 3); 
                distVector = [distX, distY, distZ]; 
                dist = sqrt(distVector(1)^2 + distVector(2)^2 + distVector(3)^2); % as L2 norm
                
                alpha = -(1i*k*dist + 1)/(dist^2);   % constant term from the derivative calculation
                G = (1/4/pi) * exp(-1i*k*dist)/dist;    % actual Green's function   
                
                gradient = distVector(component)*alpha* G;   % gradient is analytically like this
                
                if isnan(normalPoints(ii,component))   % multiply by zero if is a NaN
                    normalPoints(ii,component) = 0;
                end
                
                G_w(ii,jj) = gradient *normalPoints(ii,component);  % scalar product with normal vector component
            end   
        end

    normalGradientComponents{component}{kk} = G_w;  % fill cell array
    end
end


for kk = 1:length(omega)
    for ii = 1:length(hologramPoints(:, 1))      % for each hologram point
        for jj = 1:length(virtualPoints(:, 1))   % for each virtual  point  
            % retrieve G_v cell with L2 norm
            G_w(ii,jj) = sqrt(normalGradientComponents{1}{kk}(ii,jj)^2 +... 
                              normalGradientComponents{2}{kk}(ii,jj)^2 +...
                              normalGradientComponents{3}{kk}(ii,jj)^2 );
        end   
    end
    G_v{kk} = G_w;
end

end

