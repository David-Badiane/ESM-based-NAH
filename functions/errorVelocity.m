function [velocityErrors] = errorVelocity(v_ex_vector, violinMesh, ...
xData, yData, measuredPressure, G_p, G_v, lambda_L, k_L,... 
omega, rho)

% ERRORVELOCITY this function calculate the error through different metrics
% [nmseTIK, nccTIK, normcTIK, reTIK] and save it into a struct

%   INPUT
%   v_ex_vector      (array)  = vector of the velocities groung truth ;
%   violinMesh       (2Darray) = mesh of the geometry to interpolate with v_ex_vector in 
%                               order to find the correspondent value of velocity in the 
%                               same indeces;
%   xData            (array)   = value of x for the interpolation;
%   yData            (array)   = value of y for the interpolation;
%   measuredPressure (array)   = hologram pressure for the given radians frequency 
%                               omega;
%   G_p              (2Darray) = Green's matrix for the pressure;
%   G_v              (2Darray) = Green's matrix for the velocity;
%   lambda_L         (double)  = Thikonov parameter;
%   k_L              (double)  = TSVD parameter;
%   omega            (double)  = radians frequency where evaluate the function;
%   rho              (double)  = air density;

%   OUPUT
%   velocityErrors   (struct) = struct where the metric values are stored;

names = {'nmseTIK' 'nccTIK' 'normcTIK' 'reTIK' 'nmseTSVD' 'nccTSVD' 'normcTSVD' 'reTSVD' };

%ERRORVELOCITY 
pX = length(unique(violinMesh(:,1)));
pY = length(unique(violinMesh(:,2)));

% compact singular value decomposition
[U,s,V] = csvd (G_p);

% exact velocity norm
normV = norm(v_ex_vector,2);
[Qs, v_TSVD, v_TIK] = reguResults( k_L, lambda_L, measuredPressure, omega, rho, G_p, G_v );
q_TIK = Qs.qTIK; 
q_TSVD = Qs.qTSVD;

% adds nan to create a mesh to interpolate
v_TIK_Fin = addNans(violinMesh, v_TIK);
% interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
v_TIK_Fin = interpGrid([violinMesh(:,1) violinMesh(:,2) v_TIK_Fin], xData.', yData.', pX, pY, false);
% make vector and delete NaNs
v_TIK_Fin = v_TIK_Fin(:);
cancelindex = find(isnan(v_TIK_Fin));
v_TIK_Fin(cancelindex) = [];

% metrics
[nmseTIK, nccTIK, normcTIK, reTIK] = metrics(v_TIK_Fin, v_ex_vector, normV);

  
% adds nan to create a mesh to interpolate
v_TSVD_Fin = addNans(violinMesh, v_TSVD);
% interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
v_TSVD_Fin = interpGrid([violinMesh(:,1) violinMesh(:,2) abs(v_TSVD_Fin)], xData.', yData.', pX, pY, false);
v_TSVD_Fin = v_TSVD_Fin(:);
% cancel nan from velocity vector
cancelindex = find(isnan(v_TSVD_Fin));
v_TSVD_Fin(cancelindex) = [];
    
[nmseTSVD, nccTSVD, normcTSVD, reTSVD] = metrics(v_TSVD_Fin, v_ex_vector, normV);

velocityErrors = struct(names{1}, nmseTIK,  names{2}, nccTIK, names{3}, normcTIK, names{4}, reTIK,...
    names{5}, nmseTSVD, names{6}, nccTSVD, names{7}, normcTSVD, names{8}, reTSVD);
end


function [nmse, ncc, normc, re] = metrics(vRegu, vGroundtruth, normV)
    % NMSE
    vGroundTruth = addNans(violinMesh, vGroundtruth); 
%      figure(1111)
%     plot(1:length(vRegu), abs(vRegu), 'lineWidth', 1.2);
%     hold on; 
%     plot(1:length(vGroundtruth), abs(vGroundtruth))
%     hold off;
% 
   
    nmse  = 10*log10(norm( abs(vGroundtruth)/max(abs(vGroundtruth)) ...
                           - abs(vRegu)/max(abs(vRegu)) )^2 / ...
                        (norm(abs(vGroundtruth)/max(abs(vGroundtruth)))^2) );
    
%     nmse = 10*log10(norm( vGroundtruth - vRegu )^2 / ...
%                       (normV)^2) ;
    %NCC
    ncc   = (abs(vRegu)'*abs(vGroundtruth)) / (norm(vRegu)*normV);
    
    %Normalized Correlation
    normc = abs(vRegu'*vGroundtruth) / (norm(vRegu)*normV);
   
    %Reconstruction Error (Relative Error)
    re    = (norm((vGroundtruth - vRegu) ,2)^2) / (normV)^2 ;
end

