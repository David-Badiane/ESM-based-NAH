function [velocityErrors, ESM_Results] = errorVelocity(v_GT_vector, violinMesh, ...
xData, yData, measuredPressure, G_p, G_v, lambda_L, k_L,... 
omega, rho)

% ERRORVELOCITY this function perform the whole ESM method,
% interpolates over the measured points( xDAta, yData) the estimated velocity,
% calculates the metrics [nmseTIK, nccTIK, normcTIK, reTIK]
% and saves it back into a struct

%   INPUT
%   v_GT_vector      (1DArray)   = vector of the velocity groundtruth;
%   violinMesh       (2DArray)   = mesh of the geometry;
%   xData            (2DArray)   = x matrix for the interpolation over the measured points;
%   yData            (2DArray)   = y matrix for the interpolation over the measured points;
%   measuredPressure (1DArray)   = hologram pressure for the given radians frequency 
%                               omega;
%   G_p              (2DArray) = Green's matrix for the pressure;
%   G_v              (2DArray) = Green's matrix for the velocity;
%   lambda_L         (double)  = Thikonov regularization parameter;
%   k_L              (double)  = TSVD regularization parameter;
%   omega            (double)  = radians frequency where evaluate the function;
%   rho              (double)  = medium density;

%   OUTPUT
%   velocityErrors   (struct) = struct where the metric values are stored;
%   ESM_Results      (struct) = struct containing ew sources weights and
%                               estimated velocities;

names_Metrics = {'nmseTIK' 'nccTIK' 'normcTIK' 'reTIK' 'nmseTSVD' 'nccTSVD' 'normcTSVD' 'reTSVD' };
names_Results = {'qTIK' 'qTSVD' 'vTIK' 'vTSVD' 'vTIK_Fin' 'vTSVD_Fin'};


% compact singular value decomposition
[U,s,V] = csvd (G_p);

% exact velocity norm

[Qs, v_TSVD, v_TIK] = reguResults( k_L, lambda_L, measuredPressure, omega, rho, G_p, G_v );
q_TIK = Qs.qTIK; 
q_TSVD = Qs.qTSVD;

v_TIK_Fin = getSameOrdering(v_TIK, violinMesh, xData, yData);   
% metrics
[nmseTIK, nccTIK, normcTIK, reTIK] = metrics(v_TIK_Fin, v_GT_vector);

  

v_TSVD_Fin = getSameOrdering(v_TSVD, violinMesh, xData, yData);   
[nmseTSVD, nccTSVD, normcTSVD, reTSVD] = metrics(v_TSVD_Fin, v_GT_vector);

velocityErrors = struct(names_Metrics{1}, nmseTIK,  names_Metrics{2}, nccTIK, names_Metrics{3}, normcTIK, names_Metrics{4}, reTIK,...
    names_Metrics{5}, nmseTSVD, names_Metrics{6}, nccTSVD, names_Metrics{7}, normcTSVD, names_Metrics{8}, reTSVD);

ESM_Results = struct(names_Results{1}, q_TIK , names_Results{2}, q_TSVD , names_Results{3}, v_TIK ,...
    names_Results{4}, v_TSVD , names_Results{5},v_TIK_Fin ,names_Results{6},v_TSVD_Fin );
end


function [nmse, ncc, normc, re] = metrics(vRegu, vGroundtruth)
    % NMSE
    %     figure(111)
%     plot(1:length(vRegu), abs(vRegu), 'lineWidth', 1.2);
%     hold on; 
%     plot(1:length(vGroundtruth), abs(vGroundtruth))
%     hold off;

%     IN CASE OF NORMALIZATION
%     vGroundtruth = abs(vGroundtruth)/max(abs(vGroundtruth));
%     vRegu = abs(vRegu)/max(abs(vRegu));
     normV = norm(vGroundtruth,2);
    
%     nmse  = 10*log10(norm( abs(vGroundtruth)/max(abs(vGroundtruth)) ...
%                            - abs(vRegu)/max(abs(vRegu)) )^2 / ...
%                         (norm(abs(vGroundtruth)/max(abs(vGroundtruth)))^2) );    
      nmse = 10*log10(norm( vGroundtruth - vRegu )^2 / ...
                       (normV)^2) ;
    %NCC
    ncc   = (abs(vRegu)'*abs(vGroundtruth)) / (norm(vRegu)*normV);
    
    %Normalized Correlation
    normc = abs(vRegu'*vGroundtruth) / (norm(vRegu)*normV);
   
    %Reconstruction Error (Relative Error)
    re    = (norm((vGroundtruth - vRegu) ,2)^2) / (normV)^2 ;
end


function vRegu_Fin = getSameOrdering(vRegu, violinMesh, xData, yData)

    pX = length(unique(violinMesh(:,1)));
    pY = length(unique(violinMesh(:,2)));
    
    % adds nan to create a mesh to interpolate
    vRegu_Fin = addNans(violinMesh, vRegu);
    % interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
    vRegu_Fin = interpGrid([violinMesh(:,1) violinMesh(:,2) abs(vRegu_Fin)], xData.', yData.', pX, pY, false);
    vRegu_Fin = vRegu_Fin(:);
    % cancel nan from velocity vector
    vRegu_Fin(isnan(vRegu_Fin)) = [];
end





