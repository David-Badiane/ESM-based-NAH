function [velocityErrors, desiredAlpha] = errorVelocity(v_ex_vector, meshPoints, xData, yData, measuredPressureN, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexesVirt)
%ERRORVELOCITY 

v_ex_vector = zData(:,:,1); %!!!!!!!! ONLY FOR TEST

%  take the parameters
alphaTIK  = linspace(rangeTIK(1), rangeTIK(2), numParamsTIK);
alphaTSVD = round(linspace(rangeTSVD(1), rangeTSVD(2), numParamsTSVD));

% preallocate variables
nmseTIK  = zeros(1, numParamsTIK);
nccTIK  = zeros(1, numParamsTIK);
nmseTSVD = zeros(1, numParamsTSVD);
nccTSVD = zeros(1, numParamsTSVD); 

% compact singular value decomposition
[U,s,V] = csvd (G_p_omega);

% exact velocity norm
normV = norm(v_ex_vector,2);

% TIK cycle for every parameters
for ii = 1 %:numParamsTIK
    q_TIK = (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressureN  , alphaTIK(ii)); % reconstructed source streghts
    v_TIK = G_v_omega*q_TIK; 
    
    % adds nan to create a mesh to interpolate, CHECK IF NAN ARE ADDED
    v_TIK_Fin = addNans(meshPoints, v_TIK);
    v_TIK_Fin = reshape(v_TIK_Fin, [pX,pY]).';
    
    % interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
    v_TIK_Fin = interGrid(V_TIK_Fin, xData, yData, false);
    
    % reshape the vector to create the right mesh (?needed)
    % V_TIK_Fin = reshape(V_TIK_Fin, [length(xData(:,1)), length(yData(1,:))]);
    
    % cancel nan from velocity vector
    cancelindex = find(isnan(V_TIK_Fin));
    v_TIK_Fin(cancelindex) = [];
    v_ex_vector(cancelindex) = [];
    
    %NMSE
    nmseTIK(ii)  = 10*log(norm(v_TIK_Fin - v_ex_vector)^2 / (normV^2));
    
    %NCC, !!!!!!! CHECK COLUMN OR ROWS VECTOR, reshape to do the abs(?nedded)
    nccTIK(ii) = (abs(reshape(v_TIK_Fin, [numel(v_TIK_Fin),1]))'*...
        abs(rehsape(v_ex_vector, [numel(v_ex_vector),1])))/...
        (norm(abs(reshape(v_TIK_Fin, [numel(v_TIK_Fin),1]),2))*...
        norm(abs(rehsape(v_ex_vector, [numel(v_ex_vector),1])),2));


end

% TSVD cycle for every parameters
for jj = 1:numParamsTSVD
    q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressureN  , alphaTSVD(jj)); % reconstructed source streghts
    v_TSVD = G_v_omega*q_TSVD; 
    
    % adds nan to create a mesh to interpolate
    v_TSVD_Fin = addNans(meshPoints, v_TSVD);
    
    % interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
    v_TSVD_Fin = interGrid(V_TSVD_Fin, xData, yData, false);
    
    % reshape the vector to create the right mesh (?needed)
    % V_TSVD_Fin = reshape(V_TSVD_Fin, [length(xData(:,1)), length(yData(1,:))]);
    
     % cancel nan from velocity vector
    cancelindex = find(isnan(V_TVSD_Fin));
    v_TVSD_Fin(cancelindex) = [];
    v_ex_vector(cancelindex) = [];
    
    %NMSE
    nmseTSVD(jj) = 10*log(norm(v_TSVD_Fin - v_ex_vector)^2 / (normV^2));
    
    %NCC
    nccTSVD(jj)  = ((real(v_TSVD_Fin)'*real(v_ex_vector)) /...
        (norm(real(v_TSVD_Fin),2)*norm(real(v_ex_vector),2)) + ...
                (imag(v_TSVD)'*imag(v_ex_vector)) /...
                (norm(imag(v_TSVD),2)*norm(imag(v_ex_vector),2)))/2;
end

