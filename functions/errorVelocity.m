function [velocityErrors, desiredAlpha] = errorVelocity(v_ex_vector, violinMesh, ...
xData, yData, measuredPressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK,...
 numParamsTSVD, omega, rho, deleteIndexesVirt, pX, pY, experimentalData)

%ERRORVELOCITY 

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
for ii = 1:numParamsTIK
    q_TIK = (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressure  , alphaTIK(ii)); % reconstructed source streghts
    v_TIK = G_v_omega*q_TIK; 
    v_TIK_Fin = v_TIK;
    if experimentalData
        % adds nan to create a mesh to interpolate
        v_TIK_Fin = addNans(violinMesh, v_TIK);
        % interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
        v_TIK_Fin = interpGrid([violinMesh(:,1) violinMesh(:,2) abs(v_TIK_Fin)], xData, yData, pX, pY, false);

        v_TIK_Fin = v_TIK_Fin(:);
    end
    % cancel nan from velocity vector
    cancelindex = find(isnan(v_TIK_Fin));
    v_TIK_Fin(cancelindex) = [];
    
    %NMSE
    nmseTIK(ii)  = 10*log(norm(v_TIK_Fin - abs(v_ex_vector))^2 / (normV^2));
    
    nccTIK(ii) = (abs(v_TIK_Fin)'*abs(v_ex_vector)) / (norm(abs(v_TIK_Fin),2)*norm(abs(v_ex_vector),2));


end

% TSVD cycle for every parameters
for jj = 1:numParamsTSVD
    
    q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressure  , alphaTSVD(jj)); % reconstructed source streghts
    v_TSVD = G_v_omega*q_TSVD; 
    
    % adds nan to create a mesh to interpolate
    v_TSVD_Fin = addNans(violinMesh, v_TSVD);
    
    % interpolate this mesh with v_ex_vector to find the points on the same indeces to do the subtraction
    v_TSVD_Fin = interpGrid([violinMesh(:,1) violinMesh(:,2) abs(v_TSVD_Fin)], xData, yData, pX, pY, false);
    
    v_TSVD_Fin = v_TSVD_Fin(:);
    
    % cancel nan from velocity vector
    cancelindex = find(isnan(v_TSVD_Fin));
    v_TSVD_Fin(cancelindex) = [];
    
    %NMSE
    nmseTSVD(jj)  = 10*log(norm(v_TSVD_Fin - abs(v_ex_vector))^2 / (normV^2));
    
    nccTSVD(jj) = (abs(v_TSVD_Fin)'*abs(v_ex_vector)) / (norm(abs(v_TSVD_Fin),2)*norm(abs(v_ex_vector),2));


end

errors = {nmseTIK, nccTIK, nmseTSVD, nccTSVD};
alphaVectors = {alphaTIK; alphaTSVD};

desiredAlpha = zeros(4,2);


names = {'nmseTIK' 'nccTIK' 'nmseTSVD' 'nccTSVD'};
namesAlpha = {'alphaTIK' 'alphaTSVD' };

% figure()
for ii = 1:length(errors)
    discriminator = mod(ii,2);
    if discriminator == 1
      [val, loc] =  min( errors{ii});   % NMSE    
    else
      [val, loc] =  max( errors{ii});   %NCC
    end
    desiredAlpha(ii,:) = [val,loc];
    
    if ii == 1 || ii == 2
        alphaIndex = 1;
    else 
        alphaIndex = 2;
    end

%     subplot (2,2, ii)
%     plot(alphaVectors{alphaIndex}, errors{ii});
%     hold on
%     stem(alphaVectors{alphaIndex}(loc), val);
%     xlabel(namesAlpha{alphaIndex});
%     ylabel(names{ii});   
    desiredAlpha(ii,2) = alphaVectors{alphaIndex}(loc);

end   

velocityErrors = struct(names{1}, nmseTIK,  names{2}, nccTIK, names{3}, nmseTSVD, names{4}, nccTSVD);
end

