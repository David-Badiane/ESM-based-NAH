function [pressureErrors, desiredAlpha] = plotErrorPressure(measuredPressure, G_p_omega , measuredPressureN , omega , rho , rangeTIK , rangeTSVD, numParamsTIK, numParamsTSVD)


alphaTIK  = linspace(rangeTIK(1), rangeTIK(2), numParamsTIK);
alphaTSVD = round(linspace(rangeTSVD(1), rangeTSVD(2), numParamsTSVD));

nmseTIK  = zeros(1, numParamsTIK);
nccTIK  = zeros(1, numParamsTIK);
nmseTSVD = zeros(1, numParamsTSVD);
nccTSVD = zeros(1, numParamsTSVD);


for j = 1:numParamsTIK                 %norms are all norm-2
    q_TIK = (1/(1i*omega*rho)).*Tikhonov_SVD(G_p_omega , measuredPressureN  , alphaTIK(j)); % reconstructed source streght
    p_TIK = 1i*omega*rho*G_p_omega*q_TIK; % reconstructed velocity with Tikhonov      
    normP = norm(measuredPressure,2);   
    nmseTIK(j)  = 10*log(norm(p_TIK - measuredPressure)^2 / (normP^2));
    nccTIK(j)   = (p_TIK'*measuredPressure) / (norm(p_TIK,2)*normP);    
end


for j = 1:numParamsTSVD
    q_TSVD = (1/(1i*omega*rho)).*TSVD(G_p_omega, measuredPressureN  , alphaTSVD(j)); 
    p_TSVD = 1i*omega*rho*G_p_omega*q_TSVD;
    normP = norm(measuredPressure,2);
    nmseTSVD(j) = 10*log(norm(p_TSVD - measuredPressure)^2 / (normP^2));
    nccTSVD(j)  = (p_TSVD'*measuredPressure) / (norm(p_TSVD,2)* normP);
end

errors = {nmseTIK, nccTIK, nmseTSVD, nccTSVD};
alphaVectors = {alphaTIK; alphaTSVD};

desiredAlpha = zeros(4,2);
names = {'nmseTIK' 'nccTIK' 'nmseTSVD' 'nccTSVD'};
namesAlpha = {'alphaTIK' 'alphaTSVD' };

for ii = 1:length(errors)
    discriminator = mod(ii,2);
    if discriminator == 1
      [val, loc] =  min( errors{ii});   % NMSE    
    else
      [val, loc] =  max( errors{ii});   %NCC
    end
    desiredAlpha(ii,:) = [val,loc];
    
    figure(661)
    
    if ii == 1 || ii == 2
        alphaIndex = 1;
    else 
        alphaIndex = 2;
    end
    
    subplot (2,2, ii)
    plot(alphaVectors{alphaIndex}, errors{ii});
    hold on
    stem(alphaVectors{alphaIndex}(loc), val);
    xlabel(namesAlpha{alphaIndex});
    ylabel(names{ii});
    desiredAlpha(ii,1) = alphaVectors{alphaIndex}(loc);
end   

pressureErrors = struct(names{1}, nmseTIK,  names{2}, nccTIK, names{3}, nmseTSVD, names{4}, nccTSVD);
desiredAlpha = desiredAlpha(:,1);

end

