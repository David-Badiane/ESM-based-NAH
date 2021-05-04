function [velocityErrors, desiredAlpha] = plotErrorVelocity(v_ex_vector, measuredPressureN, G_p_omega, G_v_omega, rangeTIK, rangeTSVD, numParamsTIK, numParamsTSVD, omega, rho, deleteIndexesVirt)

alphaTIK  = linspace(rangeTIK(1), rangeTIK(2), numParamsTIK);
alphaTSVD = round(linspace(rangeTSVD(1), rangeTSVD(2), numParamsTSVD));

nmseTIK  = zeros(1, numParamsTIK);
nccTIK  = zeros(1, numParamsTIK);
nmseTSVD = zeros(1, numParamsTSVD);
nccTSVD = zeros(1, numParamsTSVD);
[U,s,V] = csvd (G_p_omega);

for j = 1:numParamsTIK                 %norms are all norm-2
    q_TIK = (1/(1i*omega*rho)).*tikhonov(U,s,V, measuredPressureN  , alphaTIK(j)); % reconstructed source streght
    v_TIK = G_v_omega*q_TIK; 
    
    normV = norm(v_ex_vector ,2);    
    nmseTIK(j)  = 10*log(norm(v_TIK - v_ex_vector)^2 / (normV^2));
    nccTIK(j)   = abs((v_TIK'*v_ex_vector) / (norm(v_TIK,2)*normV));    
end


for j = 1:numParamsTSVD
    q_TSVD = (1/(1i*omega*rho)).*tsvd (U,s,V, measuredPressureN  , alphaTSVD(j)); 
    v_TSVD = G_v_omega*q_TSVD;   
    v_TSVD = abs(v_TSVD);
    
    normV = norm(v_ex_vector ,2);
    nmseTSVD(j) = 10*log(norm(v_TSVD - v_ex_vector)^2 / (normV^2));
    nccTSVD(j)  = abs(v_TSVD'*v_ex_vector) / (norm(v_TSVD,2)* normV);
end

errors = {nmseTIK, nccTIK, nmseTSVD, nccTSVD};
alphaVectors = {alphaTIK; alphaTSVD};

desiredAlpha = zeros(4,2);
names = {'nmseTIK' 'nccTIK' 'nmseTSVD' 'nccTSVD'};
namesAlpha = {'alphaTIK' 'alphaTSVD' };

figure()
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
    
    subplot (2,2, ii)
    plot(alphaVectors{alphaIndex}, errors{ii});
    hold on
    stem(alphaVectors{alphaIndex}(loc), val);
    xlabel(namesAlpha{alphaIndex});
    ylabel(names{ii});
    
    desiredAlpha(ii,2) = alphaVectors{alphaIndex}(loc);
end   

velocityErrors = struct(names{1}, nmseTIK,  names{2}, nccTIK, names{3}, nmseTSVD, names{4}, nccTSVD);


end

