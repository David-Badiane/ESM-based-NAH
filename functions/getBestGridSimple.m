function [lossFx] = getBestGridSimple( pressure, hologramPoints, normalPoints, violinMesh , omega,...
                           zVal, xData, yData , v_ex_vector,...
                           rho, pX, pY,gridTablesNames, transposeGrids, plotData, experimentalData )     

% choose virtual points grid   
virtualPtsFilename = 'VP_1.csv';
virtualPoints = table2array(readtable(virtualPtsFilename)) ;
if transposeGrids
   virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), virtualPoints(:,3)]; 
end

% tuple of data with metrics
ZreguDataLine = zeros(1,length(gridTablesNames)-1); 
virtualPoints( ~isnan(virtualPoints(:,3)), 3) = zVal;
% compute Green functions matrix ( hologram 2 virtual points ) 
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , omega );
G_p_omega = G_p{1}; 

% compute gradient Green functions matrix 
% ( virtual points 2 reconstruction plane = violin )
[G_v] = normalGradient(virtualPoints, violinMesh , omega, normalPoints);
G_v_omega = G_v{1};

% APPROACH 1 ) L curve solution       
% 1) Inverse - individuation of regularization parameter (lambda) 
[U,s,V] = csvd (G_p_omega);
lambda_l = l_curve (U,s,pressure);
k_l = l_curve (U,s,pressure,'tsvd');

[velocityErrorsL, desiredAlphaL] = errorVelocity(v_ex_vector, violinMesh, xData, yData, pressure,...
    G_p_omega, G_v_omega, [lambda_l lambda_l],...
    [k_l k_l], 1, 1, omega, rho, deleteIndexesVirt, pX, pY, experimentalData);


% store results
ZreguDataLine(1) =  zVal;                           % 'zVal'
ZreguDataLine(2) =  lambda_l;                            % 'lambda_L' 
ZreguDataLine(3) =  k_l;                                 % 'k_L' 
ZreguDataLine(4) =  velocityErrorsL.nmseTSVD;            % nmseTSVD_L
ZreguDataLine(5) =  velocityErrorsL.nccTSVD;             % nccTSVD_L
ZreguDataLine(6) =  velocityErrorsL.nmseTIK;             % nmseTIK_L
ZreguDataLine(7) =  velocityErrorsL.nccTIK;              % nccTIK_L


lossTSVD = ZreguDataLine(4) + 10*log10( 1 - ZreguDataLine(5));
lossTIK  = ZreguDataLine(6) + 10*log10( 1 - ZreguDataLine(6));
lossFx = mean([lossTIK],2);

reg = array2table(ZreguDataLine, 'variableNames', gridTablesNames)
if plotData
    [LQs, Lv_TSVD, Lv_TIK] = reguResults(ZreguDataLine(3), ZreguDataLine(2), pressure, omega, rho, G_p_omega, G_v_omega, virtualPoints );
    titleStr = ['L curve  f = ',num2str(omega/2/pi), ' Hz   \\  z = ', num2str(ZreguDataLine(1))];
    figureNum = 660;
    if experimentalData 
        reguFiguresExperimental(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , figureNum);
    else
        reguFiguresSynthetic(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , figureNum);                    
    end
end

end 
