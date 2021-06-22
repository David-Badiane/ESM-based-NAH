function [lossFx] = getBestGridSimple( pressure, hologramPoints, normalPoints, violinMesh , omega,...
                           controlParams, xData, yData , v_ex_vector,...
                           rho,gridTablesNames, transposeGrids, plotData, experimentalData )     

              zVal = controlParams(1);
              scale = controlParams(2);
% choose virtual points grid   
virtualPtsFilename = 'VP_1.csv';
virtualPoints = scale*readmatrix(virtualPtsFilename) ;
if transposeGrids
   virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), virtualPoints(:,3)]; 
end

% tuple of data with metrics
ZreguDataLine = zeros(1,length(gridTablesNames)-1); 
virtualPoints( ~isnan(virtualPoints(:,3)), 3) = zVal;
% compute Green functions matrix ( hologram 2 virtual points ) 
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , omega );
G_p = G_p{1}; 

% compute gradient Green functions matrix 
% ( virtual points 2 reconstruction plane = violin )
[G_v] = normalGradient(virtualPoints, violinMesh , omega, normalPoints);
G_v = G_v{1};

% APPROACH 1 ) L curve solution       
% 1) Inverse - individuation of regularization parameter (lambda) 
[U,s,V] = csvd (G_p);
lambda_L = l_curve (U,s,pressure);
k_L = l_curve (U,s,pressure,'tsvd');

[velocityErrors] = errorVelocity(v_ex_vector, violinMesh, ...
                xData, yData, pressure, G_p, G_v, lambda_L, k_L, omega, rho);
% store results
ZreguDataLine(1) =  zVal;                                % 'zVal'
ZreguDataLine(2) =  lambda_L;                            % 'lambda_L' 
ZreguDataLine(3) =  k_L;                                 % 'k_L' 
ZreguDataLine(4) =  velocityErrors.nmseTIK;             % nmseTIK_L
ZreguDataLine(5) =  velocityErrors.nccTIK;              % nccTIK_L
ZreguDataLine(6) =  velocityErrors.normcTIK;            % normcTIK_L
ZreguDataLine(7) =  velocityErrors.reTIK;               % reTIK_L
ZreguDataLine(8) =  velocityErrors.nmseTSVD;            % nmseTSVD_L
ZreguDataLine(9) =  velocityErrors.nccTSVD;             % nccTSVD_L
ZreguDataLine(10) =  velocityErrors.normcTSVD;          % normcTSVD_L
ZreguDataLine(11) =  velocityErrors.reTSVD;             % reTSVD_L

lossTSVD = ZreguDataLine(4) + 10*log10( 1 - ZreguDataLine(5));
lossTIK  = ZreguDataLine(8) + 10*log10( 1 - ZreguDataLine(9));
lossFx = mean([lossTIK],2);

reg = array2table(ZreguDataLine, 'variableNames', gridTablesNames(2:end))

if plotData
    [LQs, Lv_TSVD, Lv_TIK] = reguResults(ZreguDataLine(3), ZreguDataLine(2), pressure, omega, rho, G_p, G_v, virtualPoints );
    titleStr = ['L curve  f = ',num2str(omega/2/pi), ' Hz   \\  z = ', num2str(ZreguDataLine(1)), ' scale = ' num2str(scale)];
    figureNum = 660;
    if experimentalData 
        reguFiguresExperimental(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr, pressure  , figureNum);
    else
        reguFiguresSynthetic(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  ,  figureNum);                    
    end
end

end 
