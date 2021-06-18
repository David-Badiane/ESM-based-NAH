function [reguData, ZreguDatas] = getBestGrid(nEqSourceGrids, pressure, hologramPoints, normalPoints, violinMesh , omega,...
                           nZpoints, zCenter, zSearch, xData, yData , v_ex_vector,...
                           rho, pX, pY,gridTablesNames, transposeGrids, plotData, experimentalData )     

reguData = [];
ZreguDatas = cell(nEqSourceGrids,1);
zVals = linspace(zCenter*(1- zSearch),zCenter *(1+ zSearch), nZpoints );

for jj = 1:nEqSourceGrids
        
       % choose virtual points grid   
       virtualPtsFilename = ['VP_', int2str(jj), '.csv'];
       
       virtualPoints = table2array(readtable(virtualPtsFilename)) ;
       % deleteIndexes = find(isnan(virtualPoints(:,3)));
       if transposeGrids
           virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), virtualPoints(:,3)]; 
       end
       ZreguData = [];
       
       % iterate along z
       for zz = 1:nZpoints 
            % tuple of data to choose best z
            ZreguDataLine = zeros(1,length(gridTablesNames)-1); 
            virtualPoints( ~isnan(virtualPoints(:,3)), 3) = zVals(zz);
            
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
                [k_l k_l], 1, 1, omega, rho, deleteIndexesVirt, pX, pY);
            
            
%             % APPROACH 2) metrics parametrization
            rangeTIK = [0,100]; % range of value for the regularization parameter
            TSVDup = min([length(G_p_omega(1,:)), length(G_p_omega(:,1))]);
            rangeTSVD = [1,TSVDup]; % range of value for the regularization parameter
            numParamsTIK = 1e2;
            numParamsTSVD = TSVDup-1;

            [velocityErrors, desiredAlpha] = errorVelocity(v_ex_vector, violinMesh,...
                xData, yData, pressure, G_p_omega, G_v_omega, rangeTIK, rangeTSVD,...
                numParamsTIK, numParamsTSVD, omega, rho, deleteIndexesVirt, pX, pY);

            % store results
            ZreguDataLine(1) =  zVals(zz);                           % 'zVal'
            ZreguDataLine(2) =  lambda_l;                            % 'lambda_L' 
            ZreguDataLine(3) =  k_l;                                 % 'k_L' 
            ZreguDataLine(4) =  desiredAlphaL(3,1);                  % nmseTSVD_L
            ZreguDataLine(5) =  desiredAlphaL(4,1);                  % nccTSVD_L
            ZreguDataLine(6) =  desiredAlphaL(1,1);                  % nmseTIK_L
            ZreguDataLine(7) =  desiredAlphaL(2,1);                  % nccTIK_L
% 
            ZreguDataLine(8) =  desiredAlpha(3,2);                   % lambda_nmse_M 
            ZreguDataLine(9) =  desiredAlpha(4,2);                   % k_nmse_M
            ZreguDataLine(10) = desiredAlpha(1,2);                   % k_ncc_M 
            ZreguDataLine(11) = desiredAlpha(2,2);                   % lambda_ncc_M
            ZreguDataLine(12) = desiredAlpha(3,1);                   % nmseTSVD_M
            ZreguDataLine(13) = desiredAlpha(4,1);                   % nccTSVD_M
            ZreguDataLine(14) = desiredAlpha(1,1);                   % nmseTIK_M
            ZreguDataLine(15) = desiredAlpha(2,1);                   % nccTIK_M 
            
            ZreguData = [ZreguData; ZreguDataLine];
            
       end

    % evaluate ZreguData(:, 4:7) to choose best z
    % accordingly to loss function L = ( NMSE + 20LOG(1-NCC) )
    lossTSVD = ZreguData(:,9)/max(abs(ZreguData(:,9))) + 10*log10( 1 - ZreguData(:,10))/ max(abs(10*log10( 1 - ZreguData(:,10))));
    lossTIK  = ZreguData(:,8)/max(abs(ZreguData(:,8))) + 10*log10( 1 - ZreguData(:,11))/ max(abs(10*log10( 1 - ZreguData(:,11))));
    lossFx = mean([lossTSVD 2*lossTIK],2);
    [zminVal, zminLoc] = min(lossFx);
    ZtempTable = array2table( ZreguData , 'VariableNames',gridTablesNames(2:end))
    reguDataLine = [jj, ZreguData(zminLoc,:)];
    reguData = [reguData; reguDataLine];
    
    reg = array2table(reguDataLine, 'variableNames', gridTablesNames)
    if plotData
       [LQs, Lv_TSVD, Lv_TIK] = reguResults( reguDataLine(4), reguDataLine(3), pressure, omega, rho, G_p_omega, G_v_omega );
            %[LQs, Lv_TSVD, Lv_TIK] = reguResults( ZreguDataLine(10), ZreguDataLine(11), pressure, omega, rho, G_p_omega, G_v_omega);
        titleStr = ['L curve   z = ', num2str(ZreguDataLine(2)), ' ' , virtualPtsFilename];
        figureNum = 660;
        if experimentalData 
            reguFiguresExperimental(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , figureNum);
        else
            reguFiguresSynthetic(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , figureNum);                    
        end
    end

    ZreguDatas{jj} = ZreguData;
end 
end