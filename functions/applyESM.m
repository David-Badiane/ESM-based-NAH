function [lossFx, ESM_metrics_table , ESM_Results] = applyESM( controlParams, boundZ, pressure, hologramPoints, normalPoints, violinMesh , omega,...
                            xData, yData , v_GT_vector, virtualPtsFilename,...
                            gridTablesNames, plotData, experimentalData )   

% applyESM - this function applies the whole ESM method, saving metrics and results

% INPUTS

% controlParams  ([3x1] Array) = control Params of the virtual points grids
%                                to be driven for minimization 
%                                [virtualPts z distance  scaleX   scaleY]
% boundZ         (double)      = to bound Z to a minimum value for minimization
% pressure       (1DArray)     = pressure vector
% hologramPoints (2DArray)     = hologram points coordinates
% normalPoints   (2DArray)     = normal points coordinates
% violinMesh     (2DArray)     = violin points coordinates
% omega          (double)      = radians frequency
% xData          (2DArray)     = X query measured points to interpolate
% yData          (2DArray)     = Y query measured points to interpolate
% v_GT_vector    (1DArray)     = velocity groundtruth vector
% virtualPtsFilename (string)  = string of the filename of virtual points
% gridTablesNames    (cell)    = names for the table
% plotData           (boolean) = to plot figures
% experimentalData

% OUTPUTS

% lossFx            (double) = loss function over which we minimize value
% ESM_metrics_table (table)  = table containing the metrics
% ESM_Results       (struct) = struct containing the results 
%                              virtual points weights, estimated veloicies,
%                              velocities surfaces

rho = 1.2;
zVal = controlParams(1);
scaleX = controlParams(2);
scaleY = controlParams(3);

% choose virtual points grid   
virtualPoints = readmatrix(virtualPtsFilename) ;
virtualPoints = [scaleX*virtualPoints(:,1)  scaleY*virtualPoints(:,2)  virtualPoints(:,3)];

if ~experimentalData
   virtualPoints = [virtualPoints(:,2), virtualPoints(:,1), virtualPoints(:,3)]; 
end

% tuple of data with metrics
ESM_metrics = zeros(1,length(gridTablesNames)-1);

% set z coordinate of the virtual points
% if boundZ is different than zero it is for minimizing and cap z to a
% minimum value
virtualPoints( ~isnan(virtualPoints(:,3)), 3) = - max([abs(zVal),boundZ]);

% compute Green functions matrix ( hologram --> virtual points ) 
[G_p, deleteIndexesVirt] = Green_matrix(hologramPoints , virtualPoints , omega ); G_p = G_p{1}; 

% compute gradient Green functions matrix 
% ( virtual points --> reconstruction plane = violin )
[G_v] = normalGradient(virtualPoints, violinMesh , omega, normalPoints); G_v = G_v{1};

% REGULARIZATION - L curve        
[U,s,V] = csvd (G_p);
lambda_L = l_curve (U,s,pressure);   % Tikhonov reg params
k_L = l_curve (U,s,pressure,'tsvd'); % Truncated-SVD reg param

% Compute 
[velocityErrors, ESM_Results] = errorVelocity(v_GT_vector, violinMesh,xData, yData, pressure,...
    G_p, G_v, lambda_L, k_L, omega, rho, experimentalData);
% store results
ESM_metrics(1) =  zVal;                                % 'zVal'
ESM_metrics(2) =  lambda_L;                            % 'lambda_L' 
ESM_metrics(3) =  k_L;                                 % 'k_L' 
ESM_metrics(4) =  velocityErrors.nmseTIK;             % nmseTIK_L
ESM_metrics(5) =  velocityErrors.nccTIK;              % nccTIK_L
ESM_metrics(6) =  velocityErrors.normcTIK;            % normcTIK_L
ESM_metrics(7) =  velocityErrors.reTIK;               % reTIK_L
ESM_metrics(8) =  velocityErrors.nmseTSVD;            % nmseTSVD_L
ESM_metrics(9) =  velocityErrors.nccTSVD;             % nccTSVD_L
ESM_metrics(10) =  velocityErrors.normcTSVD;          % normcTSVD_L
ESM_metrics(11) =  velocityErrors.reTSVD;             % reTSVD_L

lossTSVD = ESM_metrics(8) + 10*log10( 1 - ESM_metrics(9));
lossTIK  = ESM_metrics(4) + 10*log10( 1 - ESM_metrics(5));
lossFx = mean([lossTIK],2);

ESM_metrics_table = array2table(ESM_metrics, 'variableNames', gridTablesNames(2:end))

if plotData
    v_TSVD = ESM_Results.vTSVD;
    v_TIK = ESM_Results.vTIK;
    
    titleStr = ['  f = ',num2str(omega/2/pi), ...
                ' Hz   \\  RD = ', num2str(abs(round(virtualPoints(1,3),4))) , '  m' ];
    figureNum = 660;
    if experimentalData 
        reguFiguresExperimental(violinMesh, v_TSVD, v_TIK,  v_GT_vector,virtualPoints, titleStr, pressure  , figureNum);
    else
        reguFiguresSynthetic(violinMesh, v_TSVD, v_TIK,  v_GT_vector,virtualPoints, titleStr  ,  figureNum);                    
    end
end

end 
