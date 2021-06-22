function [reguData, ZreguDatas] = getBestGrid(nEqSourceGrids, pressure, hologramPoints, normalPoints, violinMesh , omega,...
                           nZpoints, zCenter, zSearch, xData, yData , v_ex_vector,...
                           rho, gridTablesNames, transposeGrids, plotData, experimentalData )     
% GETBESTGRID find the best virtual points grid
                     
%   INPUTS
%   nEqSourceGrids   (double)  = number of grids of equivalent sources;
%   pressure         (array)   = vector of pressure measured;
%   hologramPoints   (2Darray) = mesh of the hologram points;
%   normalPoints     (2Darray) = matrix of components of normal vector of the reconstrucet surface;
%   violinMesh       (2Darray) = mesh of the geometry
%   omega            (double)  = radians frequency where evaluate the function;
%   nZpoints         (double)  = number of z point where to search
%   zCenter          (double)  = where the grid is centered;
%   zSearch          (double)  = percentage of the deviation;
%   xData            (array)   = value of x for the interpolation;
%   yData            (array)   = value of y for the interpolation;
%   v_ex_vector      (array)   = vector of the velocities groung truth ;
%   rho              (double)  = air density;
%   gridTablesNames  (cell)    = where the name of the grid and its values are stored;
%   transposeGrids   (boolean) = choose if the grid have to be transposed;
%   plotData         (boolean) = choose if the plot have to be displayed
%   experimentalData (boolean) = choose if the plot shown;

%   OUPUTS
%   reguData         (2Darray) = best grid regularization data;
%   ZreguDatas       (2Darray) = regularization data for each z searched;           
                       
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
            ZreguDataLine(1) =  zVals(zz);                           % 'zVal'
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
       
            ZreguData = [ZreguData; ZreguDataLine];
            
       end

    % evaluate ZreguData(:, 4:7) to choose best z
    % accordingly to loss function L = ( NMSE + 20LOG(1-NCC) )
    lossTSVD = ZreguData(:,4)/max(abs(ZreguData(:,4))) + 10*log10( 1 - ZreguData(:,5))/ max(abs(10*log10( 1 - ZreguData(:,5))));
    lossTIK  = ZreguData(:,6)/max(abs(ZreguData(:,6))) + 10*log10( 1 - ZreguData(:,7))/ max(abs(10*log10( 1 - ZreguData(:,7))));
    lossFx = mean([lossTSVD 2*lossTIK],2);
    [zminVal, zminLoc] = min(lossFx);
    ZtempTable = array2table( ZreguData , 'VariableNames', gridTablesNames(2:end))
    reguDataLine = [jj, ZreguData(zminLoc,:)];
    reguData = [reguData; reguDataLine];
    
    reguTable = array2table(reguDataLine, 'variableNames', gridTablesNames)
    
    if plotData
        figureNum = 660;
       [LQs, Lv_TSVD, Lv_TIK] = reguResults( reguDataLine(4), reguDataLine(3), pressure, omega, rho, G_p, G_v, virtualPoints );
        titleStr = ['L curve  f = ',num2str(omega/2/pi), '  \\  z = ', num2str(reguDataLine(2)), ' \\ ' , virtualPtsFilename];
        if experimentalData 
            reguFiguresExperimental(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , figureNum);
        else
            reguFiguresSynthetic(violinMesh, Lv_TSVD, Lv_TIK,  v_ex_vector,virtualPoints, titleStr  , pressure, figureNum);                    
        end
    end

    ZreguDatas{jj} = ZreguData;
end 
end