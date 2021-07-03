function [] = reguFiguresExperimental(violinMesh, v_TSVD, v_TIK, v_GT_vector, virtualPoints, titleStr, measuredPressure,figureNum)
%reguFiguresExperimental - plots estimation figures for experimental data

%   INPUTS

%   violinMesh      (2DArray) = matrix with the violin Points (nPoints x 3) 
%   v_TSVD          (1DArray) = estimated velocity TSVD
%   v_TIK           (1DArray) = estimated velocity TIK
%   v_GT_vector     (1DArray) = groundtruth velocity vector
%   virtualPoints   (2DArray) = matrix with virtual Points (nPoints x 3)
%   titleStr        (string)  = string of the title of the images
%   measuredPressure(1Darray) = hologram pressure vector for the given omega;
%   figureNum       (double)  = figure number

% OUTPUTS
% ~

    velocityFilename = 'velocity_Data';
    pressureFilename = 'pressure_Data';
    pressureData = readmatrix([pressureFilename,'.csv']);


    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pY, pX]).'; 
    Y =  reshape(violinMesh(:,2), [pY, pX]).';

    v_TSVDF = addNans(violinMesh, v_TSVD); 
    v_TIKF = addNans(violinMesh, v_TIK);

[XX,YY,surfV] = getVelocityGroundtruth(v_GT_vector, velocityFilename, figureNum-1);
    
    surf_vTSVD = reshape(v_TSVDF, [pY,pX]).';
    surf_vTIK = reshape(v_TIKF, [pY,pX]).';  

    figure(figureNum+1)
    
    subplot 311
    surf(YY,-XX, abs(flip(surfV,2)),'EdgeColor','none'); view(2);
    title('Exact velocity'); colorbar;
    xlabel('x [m]')
    ylabel('y [m]')
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;
    
    subplot 312
    surf(Y,-X, abs(flip(surf_vTSVD,2)),'EdgeColor','none'); view(2);
    title('TSVD velocity'); colorbar;
    xlabel('x [m]')
    ylabel('y [m]')
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;
    
    subplot 313
    surf(Y,-X, abs(flip(surf_vTIK,2) ),'EdgeColor','none'); view(2);
    title('Tik velocity'); colorbar;
    sgtitle(titleStr,'fontSize', 18 );  

    xlabel('x [m]')
    ylabel('y [m]')
    
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;
     
     
    a=1; 
%     subplot 144
%     xPress = reshape(pressureData(:,1), [8,8]).';
%     yPress = reshape(pressureData(:,2), [8,8]).';
%     zPress = reshape(measuredPressure, [8,8]).';
%     
%     surf(xPress, yPress, abs(zPress)); view(2);

%     figure(figureNum+2)
%     plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3), 'o');
%     hold on 
%     plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 ); 
%     xlabel('x'); ylabel('y');
%     hold off;
end

