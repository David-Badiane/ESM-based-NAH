function [] = reguFiguresSynthetic(violinMesh, v_TSVD, v_TIK, v_GT_vector, virtualPoints, titleStr, figureNum)
% reguFiguresSynthetic - plots estimation figures for synthetic data

%   INPUTS

%   violinMesh      (2DArray) = matrix with the violin Points (nPoints x 3) 
%   v_TSVD          (1DArray) = estimated velocity TSVD
%   v_TIK           (1DArray) = estimated velocity TIK
%   v_GT_vector     (1DArray) = groundtruth velocity vector
%   virtualPoints   (2DArray) = matrix with virtual Points (nPoints x 3)
%   titleStr        (string)  = string of the title of the images
%   figureNum       (double)  = figure number

% OUTPUTS
% ~

    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pY, pX]).'; 
    Y =  reshape(violinMesh(:,2), [pY, pX]).';

    v_TSVD_F = addNans(violinMesh, v_TSVD); 
    v_TIK_F = addNans(violinMesh, v_TIK);
    v_ex = addNans(violinMesh, v_GT_vector);

    surfVelRecTSVD = reshape( v_TSVD_F, [pY, pX]).'; 
    surfVelRecTIK = reshape( v_TIK_F , [pY, pX]).'; 
    surfVel = reshape( v_ex , [ pY, pX]).';

    figure(figureNum) 
        
    subplot 311
    surf(X, Y, abs(surfVel),'EdgeColor','none'); view(2);
    xlabel('x [m]')
    ylabel('y [m]')
    title('Exact velocity')
    colorbar 
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;
    
    subplot 312
    surf(X, Y, abs(surfVelRecTSVD),'EdgeColor','none');view(2);
    xlabel('x [m]')
    ylabel('y [m]')
    title('TSVD velocity')
    colorbar 
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;
    
    subplot 313
    surf(X, Y, abs(surfVelRecTIK),'EdgeColor','none');view(2);
    xlabel('x [m]')
    ylabel('y [m]')
    title('Tik velocity')
    colorbar 
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    ax.FontSize = 12;    sgtitle(titleStr,'fontSize', 18 );  

    
    
    figure(figureNum+1)
    
    plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 );
    hold on 
    plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3), 'o');
    
    xlabel('x'); ylabel('y');
    hold off;  
    pause(0.01);
    
    figure(figureNum+2)
    plot3(violinMesh(:,1), violinMesh(:,2), abs(v_TIK_F));
    
end

