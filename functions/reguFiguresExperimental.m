function [] = reguFiguresExperimental(violinMesh, v_TSVD, v_TIK, v_ex_vector, virtualPoints, titleStr, measuredPressure,figureNum)
%REGUFIGURES Summary of this function goes here
%   Detailed explanation goes here

%   INPUTS
%   violinMesh        (2Darray) = ;
%   v_TSVD            (2Darray) = ;
%   v_TIK             (2Darray) = ;
%   v_ex_vector       (2Darray) = ;
%   virtualPoints     (2Darray) = ;
%   titleStr          (2Darray) = ;
%   measuredPressure  (2Darray) = ;
%   figureNum         (2Darray) = ;

    velocityFilename = 'velocity_Data';
    pressureFilename = 'pressure_Data';
    pressureData = readmatrix([pressureFilename,'.csv']);


    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pY, pX]).'; 
    Y =  reshape(violinMesh(:,2), [pY, pX]).';

    v_TSVDF = addNans(violinMesh, v_TSVD); 
    v_TIKF = addNans(violinMesh, v_TIK);
%     v_ex = addNans(violinMesh, v_ex_vector);

[XX,YY,surfV] = getVelocityGroundtruth(v_ex_vector, velocityFilename, figureNum-1);
%     surf_v = reshape(v_ex, [pY,pX]).';
    
    surf_vTSVD = reshape(v_TSVDF, [pY,pX]).';
    surf_vTIK = reshape(v_TIKF, [pY,pX]).';
    
    
    figure(figureNum+1)
    subplot 141
    surf(XX,YY, abs(surfV)/max(abs(surfV), [], 'all') ); view(2);
    title('velocity groundTruth'); %colorbar;
    
    subplot 142
    surf(X,Y, abs(surf_vTSVD)/max(abs(surf_vTSVD), [], 'all') ); view(2);
    title('TSVD velocity'); %colorbar;
%     caxis([0 max(abs(surfV(:)))])

    subplot 143
    surf(X,Y, abs(surf_vTIK )/max(abs(surf_vTIK), [], 'all')); view(2);
    title('Tik velocity'); %colorbar;
    sgtitle(titleStr);  
%     caxis([0 max(abs(surfV(:)))])

    subplot 144
    xPress = reshape(pressureData(:,1), [8,8]).';
    yPress = reshape(pressureData(:,2), [8,8]).';
    zPress = reshape(measuredPressure, [8,8]).';
    
    surf(xPress, yPress, abs(zPress)); view(2);

    figure(figureNum+2)
    plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3), 'o');
    hold on 
    plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 ); 
    xlabel('x'); ylabel('y');
    hold off;
end

