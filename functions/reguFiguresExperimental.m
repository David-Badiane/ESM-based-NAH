function [] = reguFiguresExperimental(violinMesh, v_TSVD, v_TIK, v_ex_vector, virtualPoints, titleStr, figureNum)
%REGUFIGURES Summary of this function goes here
%   Detailed explanation goes here
    
    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pY, pX]).'; 
    Y =  reshape(violinMesh(:,2), [pY, pX]).';

    v_TSVDF = addNans(violinMesh, v_TSVD); 
    v_TIKF = addNans(violinMesh, v_TIK);
%     v_ex = addNans(violinMesh, v_ex_vector);
    [XX,YY,surfV] = getVelocityGroundtruth(v_ex_vector);
%     surf_v = reshape(v_ex, [pY,pX]).';
    
    surf_vTSVD = reshape(v_TSVDF, [pY,pX]).';
    surf_vTIK = reshape(v_TIKF, [pY,pX]).';
    
    
    figure(figureNum+1)
    subplot 131
    surf(XX,YY, abs(surfV)); view(2);
    title('velocity groundTruth');
    subplot 132
    surf(X,Y, abs(surf_vTSVD )); view(2);
    title('TSVD velocity')
    subplot 133
    surf(X,Y, abs(surf_vTIK )); view(2);
    title('Tik velocity')
    sgtitle(titleStr);  
  
    figure(figureNum+2)
    plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3), 'o');
    hold on 
    plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 ); 
    xlabel('x'); ylabel('y');
    hold off;
    
    

end

