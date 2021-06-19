function [] = reguFiguresSynthetic(violinMesh, v_TSVD, v_TIK, v_ex_vector, virtualPoints, titleStr, figureNum)
%REGUFIGURES Summary of this function goes here
%   Detailed explanation goes here
    
    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pY, pX]).'; 
    Y =  reshape(violinMesh(:,2), [pY, pX]).';

    v_TSVD = addNans(violinMesh, v_TSVD); 
    v_TIK = addNans(violinMesh, v_TIK);
    v_ex = addNans(violinMesh, v_ex_vector);

    surfVelRecTSVD = reshape( v_TSVD, [pY, pX]).'; 
    surfVelRecTIK = reshape( v_TIK , [pY, pX]).'; 
    surfVel = reshape( v_ex , [pY, pX]).';

    
    figure(figureNum) 
        
    subplot 311
    surf(X, Y, abs(surfVel)); view(2);
    title('Exact velocity')
    subplot 312
    surf(X, Y, abs(surfVelRecTSVD));view(2);
    title('TSVD velocity')
    subplot 313
    surf(X, Y, abs(surfVelRecTIK));view(2);
    title('Tik velocity')
    sgtitle(titleStr);  
  
    figure(figureNum+1)
    plot3(violinMesh(:,1), violinMesh(:,2), violinMesh(:,3), 'o');
    hold on 
    plot3(virtualPoints(:,1), virtualPoints(:,2), virtualPoints(:,3), '.', 'markerSize', 8 ); 
    xlabel('x'); ylabel('y');
    hold off;  
    pause(0.01);
end

