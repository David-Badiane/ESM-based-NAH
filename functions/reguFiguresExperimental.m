function [] = reguFiguresExperimental(violinMesh, v_TSVD, v_TIK, v_ex_vector, virtualPoints, titleStr, figureNum)
%REGUFIGURES Summary of this function goes here
%   Detailed explanation goes here
    
    pX = length(unique(violinMesh(:, 1)));
    pY = length(unique(violinMesh(:, 2))); 
   
    X =  reshape(violinMesh(:,1), [pX, pY]); 
    Y =  reshape(violinMesh(:,2), [pX, pY]);

    v_TSVD = addNans(violinMesh, v_TSVD); 
    v_TIK = addNans(violinMesh, v_TIK);
    %v_ex = addNans(violinMesh, v_ex_vector);

    
    figure(figureNum) 
        
  
    subplot 121
    plot3(violinMesh(:,1), violinMesh(:,2), abs(v_TSVD));
    title('TSVD velocity')
    subplot 122
    plot3(violinMesh(:,1), violinMesh(:,2), abs(v_TIK));
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

