uniqueY = unique(pts(:,2));
uniqueX = unique(pts(:,1));

for ii = 1:length(uniqueY)
    idxs = find(pts(:,2) == uniqueY(ii));
    x = pts(idxs,1);
    y = pts(idxs,2);
    z = pts(idxs,3);
    idxsPlot = find(~isnan(z));
    idxs = idxs(idxsPlot);
    
    if ~isempty(idxsPlot)
    figure(1)
    plot(x(idxsPlot),z(idxsPlot), '.');
    pause(0.1)
    
    figure(2)
    plot3(pts(:,1),pts(:,2),pts(:,3), '.')
    hold on;
    plot3(x,y,z, '.', 'markerSize', 5);
    hold off;
    pause(0.1)


%     prompt = 'IdxLeft? ';
%     idxL = input(prompt);
%     prompt = 'vals ? ';
%     valL = input(prompt);
% 
%     
%     prompt = ['IdxRight? end = ', num2str(length(idxs))];
%     idxR = input(prompt);
%     prompt = 'vals ? ';
%     valR = input(prompt);
% 
%     if ~isempty(idxL)
%         pts(idxs(idxL),3) = valL;
%     end
%     
%     if ~isempty(idxR)
%         pts(idxs(idxR),3) = valR;
%     end    
     end
end