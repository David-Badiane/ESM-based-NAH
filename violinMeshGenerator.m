%% Import grid ( downsampling or the already computed one )
baseFolder = pwd;
violinMeshesFolder = [baseFolder,'\violinMeshes'];
addpath(genpath('functions'));

readSTL = false;

if readSTL
    fv = stlread('NAH_ESM_mesh_A.stl');
    points = sortrows(fv.Points);
    [pts] = downsampling_regular(outMatrix, 128, 128);
else
    pts = readmatrix('grid128x128Fin.csv');
end

%% SYMMETRIZE A SCANNED GRID (asymmetries may arise because of sampling)

uniqueY = unique(pts(:,2));
uniqueX = unique(pts(:,1));

for ii = 1:length(uniqueY)
    idxsViolin = find(~ isnan(pts(:,3)));
    
    idxs = find(pts(:,2) == uniqueY(ii));
    x = pts(idxs,1);
    y = pts(idxs,2);
    z = pts(idxs,3);
    idxsPlot = find(~isnan(z));
    idxs = idxs(idxsPlot);
    
    if ~isempty(idxsPlot)
    figure(1)
    plot(x(idxsPlot),y(idxsPlot), '.');
    pause(0.1)
    
    figure(2)
    plot(pts(idxsViolin,1),pts(idxsViolin,2), '.')
    hold on;
    plot(x(idxsPlot),y(idxsPlot), '.', 'markerSize', 5);
    hold off;
    pause(0.1)

    disp([x(idxsPlot(1)), x(idxsPlot(end))])
      
    prompt = 'IdxLeft? ';
    idxL = input(prompt);

  
    prompt = ['IdxRight? end = ', num2str(length(idxs))];
    idxR = input(prompt);


    if ~isempty(idxL)
        pts(idxs(idxL),3) = NaN;
    end
    
    if ~isempty(idxR)
        pts(idxs(idxR),3) = NaN;
    end    
    end
end

%% Downsample a mesh - (violinMeshes generator);
cd(violinMeshesFolder)
nrows = 85; ncols = 55; 
fileName = ['grid',int2str(nrows),'x',int2str(ncols)];
saveData = true;
[outMatrix] = downsampling_regular(pts, nrows, ncols, fileName, saveData);
cd(baseFolder);
