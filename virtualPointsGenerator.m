% folders
baseFolder = pwd;
virtualPointsFolder = [baseFolder,'\VP_Grids'];
addpath(genpath('functions'));
addpath(genpath('Data'));
addpath('violinMeshes');
addpath(virtualPointsFolder)

%% Virtual Points generator

% ATTENTION: measures in mm

pts = table2array(readtable('grid128x128Fin.csv'));

edgeX = max(pts(:,1)) - min(pts(:,1));
edgeY = max(pts(:,2)) - min(pts(:,2));

zVal = 0; % <-- lattice
nGrids = 10;

% from last file automatically
filesList = ls(virtualPointsFolder);
filesList(1:2,:) = [];
fileNums = [];

for ii = 1:length(filesList(:,1))
    idx = find(filesList(ii,:)== '.');
    fileNums = [fileNums; eval(filesList(ii,4:idx))];
end

fileNums = max(fileNums);
fileNums = sort(fileNums);

% all already set and debugged

for ii = 1:nGrids
    disp('');
    disp(' contr = 0 rectangular')
    disp('contr = 1, circular grids')
    disp('contr = 2, ellipsoidal grids')
    disp('contr = 3, circular + border')
    disp('contr = 4, ellipsoidal + border')
    disp('contr = 5, inner + border')
    disp('contr = 6, BORDER ONLY')
    disp('contr = 7, inner only')
    disp('contr = 8, rect + border + inner')
    disp('');
    disp(['edgeX: ', num2str(edgeX),' edgeY: ', num2str(edgeY)]); 
    
    controller = input('choose kind of grid(0-8) :');
    genVirtualPoints(pts,['VP_',int2str(ii)], controller, zVal,virtualPointsFolder);
end
