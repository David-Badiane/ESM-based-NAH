function genVirtualPoints(pts, fileName, controller, zVal,virtualPtsFolder, saveData)

% INPUTS
% pts  [nPts x 3]  (2DArray) = points of the mesh
% fileName          (string) = name of the .csv file containing the points   
%
%   ex ['VP_1'] {.csv NOT NEEDED}
%
% controller        (double) = 0, gen rectangular grids
% controller = 1, gen circular grids
% controller = 2, gen ellipsoidal grids
% controller = 3, gen circular + border
% controller = 4, gen ellipsoidal + border
% controller = 5, gen inner + border
% CONTROLLER = 6, GEN BORDER ONLY
% controller = 7, gen inner only
% controller = 8, gen rect + border + inner
% zVal             (double) = value of the z coordinate of the grid
% virtualPtsFolder (string) = filepath to the folder containing the virtual points
% saveData        (boolean) = if true creates the csv file of name fileName

% OUTPUT
% ~

baseFolder = pwd;
cd(virtualPtsFolder)
switch controller
    case 0 % gen rectangular grids
        xRect = input('input rectangle x edge: ');
        yRect = input('input rectangle y edge: ');
        xCenter = input('input  x center: ');
        yCenter = input('input  y center: ');

        [virtualPts] = rectVirtPoints(pts, xRect, yRect, xCenter, yCenter, zVal);
        
    case 1
        
        radius = input('input circle radius: ');
        innerRadius = input('input inner radius: ');        
        [virtualPts] = ellipseVirtualPoints(pts, radius, radius,innerRadius, innerRadius, zVal);
         
    case 2
        radius_x = input('input ellpise radius x: ');
        radius_y = input('input ellpise radius y: ');
        innerRadius_x = input('input inner radius x: ');
        innerRadius_y = input('input inner radius y: ');

        [virtualPts] = ellipseVirtualPoints(pts, radius_x, radius_y,innerRadius_x, innerRadius_y, zVal);
        
    case 3
        disp("circle");
        radius = input('input radius: ');
        innerRadius = input('input inner radius: ');     
        [virtualPts] = ellipseVirtualPoints(pts, radius, radius,innerRadius, innerRadius, zVal);
         
        disp("border");
        xBorder = input('how much of the border along x? (0 to 0.5): ');
        yBorder = input('how much of the border along y? (0 to 0.5): ');

        [border] = borderVirtualPoints(pts,  xBorder, yBorder, zVal);
        idxs = find(~isnan(border(:,3)));
        virtualPts(idxs,3) = border(idxs,3);
   
    case 4
        disp('ellipse')
        radius_x = input('input ellipse radius x: ');
        radius_y = input('input ellipse radius y: ');
        innerRadius_x = input('input inner radius x: ');
        innerRadius_y = input('input inner radius y: ');
  
        [virtualPts] = ellipseVirtualPoints(pts, radius_x, radius_y,innerRadius_x, innerRadius_y, zVal);
        
        disp("border");
        xBorder = input('how much of the border along x? (0 to 0.5): ');
        yBorder = input('how much of the border along y? (0 to 0.5): ');

        [border] = borderVirtualPoints(pts,  xBorder, yBorder, zVal);
        idxs = find(~isnan(border(:,3)));
        virtualPts(idxs,3) = border(idxs,3);

    case 5
        
   xBorder = input('how much interior to the border along x? (0 to 0.5): ');
   yBorder = input('how much interior to the border along y? (0 to 0.5): ');
   
   [virtualPts] = innerVirtualPoints(pts,  xBorder, yBorder, zVal);

    disp("border");
    xBorder = input('how much of the border along x? (0 to 0.5): ');
    yBorder = input('how much of the border along y? (0 to 0.5): ');

    [border] = borderVirtualPoints(pts,  xBorder, yBorder, zVal);
    idxs = find(~isnan(border(:,3)));
    virtualPts(idxs,3) = border(idxs,3); 
    
    case 6
        disp("border");
        xBorder = input('how much of the border along x? (0 to 0.5): ');
        yBorder = input('how much of the border along y? (0 to 0.5): ');

        [virtualPts] = borderVirtualPoints(pts,  xBorder, yBorder, zVal);
        
   case 7
   disp("inner"); 
   xBorder = input('how much interior to the border along x? (0 to 0.5): ');
   yBorder = input('how much interior to the border along y? (0 to 0.5): ');
    
   [virtualPts] = innerVirtualPoints(pts,  xBorder, yBorder, zVal);
   
    case 8
    disp(' ' ); disp('rectangle'); disp(' ' );
    xRect = input('input rectangle x edge: ');
    yRect = input('input rectangle y edge: ');
    xCenter = input('input  x center: ');
    yCenter = input('input  y center: ');

    [virtualPts] = rectVirtPoints(pts, xRect, yRect, xCenter, yCenter, zVal); 
    meshPts = table2array(readtable('grid128x128Fin.csv'));
    virtualPts(~isnan(meshPts(:,3)), 3 ) = nan;
    sparser(virtualPts);
    scaler(virtualPts);
  
    disp("border");
    xBorder = input('how much of the border along x? (0 to 0.5): ');
    yBorder = input('how much of the border along y? (0 to 0.5): ');
    [border] = borderVirtualPoints(pts,  xBorder, yBorder, zVal);
    scaler(border); sparser(border);
    
    add2VPgrid = input('add to VP grid? (0, 1): ');
    
    if add2VPgrid ~= 0
        idxs = find(~isnan(border(:,3)));
        virtualPts(idxs,3) = border(idxs,3); 
    end
    
    disp("inner"); 
    xBorder = input('how much interior to the border along x? (0 to 0.5): ');
    yBorder = input('how much interior to the border along y? (0 to 0.5): ');
    
    [inner] = innerVirtualPoints(pts,  xBorder, yBorder, zVal);
    scaler(inner); sparser(inner);
    
    add2VPgrid = input('add to VP grid? (0, 1): ');
    if add2VPgrid ~= 0
        idxs = find(~isnan(inner(:,3)));
        virtualPts(idxs,3) = inner(idxs,3); 
    end   
end


        figure(1)
        plot3(virtualPts(:,1), virtualPts(:,2), virtualPts(:,3), '.');
        
        sparser(virtualPts);
        scaler(virtualPts);
        
        writeMat2File(virtualPts, [fileName,'.csv'], {'x' 'y' 'z'}, 3, true);
        
        downsampling = input('resample ? ( 0 or 1)');        
        if downsampling == 1
            nrows = input('nrows: ');
            ncols = input('ncols: ');

            [virtualPts] = downsampling_regular(virtualPts, nrows, ncols, fileName, saveData);
        end       
        cd(baseFolder);
end


function [] = sparser(virtualPts)        
    check = input('make sparser ? (0 or 1): ');
        while check ~= 0
            sparseMode = input(' how to ? (0 = modulo, 1 = random, 2 = modulo rand) :');
            if sparseMode == 1
                xCut =1; yCut = 1;
            else
                xCut = input('modulo x: ');
                yCut = input('modulo y: ');
            end         
            [virtualPts] = sparserVirtualPoints(virtualPts,sparseMode, xCut, yCut);
            check = input('still sparser? (0 or 1): ');
            figure(1)
            plot3(virtualPts(:,1), virtualPts(:,2), virtualPts(:,3), '.'); 
            view(2);
        end
end

function [] = scaler(virtualPts)
    check = input('scale ? (0 or 1): ');
    while check ~= 0
        xScale = input('x factor: ');
        yScale = input('y factor: ');
        virtualPts(:,1) = xScale * virtualPts(:,1); 
        virtualPts(:,2) = yScale * virtualPts(:,2);
        check = input('still scale? (0 or 1): ');
        figure(1)
        plot3(virtualPts(:,1), virtualPts(:,2), virtualPts(:,3), '.'); 
        view(2);        
    end  
end