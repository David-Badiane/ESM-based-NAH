function genVirtualPoints(pts, fileName, controller, zVal,virtualPtsFolder)

% controller = 0, gen rectangular grids
% controller = 1, gen circular grids
% controller = 2, gen ellipsoidal grids
% controller = 3, gen circular + border
% controller = 4, gen ellipsoidal + border
% controller = 5, gen inner + border
% CONTROLLER = 6, GEN BORDER ONLY
% controller = 7, gen inner only
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
             
   xBorder = input('how much interior to the border along x? (0 to 0.5): ');
   yBorder = input('how much interior to the border along y? (0 to 0.5): ');
    
   [virtualPts] = innerVirtualPoints(pts,  xBorder, yBorder, zVal);
end


        figure(1)
        plot3(virtualPts(:,1), virtualPts(:,2), virtualPts(:,3), '.');
        
        check = input('make sparser ? (0 or 1): ');
        while check == 1
            controller = input(' how to ? (0 = modulo, 1 = random, 2 = modulo rand) :');
            if controller == 1
                xCut =1; yCut = 1;
            else
                xCut = input('modulo x: ');
                yCut = input('modulo y: ');
            end         
            [virtualPts] = sparserVirtualPoints(virtualPts,controller, xCut, yCut);
            check = input('still sparser? (0 or 1): ');
        end
        
        check = input('scale ? (0 or 1): ');
        while check == 1

            xScale = input('x factor: ');
            yScale = input('y factor: ');
            virtualPts(:,1) = xScale * virtualPts(:,1); 
            virtualPts(:,2) = yScale * virtualPts(:,2); 

            check = input('still scale? (0 or 1): ');
        end
        
        
        
        
        writeMat2File(virtualPts, [fileName,'.csv'], {'x' 'y' 'z'}, 3, true);
        cd(baseFolder);
end

        