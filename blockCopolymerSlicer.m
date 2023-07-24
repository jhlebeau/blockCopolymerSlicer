%Each function in this program contains a description of what it does, 
%please read through the program to get a better understanding of what each
%function does and the overall capabilites of the program.

% Version 1.0.0 (July 24, 2023)


function blockCopolymerSlicer
% This is equivalent to a main() or runner() function, the code in this
% fuction will be executed when the program is run.

    [MIN_BLOCK_VOL_FRAC, SLICES, SLICE_DEPTH, NUM_CELLS, INVERT_COLOR, ...
     MINORITY_BLOCK_COLOR_1, MINORITY_BLOCK_COLOR_2, ...
     MAJORITY_BLOCK_COLOR, NETWORK_TYPE] = setup(); 
    % Calling setup will save default values into these named variables.
    % They can be used, but do not have to be. 
    
    filePath = 'C:\'; % Modify this filePath to save figures to your desired location
    
    %Below are two examples to show some functionality of this program. 

    plotSlice(1,1,2,0.4,5,0.25,'double gyroid',false);
    %This plots a slice of the double gyroid network on the (112) plane.
    %The volume fraction of the minority network is 0.4, five unit cells
    %are plotted in the x and y directions, and the slice occurs 1/4 of the
    %way through a unit cell (depth = 0.25). The minority networks are 
    % white and the majority network is black (invertColor = false).

    plotSliceColorNetworks(1,1,1,0.42,5,0.55,'dd',MINORITY_BLOCK_COLOR_1,MINORITY_BLOCK_COLOR_2,MAJORITY_BLOCK_COLOR);
    %This plots a slice of the double diamond network on the (111) plane.
    %The volume fraction of the minority network is 0.42, five unit cells
    %are plotted in the x and y directions, and the slice occurs near the
    %center of a unit cell (depth = 0.55). The networks are colored using
    %the default colors located in the setup() function below.   

end


function [MIN_BLOCK_VOL_FRAC, SLICES, SLICE_DEPTH, NUM_CELLS, INVERT_COLOR, ...
     MINORITY_BLOCK_COLOR_1, MINORITY_BLOCK_COLOR_2, ...
     MAJORITY_BLOCK_COLOR, NETWORK_TYPE] = setup()
% This function can be useful to assign variables to values to be used in
% calling other functions in the program. It is not necessary to use this
% function, and values can be directly input as parameters when calling
% functions in the blockCopolymerSlicer function. Each variable is
% described below:

    MIN_BLOCK_VOL_FRAC = 0.4; 
    %Input minority block volume fraction.
    
    SLICES = 51; 
    % This is the number of slices made within a unit cell when using 
    % exportSlices.
    % Default: 51
    
    SLICE_DEPTH = 0.25; 
    % This is the position within a unit cell where the slice will be made for
    % plotSlice and similar functions that make a single slice. It is measured
    % from the bottom of the unit cell.
    % Give slice depth as a fraction of a unit cell (must be between 0 and 1).
    
    NUM_CELLS = 5; 
    % Number of unit cells plotted in x and y directions. Minimum = 1.
    
    MINORITY_BLOCK_COLOR_1 = [0.5,0.5,1]; 
    % Color for one minority block network, given as an RGB triplet in a 1x3 vector.
    % Used in any of the functions that color the networks separately.
    % NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255].
    % Default: [0.5,0.5,1] (light blue)
    
    MINORITY_BLOCK_COLOR_2 = [1,0,0]; 
    % Color for other minority block network, given as an RGB triplet in a 1x3 vector.
    % Used in any of the functions that color the networks separately.
    % NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255].
    % Default: [1,0,0] (red)
    
    MAJORITY_BLOCK_COLOR = [0,0,0];
    % Color for the majority block network, given as an RGB triplet in a 1x3 vector. 
    % Used in any of the functions that color the networks separately.
    % NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255].
    % Default: [0,0,0] (black)
    
    NETWORK_TYPE = 'double gyroid';
    % Which network is going to be sliced, either double diamond or 
    % double gyroid. Not case sensitive.
    % For double diamond, the program recognizes: "diamond", "double
    % diamond", and "dd"
    % For double gyroid, the program recognizes: "gyroid", "double gyroid",
    % and "dg"

    INVERT_COLOR = false;
    % Invert black and white in plotSlice, exportSlice, and plotAndSaveSlice.
    % This is a logical value (true/false) that inverts the color of the plot. 
    % When false, minority networks are white and the majority network is 
    % black. When true, minority networks are black and the majority
    % network is white.

end


function type = interpretNetworkType(networkType)
% THIS FUNCTION SHOULD NOT BE CALLED DIRECTLY
% This function interprets the network type input by the user. It could 
% be modified to accept more phrases.
% Inputs:
% networkType - A string that identifies the network type to be sliced
% Outputs:
% type - A string that identifies the network type to be sliced using only
% 'diamond' or 'gyroid' so that the program can interpret this later on

    inputType = lower(networkType);
    
    if strcmp(inputType,'diamond') || strcmp(inputType,'double diamond') || strcmp(inputType,'dd')

        type = 'diamond';

    elseif strcmp(inputType,'gyroid') || strcmp(inputType,'double gyroid') || strcmp(inputType,'dg')

        type = 'gyroid';

    else

        error("Network type could not be recongized. For double diamond, the program " + ...
            "recognizes: diamond, double diamond, and dd. For double gyroid, the program " + ...
            "recognizes: gyroid, double gyroid, and dg. It is not case sensitive.");

    end


end


function t = calculateT(minBlockVolFrac,type)
% THIS FUNCTION SHOULD NOT BE CALLED DIRECTLY
% This function converts the user input minority block volume fraction 
% into a constant to be used in the level set equations.
% Inputs:
% minBlockVolFrac - The volume fraction of the minority block
% type - A string that identifies the network type to be sliced, either
% 'diamond' or 'gyroid'.
% Outputs:
% t - A constant that is used in the level set equation to produce a plot.

    if minBlockVolFrac > 0.5 || minBlockVolFrac < 0

        error("The minority block volume fraction must be between 0 and 0.5, although values " + ...
            "should be in the range 0.33 to 0.42, which is where the double diamond and " + ...
            "double gyroid typically occur.")

    elseif minBlockVolFrac > 0.45 || minBlockVolFrac < 0.3

        warning("It is recommended that the minority block volume fraction is between " + ...
            "0.33 and 0.42, which is where the double diamond and double gyroid typically occur.")

    end


    if strcmp(type,'diamond')

        if minBlockVolFrac < 0.11

            warning("This simulation produces inaccurate results for the " + ...
                "double diamond network if the minority block volume fraction" + ...
            " is set below 0.11. It is recommended that the minority " + ...
            "block volume fraction is between 0.33 and 0.42, which is " + ...
            "where the double diamond and double gyroid typically occur.")

        end

        t = -1.1784*minBlockVolFrac + 1.1971;

    else 

        t = -1.4355*minBlockVolFrac + 1.4888;

    end

end


function [uHat,vHat,wHat] = rotate001toNewAxis(newAxis)
% THIS FUNCTION SHOULD NOT BE CALLED DIRECTLY
% This function rotates the coordinate system such that the 3D slicing 
% plane can be properly displayed in 2D without distortion.
% Inputs:
% newAxis - The normal vector to the plane that is being sliced, given as a
% vertical 3x1 vector.
% Outputs:
% uHat - This is the modified x axis, x', in the new coordinate system.
% vHat - This is the modified y axis, y', in the new coordinate system.
% wHat - This is the modified z axis, z', in the new coordinate system.

    if newAxis == [0;0;0]

        error("(000) is Not a Valid Direction, Please Use a Non-zero Value.");

    elseif newAxis(1) == 0 && newAxis(2) == 0 

        uHat = [1;0;0];
        vHat = [0;1;0];
        wHat = [0;0;1];

    else

        n = [0;0;1];
        normN = norm(n);
        nHat = n/normN;
    
        k = newAxis;
        normK = norm(k);
        kHat = k/normK;
    
        b = cross(kHat,nHat);
        bHat = b/norm(b);

        theta = acosd(dot(kHat,nHat));
    
        q0 = cosd(theta/2);
        q1 = sind(theta/2)*bHat(1);
        q2 = sind(theta/2)*bHat(2);
        q3 = sind(theta/2)*bHat(3);
    
        Q = zeros(3,3);
        Q(1,1) = -(q0*q0+q1*q1-q2*q2-q3*q3);
        Q(1,2) = -2*(q1*q2-q0*q3);
        Q(1,3) = -2*(q1*q3+q0*q2);
        Q(2,1) = -2*(q1*q2+q0*q3);
        Q(2,2) = -(q0*q0-q1*q1+q2*q2-q3*q3);
        Q(2,3) = -2*(q2*q3-q0*q1);
        Q(3,1) = 2*(q1*q3-q0*q2);
        Q(3,2) = 2*(q2*q3+q0*q1);
        Q(3,3) = q0*q0-q1*q1-q2*q2+q3*q3;
    
        uHat = Q*[1;0;0];
        vHat = Q*[0;1;0];
        wHat = Q*[0;0;1];
    
        % x' = (u(1)*x+u(2)*y+u(3)*z);
        % y' = (v(1)*x+v(2)*y+v(3)*z);
        % z' = (w(1)*x+w(2)*y+w(3)*z);

    end

end


function plotSlice(index1,index2,index3,minBlockVolFrac,cells,sliceDepth,networkType,invertColor)
% Produces an image of the double diamond or double gyroid sliced on the
% given plane. 
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% sliceDepth - The position of the slice within the unit cell, given as a
% fraction or decimal between 0 and 1.
% networkType - The kind of network, either double diamond or double
% gyroid, given as a string.
% invertColor - Logical value (true/false) that inverts the color of the plot. 
% When false, minority networks are white and the majority network is 
% black. When true, minority networks are black and the majority
% network is white.
% Outputs:
% None

    if sliceDepth > 1 || sliceDepth < 0

        error("Slice Depth Must be a Fraction or Decimal With Value in the Range 0 to 1.")

    elseif ~islogical(invertColor)

        error("invertColor must be a logical value (either true or false)." + ...
            " Do not use quotation marks.");

    end

    network = interpretNetworkType(networkType);

    t = calculateT(minBlockVolFrac,network);

    newAxis = [index1;index2;index3];

    [u,v,w] = rotate001toNewAxis(newAxis);

    range = sqrt(index1*index1+index2*index2+index3*index3);

    z = sliceDepth*range;

    if cells >=1

        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;  
    [x,y] = meshgrid(points);  

    if strcmp(network,'diamond')

        value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));

    else

        value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end

    idx = double(abs(value)<t);
    
    figure();
    x0=100;
    y0=70;
    width=750;
    height=700;
    set(gcf,'position',[x0,y0,width,height])
    hold on;
    
    pcolor(x,y,idx);

    if invertColor

        map = [[0,0,0];[1,1,1]];

    else 
    
        map = [[1,1,1];[0,0,0]];

    end

    colormap(map);
    shading interp;

    if strcmp(network,'diamond')

        if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
    
            init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
        
        elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
    
            init = 'Double Diamond Slice on (%d%d%d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        else
    
            init = 'Double Diamond Slice on (%d %d %d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        end
    
    else

        if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))

            init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
    
            init = 'Double Gyroid Slice on (%d%d%d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        else
    
            init = 'Double Gyroid Slice on (%d %d %d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        end

    end

    title(str);
    xlim([0, cells]);
    ylim([0, cells]);

    hold off;

end


function exportSlices(index1,index2,index3,minBlockVolFrac,cells,numSlices,filePath,networkType,invertColor)
% Produces a series of images of the double diamond or double gyroid 
% sliced on the given plane at varying heights in the unit cells, and then 
% directly exports the figures to the specfied folder as '.tiff' images.
% WARNING: This function creates images named 1,2,3...etc. If there are
% already files with these name in the destination folder, they will be
% overwritten. Running the function twice with the same destination will
% lead to the second set of images overwriting the first.
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% numSlices - Number of slices the unit cell is divided into.
% filePath - The file directory where the images will be saved, should be a
% folder.
% networkType - The kind of BCP network, either double diamond or double
% gyroid.
% invertColor - Boolean value that inverts the color of the plot. 
% When false, minority networks are white and the majority network is 
% black. When true, minority networks are black and the majority
% network is white.
% Outputs:
% None

    if ~islogical(invertColor)
    
        error("invertColor must be a logical value (either true or false)." + ...
            " Do not use quotation marks.");

    elseif numSlices < 3

        error("Must use at least 3 slices. For single slices, use plotAndSaveSlice.");
    
    end
    
    network = interpretNetworkType(networkType);
    
    t = calculateT(minBlockVolFrac,network);
    
    newAxis = [index1;index2;index3];
    
    [u,v,w] = rotate001toNewAxis(newAxis);
    
    range = sqrt(index1*index1+index2*index2+index3*index3);
    
    i = 1;
    
    step = range/(numSlices-1);
    
    slices = 0:step:range;
    
    if cells >=1
    
        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;    
    [x,y] = meshgrid(points);  
    
    for z = slices

        if strcmp(network,'diamond')

            value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
            sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
            +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
            .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
            .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
            +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
            .*sin(pi*(w(1)*x+w(2)*y+w(3)*z)); 
    
        else
    
            value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
            sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
            sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));
    
        end
    
        idx = double(abs(value)<t);
        
        figure();
        x0=100;
        y0=70;
        width=750;
        height=700;
        set(gcf,'position',[x0,y0,width,height])
        hold on;
        
        pcolor(x,y,idx);
      
        if invertColor

            map = [[0,0,0];[1,1,1]];
    
        else 
        
            map = [[1,1,1];[0,0,0]];

        end

        colormap(map);
        shading interp;

        if strcmp(network,'diamond')

            if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
        
                init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
            
            elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
        
                init = 'Double Diamond Slice on (%d%d%d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            else
        
                init = 'Double Diamond Slice on (%d %d %d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            end
        
        else
    
            if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
    
                init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
        
                init = 'Double Gyroid Slice on (%d%d%d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            else
        
                init = 'Double Gyroid Slice on (%d %d %d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            end

        end
    
        title(str);
        xlim([0, cells]);
        ylim([0, cells]);
    
        hold off;
    
        imagewd = getframe(gcf); 
        
        specificName = '\\%d.tiff';
        newString = sprintf(specificName,i);
        fileLocation = strcat(filePath,newString);
        imwrite(imagewd.cdata, fileLocation);
    
        i= i+1;
    
        seq = mod(i,5);
    
        if seq == 0
    
            close all;
    
        end
    
    end
    
    close all;    

end


function plotAndSaveSlice(index1,index2,index3,minBlockVolFrac,cells,sliceDepth,filePath,networkType,invertColor)
% Produces an image of the double diamond or double gyroid sliced on the 
% given plane. It then directly exports the figures to the specfied folder 
% as a '.tiff' image. 
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% sliceDepth - The position of the slice within the unit cell, given as a
% fraction or decimal between 0 and 1.
% networkType - The kind of BCP network, either double diamond or double
% gyroid.
% invertColor - Logical value (true/false) that inverts the color of the plot. 
% When false, minority networks are white and the majority network is 
% black. When true, minority networks are black and the majority
% network is white.
% Outputs:
% None

    if sliceDepth > 1 || sliceDepth < 0

        error("Slice Depth Must be a Fraction or Decimal With Value in the Range 0 to 1.")

    elseif ~islogical(invertColor)
    
        error("invertColor must be a logical value (either true or false)." + ...
            " Do not use quotation marks.");

    end

    network = interpretNetworkType(networkType);

    t = calculateT(minBlockVolFrac,network);

    newAxis = [index1;index2;index3];

    [u,v,w] = rotate001toNewAxis(newAxis);

    range = sqrt(index1*index1+index2*index2+index3*index3);

    z = sliceDepth*range;

    if cells >=1

        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;  
    [x,y] = meshgrid(points);  % get 2-D mesh for x and y

    if strcmp(network,'diamond')

        value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));

    else

        value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end
    
    idx = double(abs(value)<t);
    
    figure();
    x0=100;
    y0=70;
    width=750;
    height=700;
    set(gcf,'position',[x0,y0,width,height])
    hold on;
    
    pcolor(x,y,idx);

    if invertColor

        map = [[0,0,0];[1,1,1]];

    else 
    
        map = [[1,1,1];[0,0,0]];

    end

    colormap(map);
    shading interp;

    if strcmp(network,'diamond')

        init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
        str = sprintf(init,index1,index2,index3,z/range);

    else

        init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
        str = sprintf(init,index1,index2,index3,z/range);

    end

    title(str);
    xlim([0, cells]);
    ylim([0, cells]);

    hold off;

    imagewd = getframe(gcf); 
    
    fileType = '.tiff';
    backslash = '\';
    newString = strcat(backslash,str,fileType);
    fileLocation = strcat(filePath,newString);
    imwrite(imagewd.cdata, fileLocation);

end


function plotSliceColorNetworks(index1,index2,index3,minBlockVolFrac,cells,sliceDepth,networkType,minorityColor1,minorityColor2,majorityColor)
% Produces an image of the double diamond or double gyroid sliced on the 
% given plane and colors each network a custom color.
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% sliceDepth - The position of the slice within the unit cell, given as a
% fraction or decimal between 0 and 1.
% networkType - The kind of BCP network, either double diamond or double
% gyroid.
% minorityColor1 - Color for first minority block network, given as an 
% RGB triplet in a 1x3 vector.
% minorityColor2 - Color for second minority block network, given as an 
% RGB triplet in a 1x3 vector.
% majorityColor - Color for the majority block network, given as an 
% RGB triplet in a 1x3 vector.
% NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255]. 
% Outputs:
% None

    if sliceDepth > 1 || sliceDepth < 0

        error("Slice Depth Must be a Fraction or Decimal With Value in the Range 0 to 1.");

    end

    network = interpretNetworkType(networkType);

    t = calculateT(minBlockVolFrac,network);

    newAxis = [index1;index2;index3];

    [u,v,w] = rotate001toNewAxis(newAxis);

    range = sqrt(index1*index1+index2*index2+index3*index3);

    z = sliceDepth*range;

    if cells >=1

        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;  
    [x,y] = meshgrid(points);  % get 2-D mesh for x and y

    if strcmp(network,'diamond')

        value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));

    else

        value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end

    idx = double(value>t);
    idx = idx - double(value<-t);
    
    figure();
    x0=100;
    y0=70;
    width=750;
    height=700;
    set(gcf,'position',[x0,y0,width,height])
    hold on;
    
    pcolor(x,y,idx);

    map = [minorityColor1;majorityColor;minorityColor2];

    colormap(map);
    shading interp;

    if strcmp(network,'diamond')

        if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
    
            init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
        
        elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
    
            init = 'Double Diamond Slice on (%d%d%d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        else
    
            init = 'Double Diamond Slice on (%d %d %d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        end
    
    else

        if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))

            init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
    
            init = 'Double Gyroid Slice on (%d%d%d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        else
    
            init = 'Double Gyroid Slice on (%d %d %d) Plane, Depth = %1.4f';
            str = sprintf(init,index1,index2,index3,z/range);
    
        end

    end

    title(str);
    xlim([0, cells]);
    ylim([0, cells]);

    hold off;

end


function exportSlicesColorNetworks(index1,index2,index3,minBlockVolFrac,cells,numSlices,networkType,minorityColor1,minorityColor2,majorityColor,filePath)
% Produces a series of images of the double diamond or double gyroid
% sliced on the given plane at varying heights in the unit cells and 
% colors each network a custom color. It then directly exports the figures
% to the specfied folder as '.tiff' images.
% WARNING: This function creates images named 1,2,3...etc. If there are
% already files with these name in the destination folder, they will be
% overwritten. Running the function twice with the same destination will
% lead to the second set of images overwriting the first.
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% numSlices - Number of slices the unit cell is divided into.
% networkType - The kind of BCP network, either double diamond or double
% gyroid.
% minorityColor1 - Color for first minority block network, given as an 
% RGB triplet in a 1x3 vector.
% minorityColor2 - Color for second minority block network, given as an 
% RGB triplet in a 1x3 vector.
% majorityColor - Color for the majority block network, given as an 
% RGB triplet in a 1x3 vector.
% NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255]. 
% filePath - The file directory where the images will be saved, should be a
% folder.
% Outputs:
% None

    if numSlices < 3

        error("Must use at least 3 slices. To save single slices, use plotAndSaveSlice.");
    
    end
    
    network = interpretNetworkType(networkType);
    
    t = calculateT(minBlockVolFrac,network);
    
    newAxis = [index1;index2;index3];
    
    [u,v,w] = rotate001toNewAxis(newAxis);
    
    range = sqrt(index1*index1+index2*index2+index3*index3);
    
    i = 1;
    
    step = range/(numSlices-1);
    
    slices = 0:step:range;
    
    if cells >=1
    
        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;    
    [x,y] = meshgrid(points);  % get 2-D mesh for x and y
    
    for z = slices

        if strcmp(network,'diamond')

            value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
            sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
            +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
            .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
            .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
            +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
            .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));
    
        else
    
            value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
            sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
            sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));
    
        end
    
         idx = double(value>t);
         idx = idx - double(value<-t);
        
        figure();
        x0=100;
        y0=70;
        width=750;
        height=700;
        set(gcf,'position',[x0,y0,width,height])
        hold on;
        
        pcolor(x,y,idx);
      
        map = [minorityColor1;majorityColor;minorityColor2];

        colormap(map);
        shading interp;

        if strcmp(network,'diamond')

            if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
        
                init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
            
            elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
        
                init = 'Double Diamond Slice on (%d%d%d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            else
        
                init = 'Double Diamond Slice on (%d %d %d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            end
        
        else
    
            if (~(floor(index1)==index1 && floor(index2)==index2 && floor(index3)==index3))
    
                init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            elseif ((index1>=0) && (index2>=0) && (index3>=0) && (index1<10) && (index2<10) && (index3<10))
        
                init = 'Double Gyroid Slice on (%d%d%d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            else
        
                init = 'Double Gyroid Slice on (%d %d %d) Plane, Depth = %1.4f';
                str = sprintf(init,index1,index2,index3,z/range);
        
            end

        end
    
        title(str);
        xlim([0, cells]);
        ylim([0, cells]);
    
        hold off;
    
        imagewd = getframe(gcf); 
        
        specificName = '\\%d.tiff';
        newString = sprintf(specificName,i);
        fileLocation = strcat(filePath,newString);
        imwrite(imagewd.cdata, fileLocation);
    
        i= i+1;
    
        seq = mod(i,5);
    
        if seq == 0
    
            close all;
    
        end
    
    end
    
    close all;    

end


function plotAndSaveSliceColorNetworks(index1,index2,index3,minBlockVolFrac,cells,sliceDepth,networkType,minorityColor1,minorityColor2,majorityColor,filePath)
% Produces an image of the double gyroid or double diamond sliced on the 
% given plane, and colors each network a custom color. It then directly 
% exports the figure to the specfied folder as a '.tiff' image.
% Inputs:
% index1 - The first index (Miller Index) of the plane to be sliced.
% index2 - The second index (Miller Index) of the plane to be sliced.
% index3 - The third index (Miller Index) of the plane to be sliced.
% minBlockVolFrac - The volume fraction of the minority block.
% cells - The number of unit cells to be plotted in the x and y directions.
% sliceDepth - The position of the slice within the unit cell, given as a
% fraction or decimal between 0 and 1.
% networkType - The kind of BCP network, either double diamond or double
% gyroid.
% minorityColor1 - Color for first minority block network, given as an 
% RGB triplet in a 1x3 vector.
% minorityColor2 - Color for second minority block network, given as an 
% RGB triplet in a 1x3 vector.
% majorityColor - Color for the majority block network, given as an 
% RGB triplet in a 1x3 vector.
% NOTE: In Matlab, R,G, and B must be in the range [0,1], not [0,255]. 
% filePath - The file directory where the images will be saved, should be a
% folder.
% Outputs:
% None

    if sliceDepth > 1 || sliceDepth < 0

        error("Slice Depth Must be a Fraction or Decimal With Value in the Range 0 to 1.")

    end

    network = interpretNetworkType(networkType);

    t = calculateT(minBlockVolFrac,network);

    newAxis = [index1;index2;index3];

    [u,v,w] = rotate001toNewAxis(newAxis);

    range = sqrt(index1*index1+index2*index2+index3*index3);

    z = sliceDepth*range;

    if cells >=1

        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end
    
    points = 0:resolution:cells;  
    [x,y] = meshgrid(points);  % get 2-D mesh for x and y

    if strcmp(network,'diamond')

        value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z)); 

    else

        value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end
    
    idx = double(value<t);
    idx = idx - double(value<-t);
    
    figure();
    x0=100;
    y0=70;
    width=750;
    height=700;
    set(gcf,'position',[x0,y0,width,height])
    hold on;
    
    pcolor(x,y,idx);

    map = [minorityColor1;majorityColor;minorityColor2];

    colormap(map);
    shading interp;

    if strcmp(network,'diamond')

        init = 'Double Diamond Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
        str = sprintf(init,index1,index2,index3,z/range);

    else

        init = 'Double Gyroid Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f';
        str = sprintf(init,index1,index2,index3,z/range);

    end

    title(str);
    xlim([0, cells]);
    ylim([0, cells]);

    hold off;

    imagewd = getframe(gcf); 
    
    fileType = '.tiff';
    backslash = '\';
    newString = strcat(backslash,str,fileType);
    fileLocation = strcat(filePath,newString);
    imwrite(imagewd.cdata, fileLocation);

end


function twoPlaneSlice(plane1Index1,plane1Index2,plane1Index3,sliceDepth1,networkType1,minBlockVolFrac1,plane2Index1,plane2Index2,plane2Index3,sliceDepth2,networkType2,minBlockVolFrac2,cells,xCoefficient,b,invertColor1,invertColor2)
% Produces an image of the double gyroid or double diamond sliced on two
% given planes, separated by a linear grain boundary. These two planes can
% be of the same network type (ie. both double diamond) or different 
% network types.
% Inputs:
% plane1Index1 - The first index (Miller Index) of the first plane to be sliced.
% plane1Index2 - The second index (Miller Index) of the first plane to be sliced.
% plane1Index3 - The third index (Miller Index) of the first plane to be sliced.
% sliceDepth1 - The position of the first slice within the unit cell, 
% given as a fraction or decimal between 0 and 1.
% networkType1 - The kind of network for the first network, either double 
% diamond or double gyroid.
% minBlockVolFrac1 - The volume fraction of the minority block for the 
% first network.
% plane2Index1 - The first index (Miller Index) of the second plane to be sliced.
% plane2Index2 - The second index (Miller Index) of the second plane to be sliced.
% plane2Index3 - The third index (Miller Index) of the second plane to be sliced.
% sliceDepth2 - The position of the second slice within the unit cell, 
% given as a fraction or decimal between 0 and 1.
% networkType2 - The kind of network for the second network, either double 
% diamond or double gyroid.
% minBlockVolFrac2 - The volume fraction of the minority block for the 
% second network.
% cells - The total number of unit cells to be plotted in the x and y 
% directions.
% xCoefficient - The grain boundary separating the two planes is a line
% given by the equation y = m*x+b. This variable is eqivalent to m.
% b - The grain boundary separating the two planes is a line given by the 
% equation y = m*x+b. This variable is eqivalent to b.
% invertColor1 - Logical value (true/false) that inverts the color of the 
% first network. When false, minority networks are white and the majority 
% network is black. When true, minority networks are black and the majority
% network is white.
% invertColor2 - Logical value (true/false) that inverts the color of the 
% second network. When false, minority networks are white and the majority 
% network is black. When true, minority networks are black and the majority
% network is white.
% Outputs:
% None

    if sliceDepth1 > 1 || sliceDepth1 < 0 || sliceDepth2 > 1 || sliceDepth2 < 0
    
        error("Slice Depth Must be a Fraction or Decimal With Value in the Range 0 to 1.")

    elseif ~(islogical(invertColor1) && islogical(invertColor2))

        error("invertColor must be a logical value (either true or false)." + ...
            " Do not use quotation marks.");
    
    end

    network1 = interpretNetworkType(networkType1);
    network2 = interpretNetworkType(networkType2);

    tPlane1 = calculateT(minBlockVolFrac1,network1);
    tPlane2 = calculateT(minBlockVolFrac2,network2);

    if cells >=1

        resolution = round(cells)/2000;
    
    else
    
        error("Invalid Number of Cells, You Must Plot at Least One Cell.");
    
    end

    figure();
    x0=100;
    y0=70;
    width=750;
    height=700;
    set(gcf,'position',[x0,y0,width,height])
    hold on;
    
    points = 0:resolution:cells;
    [x,y] = meshgrid(points);  % get 2-D mesh for x and y

    xMinusY = xCoefficient.*x-y; % y = A*x+b
    [rowIndices1, colIndices1] = find(xMinusY > -b); % y < A*x+b
    indicies1 = [rowIndices1, colIndices1];
    [rowIndices2, colIndices2] = find(xMinusY <= -b); % y >= A*x+b
    indicies2 = [rowIndices2,colIndices2];

    newAxis = [plane1Index1;plane1Index2;plane1Index3];

    [u,v,w] = rotate001toNewAxis(newAxis);

    range = sqrt(plane1Index1*plane1Index1+plane1Index2*plane1Index2+plane1Index3*plane1Index3);

    z = sliceDepth1*range;

    if strcmp(network1,'diamond')

        value = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));

    else

        value = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end 

    idx = double(abs(value)<tPlane1);

    for row = 1:size(indicies1,1)

        idx(indicies1(row,1),indicies1(row,2)) = NaN; %set to NaN so it is not plotted

    end
    
    pcolor(x,y,idx);
    
    if invertColor1

        map = [[0,0,0];[1,1,1]];

    else 
    
        map = [[1,1,1];[0,0,0]];

    end

    colormap(map);
    shading interp;

    newAxis2 = [plane2Index1;plane2Index2;plane2Index3];

    [u,v,w] = rotate001toNewAxis(newAxis2);

    range2 = sqrt(plane2Index1*plane2Index1+plane2Index2*plane2Index2+plane2Index3*plane2Index3);

    z = sliceDepth2*range2;

    if strcmp(network2,'diamond')

        value2 = sin(pi*(u(1)*x+u(2)*y+u(3)*z)).* ...
        sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*sin(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +sin(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*cos(pi*(w(1)*x+w(2)*y+w(3)*z))+cos(pi*(u(1)*x+u(2)*y+u(3)*z)) ...
        .*sin(pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(pi*(w(1)*x+w(2)*y+w(3)*z)) ...
        +cos(pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(pi*(v(1)*x+v(2)*y+v(3)*z)) ...
        .*sin(pi*(w(1)*x+w(2)*y+w(3)*z));

    else

        value2 = sin(2*pi*(u(1)*x+u(2)*y+u(3)*z)).*cos(2*pi*(v(1)*x+v(2)*y+v(3)*z)) + ...
        sin(2*pi*(v(1)*x+v(2)*y+v(3)*z)).*cos(2*pi*(w(1)*x+w(2)*y+w(3)*z))+ ...
        sin(2*pi*(w(1)*x+w(2)*y+w(3)*z)).*cos(2*pi*(u(1)*x+u(2)*y+u(3)*z));

    end 

    idx = double(abs(value2)<tPlane2);

    for row = 1:size(indicies2,1)

        idx(indicies2(row,1),indicies2(row,2)) = NaN; %set to NaN so it is not plotted

    end
    
    pcolor(x,y,idx);
    
    if invertColor2

        map = [[0,0,0];[1,1,1]];

    else 
    
        map = [[1,1,1];[0,0,0]];

    end

    colormap(map);
    shading interp;

    if strcmp(network1,'diamond')

        type1 = 'Diamond';

    else 

        type1 = "Gyroid";

    end

    if strcmp(network2,'diamond')

        type2 = 'Diamond';

    else 

        type2 = "Gyroid";

    end

    init = 'Double %s Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f, and Double %s Slice on (%1.2f %1.2f %1.2f) Plane, Depth = %1.4f.';
    str = sprintf(init,type1,plane1Index1,plane1Index2,plane1Index3,sliceDepth1,type2,plane2Index1,plane2Index2,plane2Index3,sliceDepth2);

    title(str);
    ax = gca;
    ax.TitleFontSizeMultiplier = .6;
    xlim([0, cells]);
    ylim([0, cells]);

    hold off;

end
