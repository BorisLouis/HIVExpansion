classdef HIVCellMovie < Core.HIVLocMovie
    properties (SetAccess = 'protected')
        lamina
        nup
        lipid
        
                
    end

    methods
        
        function obj = HIVCellMovie(ext,info,varargin)
            
            obj  = obj@Core.HIVLocMovie(ext,info,varargin{1});
            
        end
        
        function [finalBW] = segmentLamina(obj)
            
            idx = contains({obj.raw.movInfo.fileName},'lamina','IgnoreCase',true);
            if sum(idx)~=1
                warning('no lamina file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'lamina');
                % Gaussian filtering
                S = 2;
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                 if size(data,3)~=1
                    sigma = [S,S,S/zFactor];
                    IMs = imgaussfilt3(data, sigma);
                    neigh = [51 51 3];
                else
                    IMs = imgaussfilt(data,S);
                    neigh = [51 51];
                end
                
                disp('DONE with filtering ------------')
                threshold = adaptthresh(uint16(IMs),0.4,'neigh',neigh);
                gBW = imbinarize(uint16(IMs),threshold);
             %   gBW = imbinarize(uint16(IMs));
                gBW = bwareaopen(gBW,10000);
                gBW = imclearborder(gBW);
                se = strel('disk',10);
                gBW = imclose(gBW,se);
                
                %gBW = imfill(gBW,'holes');
                disp('Extracting Contour')
                
                BWdata = regionprops(gBW,{'Area','Perimeter', 'pixelIDXList'});
                
                [~,maxArea] = max([BWdata.Perimeter]);
                
                finalBW = zeros(size(gBW));
                finalBW(BWdata(maxArea).PixelIdxList) = 1;
                
                %TODO BW CHECK
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.lamina.mask = finalBW;
                obj.lamina.contour = contour;
                %
            end
        end
        
        
        function [finalBW] = segmentNUP(obj)
            
            idx = contains({obj.raw.movInfo.fileName},'nup','IgnoreCase',true);
            if sum(idx)~=1
                warning('no NUP file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'nup');
                % Gaussian filtering
                S = 5;
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                sigma = [S,S,S/zFactor];
                if size(data,3)~=1
                    sigma = [S,S,S/zFactor];
                    IMs = imgaussfilt3(data, sigma);
                    neigh = [51 51 3];
                else
                    IMs = imgaussfilt(data,S);
                    neigh = [101 101];
                end
                
                disp('DONE with filtering ------------')
                threshold = adaptthresh(uint16(IMs),0.6,'neigh',neigh);
                gBW = imbinarize(uint16(IMs),threshold);
                %gBW = imbinarize(uint16(IMs),'adaptive','Sensitivity',0.4);
                %  gBW = imbinarize(uint16(IMs));
                gBW = bwareaopen(gBW,10000);
                gBW = imclearborder(gBW);
                se = strel('disk',10);
                gBW = imclose(gBW,se);
                
                BWdata = regionprops(gBW,{'Area','Perimeter','pixelIDXList'});
                
                [~,maxArea] = max([BWdata.Perimeter]);
                
                finalBW = zeros(size(gBW));
                finalBW(BWdata(maxArea).PixelIdxList) = 1;
                
                se = strel('disk',25);
                finalBW = imclose(finalBW,se);
                
                %gBW = imfill(gBW,'holes');
                disp('Extracting Contour')
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.nup.mask = finalBW;
                obj.nup.contour = contour;
            end
        end
        
        function [finalBW] = segmentLipid(obj)
            idx = contains({obj.raw.movInfo.fileName},'lipid','IgnoreCase',true);
            
            if sum(idx)~=1
                warning('no lipid file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'lipid');
                
                % Gaussian filtering
                S = 2;
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                sigma = [S,S,S/zFactor];
                IMs = imgaussfilt3(data, sigma);
                disp('DONE with filtering ------------')
                
                gBW = imbinarize(uint16(IMs),'adaptive','Sensitivity',0.8);
                gBW = bwareaopen(gBW,10000);
                se = strel('disk',20);
                gBW = imclose(gBW,se);
                
                % invert
                invertBW = imcomplement(gBW);
                invertBW = bwareaopen(invertBW,20000);
                
                for i = 1:size(invertBW,3)
                    invertBW(:,:,i) = imfill(invertBW(:,:,i),'holes');
                end
                
                % find the biggest round object common in middle plane
                BWdata = regionprops(invertBW(:,:,1),{'Area','PixelIdxList','MajorAxisLength','MinorAxisLength'});
                
                for i =1:length(BWdata)
                    currData = BWdata(i);
                    
                    calcR = (currData.MajorAxisLength+currData.MinorAxisLength)/2;
                    
                    expectedDiskArea = pi*(calcR/2)^2;
                    
                    BWdata(i).roundness = BWdata(i).Area/expectedDiskArea;
                    
                    
                end
                
                roundBW = BWdata([BWdata.roundness]>0.85);
                
                biggestRoundBW = roundBW([roundBW.Area] == max([roundBW.Area]));
                
                volBWData = regionprops3(invertBW,'VoxelIdxList');
                finalIDX = 0;
                for i = 1:size(volBWData,1)
                    currIDX = volBWData.VoxelIdxList{i};
                    
                    if all(ismember(biggestRoundBW.PixelIdxList,currIDX))
                        finalIDX = currIDX;
                        break;
                    end
                end
                
                
                if finalIDX ~=0
                else
                    error('Something went wrong in the calculation')
                end
                
                disp('Extracting Contour')
                
                finalBW = zeros(size(gBW));
                finalBW(finalIDX) = 1;
                
                for i = 1:size(finalBW,3)
                    %Need to make a second make 2pixel bigger than the first
                    se = strel('disk',20);
                    biggerMask = imdilate(finalBW(:,:,i),se);
                    
                    finalBW(:,:,i)  =  biggerMask-finalBW(:,:,i);
                    
                    
                end
                
                %TODO BW CHECK
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.lipid.mask = finalBW;
                obj.lipid.contour = contour;
            end
            
            
        end
        
        function [finalBW] = segmentRedLipid(obj)
            idx = contains({obj.raw.movInfo.fileName},'lipid','IgnoreCase',true);
            
            if sum(idx)~=1
                warning('no lipid file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'lipid');
                
                % Gaussian filtering
                S = 2;
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                sigma = [S,S,S/zFactor];
                IMs = imgaussfilt3(data, sigma);
                disp('DONE with filtering ------------')
                threshold = adaptthresh(uint16(IMs),0.8,'neigh',[301 301 3]);
                gBW = imbinarize(uint16(IMs),threshold);
                gBW = bwareaopen(gBW,10000);
%                 se = strel('disk',10);
%                 gBW = imclose(gBW,se);
                
                % invert
                invertBW = imcomplement(gBW);
                invertBW = bwareaopen(invertBW,20000);
                
                for i = 1:size(invertBW,3)
                    invertBW(:,:,i) = imclearborder(invertBW(:,:,i));
                    invertBW(:,:,i) = imfill(invertBW(:,:,i),'holes');
                end
                
                % find the biggest round object common in middle plane
                BWdata = regionprops(invertBW(:,:,1),{'Area','PixelIdxList','MajorAxisLength','MinorAxisLength'});
                
                for i =1:length(BWdata)
                    currData = BWdata(i);
                    
                    calcR = (currData.MajorAxisLength+currData.MinorAxisLength)/2;
                    
                    expectedDiskArea = pi*(calcR/2)^2;
                    
                    BWdata(i).roundness = BWdata(i).Area/expectedDiskArea;
                    
                    
                end
                
                roundBW = BWdata([BWdata.roundness]>0.85);
                
                biggestRoundBW = roundBW([roundBW.Area] == max([roundBW.Area]));
                
                volBWData = regionprops3(invertBW,'VoxelIdxList');
                finalIDX = 0;
                for i = 1:size(volBWData,1)
                    currIDX = volBWData.VoxelIdxList{i};
                    
                    if all(ismember(biggestRoundBW.PixelIdxList,currIDX))
                        finalIDX = currIDX;
                        break;
                    end
                end
                
                
                if finalIDX ~=0
                else
                    error('Something went wrong in the calculation')
                end
                
                disp('Extracting Contour')
                
                finalBW = zeros(size(gBW));
                finalBW(finalIDX) = 1;
                
                for i = 1:size(finalBW,3)
                    %Need to make a second make 2pixel bigger than the first
                    se = strel('disk',20);
                    biggerMask = imdilate(finalBW(:,:,i),se);
                    
                    finalBW(:,:,i)  =  biggerMask-finalBW(:,:,i);
                    
                end
                
                %TODO BW CHECK
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.lipid.mask = finalBW;
                obj.lipid.contour = contour;
            end      
        end
        
        function showMembrane(obj,membrane,idx)
            %get membrane data
            data = obj.getFrame(1,membrane);
            [~,contour,fitMembrane] = getMembrane(obj,membrane);
            
            if nargin<3
                idx = round(size(data,3)/2);
            end    

            if ~isempty(data)
                %get middle slice
                
                figure
                hold on
                imagesc(data(:,:,idx))
                colormap('hot')
                plot(contour{idx}{1}(:,2),contour{idx}{1}(:,1),'w','Linewidth',2);
                plot(contour{idx}{2}(:,2),contour{idx}{2}(:,1),'w','Linewidth',2);
                if ~isempty(fitMembrane)
                    plot(fitMembrane{idx}(:,2),fitMembrane{idx}(:,1),'g','Linewidth',2);
                end
                
                title([membrane ' staining'])
                xlabel('Pixels')
                ylabel('Pixels')
                axis image
                
            end
    
        end
        
        function [membrane] = getMembranePos(obj,membrane)
            %get membrane data
            data = obj.getFrame(1,membrane);
            [~,contour,~] = getMembrane(obj,membrane);
            
            if ~isempty(data)
                membranePos = cell(size(contour));
                for i = 1:length(contour)
                    %get small and large contour
                    currContour = contour{i};
                    idx = cellfun(@length,currContour);
                    sContour = currContour{idx==min(idx)};
                    lContour = currContour{idx==max(idx)};

                    currData = data(:,:,i);
                    dim = size(currData);
                    fContour = zeros(size(sContour));
                    %loop through the smallest contour
                    for j = 1: length(sContour)
                        currPoint = sContour(j,:);
                       
                        %find closest point in the largest lContour
                        closestList = sqrt((lContour(:,1)-currPoint(:,1)).^2 +...
                            ((lContour(:,2)-currPoint(:,2)).^2));
                        [~,idx] = min(closestList);
                        closestPoint = lContour(idx,:);
                        %create x-y points in order to make a line of pixel
                        x = [currPoint(1,2); closestPoint(1,2)];
                        y = [currPoint(1,1); closestPoint(1,1)];
                        
                        center = [mean(x),mean(y)];
                        %get the slope of the line
                        m = (y(2)-y(1))/(x(2)-x(1));
                        %get the origin of the line
                        b = y(1)-m*x(1);
                        %sort data point base on x axis
                        [x,idx] = sort(x);
                        %store data about the line to be used in next
                        %function
                        linEq.m = m;
                        linEq.b = b;
                        linEq.x = x;
                        linEq.y = y(idx);
                        %get line of pixel 
                        [idx,xVec,yVec] = obj.getLinePixel(linEq,dim);

                        pxLine = currData(idx);
                        %get minimum and maximum of the vector
                        %check which coordinate has the large differential
                        %to be used as domain for fitting
                        if abs(mean(diff(xVec)))< abs(mean(diff(yVec)))
                            vec2Use = yVec;
                            minMax  = [y(1) y(2)];
                        else
                            vec2Use = xVec;
                            minMax  = [x(1) x(2)];close all
                        end
           
                        guess.sig = 10;
                        guess.mu  = vec2Use(round(length(pxLine)/2));
                        guess.minMaxDomain = minMax;

                        try
                            [gPar, Fit,res] = SimpleFitting.gauss1D(pxLine,vec2Use,guess);
                            mse(1,j) = gPar(2) - median(vec2Use);
                            mse(2,j) =max(pxLine)/std(pxLine); 
                        catch
                        end
%                                           figure
%                                           plot(pxLine)
%                                           hold on
%                                           plot(Fit)
%                                           clf

                        if abs(mean(diff(xVec))) < abs(mean(diff(yVec)))
                            %then we get y position and x constant
                            fContour(j,1) = gPar(2);
                            %check that m is not infinite (x == constant)
                            if abs(m)~=inf
                                fContour(j,2) = (gPar(2)-b)/m;
                            else
                                fContour(j,2) = xVec(1);
                            end
                        else
                            %then x position from fit and y from equation
                            fContour(j,2) = gPar(2);
                            fContour(j,1) = m*gPar(2)+b;

                        end
                        if any(isnan(fContour(:)))
                            error('Something went wrong, need to check the count');
                        end

                    end
                    membranePos{i} = [smooth(fContour(:,1)) smooth(fContour(:,2))];
                end

            else
                membranePos = [];
            end
            
            obj.(lower(membrane)).membranePos = membranePos;
            
        end      
        
        function [mask,contour,membranePos] = getMembrane(obj,membrane)
            idx = contains({obj.raw.movInfo.fileName},membrane,'IgnoreCase',true);
            
            if sum(idx)~=1
                warning(['no ' membrane ' file was found, operation aborted']);
      
                mask = [];
                contour = [];
                
            else                
                
                var = lower(membrane);
                contour = obj.(var).contour;
                mask    = obj.(var).mask;
                if isfield(obj.(var),'membranePos')
                    membranePos = obj.(var).membranePos;
                else
                    membranePos = [];
                end
                
            end
        end
        
        function [partPos] = getHIVToMembraneDistance(obj,membrane)
            
            [~,~,memPos] = obj.getMembrane(membrane);
            partPos = obj.particles.SRList;
            pxSize = obj.info.pxSize;
            zPos = obj.raw.movInfo.planePos;
            zSpacing = obj.raw.movInfo.zSpacing;
            
            for i = 1:height(partPos)
                currPart = partPos(i,:);
                %check which membrane to use
                idx2Membrane = ceil(currPart.z/zSpacing);
                
                %get the membrane
                currMembrane = memPos{idx2Membrane}*pxSize;
                
                %get sign of distance to membrane by looking if particle is
                %inside or outside of it, inside will be positive, outside
                %will be negative for that, we compare the distances with
                %the center of mass
                
                CM = mean(currMembrane,1);
                %find the closest point on the membrane to the current
                %particle
                [~,idx] = min(sqrt((currPart.col-currMembrane(:,1)).^2 +(currPart.row-currMembrane(:,2)).^2));
                
                %get the distance between the current point and the CM
                dist2CM = sqrt((currPart.col-CM(1)).^2 +(currPart.row-CM(2)).^2);
                %get the distance between the closest point on the membrane
                %and the CM
                distMem2CM = sqrt((currMembrane(idx,1)-CM(1)).^2 +(currMembrane(idx,2)-CM(2)).^2);
                %substract one by the other so if the particle is closer to
                %the CM than the membrane it gives a positive distance and
                %a negative otherwise.
                
                partPos.dist2Mem(i) = distMem2CM-dist2CM;

%                 minXm = min(currMembrane(:,1));
%                 maxXm = max(currMembrane(:,1));
%                 minYm = min(currMembrane(:,2));
%                 maxYm = max(currMembrane(:,2));
%                 
%                 if and(and(minXm<currPart.col,maxXm>currPart.col),and(minYm<currPart.row,maxYm>currPart.row))
%                     mult = 1;
%                 else
%                     mult = -1;
%                 end
%  
%                partPos.dist2Mem(i) = mult*min(sqrt((currPart.col-currMembrane(:,1)).^2 +(currPart.row-currMembrane(:,2)).^2));

            end
            
            
            
            
        end
        
        function [distances] = getMembraneToMembraneDistance(obj,membrane1,membrane2)
            
            mem1 = obj.(lower(membrane1)).membranePos;
            mem2 = obj.(lower(membrane2)).membranePos;
            distances = cell(size(mem1));
            pxSize = obj.info.pxSize;
            for i = 1:length(mem1)
               %change from pixels to nanometer
               currMem1 = mem1{i}*pxSize;
               currMem2 = mem2{i}*pxSize;
               
               distances{i} = zeros(length(currMem1),1);
               
               for j =1:length(currMem1)
                   currPoint = currMem1(j,:);
                  
                   CM = mean(currMem2,1);
                   %find the closest point on the membrane to the current
                   %particle
                   [~,idx] = min(sqrt((currPoint(:,1)-currMem2(:,1)).^2 +(currPoint(:,2)-currMem2(:,2)).^2));

                   %get the distance between the current point and the CM
                   distMem12CM = sqrt((currPoint(:,1)-CM(1)).^2 +(currPoint(:,2)-CM(2)).^2);
                   %get the distance between the closest point on the membrane
                   %and the CM
                   distMem22CM = sqrt((currMem2(idx,1)-CM(1)).^2 +(currMem2(idx,2)-CM(2)).^2);
                   
%                  minXm = min(currMem2(:,1));
%                  maxXm = max(currMem2(:,1));
%                  minYm = min(currMem2(:,2));
%                  maxYm = max(currMem2(:,2));
%                 
%                     if and(and(minXm<currPoint(:,1),maxXm>currPoint(:,1)),and(minYm<currPoint(:,2),maxYm>currPoint(:,2)))
%                         mult = 1;
%                     else
%                         mult = -1;
%                     end
                distances{i}(j) = distMem22CM-distMem12CM;
                
               end
                
            end
            
            
            
            
            
        end
               
        function plotContour3D(obj,frame)
            
            contourTmp = obj.contour{frame};
            figure(2)
            hold on
            for i=1:length(contourTmp)
                if ~isempty(contourTmp{i})
                    for j=1:length(contourTmp{i})
                        plot3(contourTmp{i}{j}(:,2)*obj.info.pxSize/1000,...
                            contourTmp{i}{j}(:,1)*obj.info.pxSize/1000,...
                            repmat(i*obj.info.dZ/1000,1,...
                            length(contourTmp{i}{j}(:,1))),'k')
                    end
                end
            end
            
        end
        
        function plotCellContour(obj,frame,z)
            
            data = obj.getFrame(frame,'cells');
            contourTmp = obj.contour{frame};
            cCont = contourTmp{z};
            
            figure
            imagesc(data(:,:,z));
            
            hold on
            for i = 1:length(cCont)
                plot(cCont{i}(:,2),cCont{i}(:,1),'r','LineWidth',2)
            end
            
            
        end
        
        function renderCell3D(obj,frame)
            %compile c code for smoothing
            mex +rendering3D\smoothpatch_curvature_double.c -v
            mex +rendering3D\smoothpatch_inversedistance_double.c -v
            mex +rendering3D\vertex_neighbours_double.c -v
            
            file2Test = [obj.raw.movInfo.path2Cell filesep 'mask.mat'];
            assert(isfile(file2Test),'no mask was found, please run getCellMask first')
            
            allData = load(file2Test);
            field = fieldnames(allData);
            allData = allData.(field{1});
            
            data2Render = allData(:,:,:,frame);
            iSurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
            smoothISurface = rendering3D.smoothpatch(iSurface,0,10);
            %comnvert to px
            smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1));
            smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2));
            smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3));
            
            %% Displaying network model
            %z-coloring
            colorModel = smoothISurface.vertices(:,3)/max(smoothISurface.vertices(:,3));
            zColor = true;
            
            %Plot the network with Z coloring or unique color depending on the user
            %input
            figure
            if zColor
                p = patch('Faces',smoothISurface.faces,'Vertices',smoothISurface.vertices,'FaceVertexCData',colorModel,'FaceColor','interp');
                colormap('jet')
                p.EdgeColor = 'none';
                daspect([2 2 1])
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('Z-coloring')
            else
                p2 = patch(smoothISurface);
                p2.FaceColor = [0.4 0.4 0.4];
                p2.EdgeColor = 'none';
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('unicolor');
            end
            
            
            
        end
        
        
    end
    
    methods (Access = private)
        
        function [contour] = getContour(~,gBW)
            contour = cell(1,size(gBW,3));
            for i = 1:size(gBW,3)
                currBW = gBW(:,:,i);
                
                %  newBW = imfill(currBW,'holes');
                
                %Get the largest area
                %  cBWarea = regionprops(currBW,'Area');
                
                %Get the largest area
                % nBWarea = regionprops(newBW,'Area');
                
                % get index to hollow object (area changed after imfill)
                %idx = ismember([nBWarea.Area],[cBWarea.Area]);
                
                %kill all the other area found
                [pContour] = bwboundaries(currBW);
                contour{i} = pContour;
                
            end
        end
        function [idx,xVec,yVec] = getLinePixel(~,linEq,dim)
            m = linEq.m;
            b = linEq.b;
            x = linEq.x;
            y = linEq.y;
            
            nPixel = max([abs(diff(x)); abs(diff(y))]);
            
            pt1 = [x(1) y(1)];
            pt2 = [x(2) y(2)];
            %get vector between the two points
            V = pt2-pt1;
            
            factor_distance = 4;
            
            if diff(x)~=0
                %recalculate position of vector
                x(1) = pt1(1) - abs(V(1));
                x(2) = pt2(1) + abs(V(1));
                
                xVec = linspace(x(1),x(2),nPixel*(factor_distance+1)+1);
                xVec = unique(xVec);
                yVec = m*xVec+b;
                
            else
                yb = sort(y);
                
                yb(1) = yb(1) - abs(V(2));
                yb(2) = yb(2) + abs(V(2));
                yVec = linspace(yb(1),yb(2),nPixel*(factor_distance+1)+1);
                yVec = unique(yVec);
                if yb == y
                    
                else
                    yVec = fliplr(yVec);
                end
                
                xVec = ones(size(yVec))*x(1);
                
            end
            
            xVec = round(xVec);
            yVec = round(yVec);
            %get the index of the pixel in the image coordinates
            idx = sub2ind(dim,yVec,xVec);
            
        end
    end
end