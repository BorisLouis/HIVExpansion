classdef HIVCellMovie < Core.HIVLocMovie
    properties (SetAccess = 'protected')
        Lamina
        NUP
        Lipid
        
        
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
                sigma = [S,S,S/zFactor];
                IMs = imgaussfilt3(data, sigma);
                disp('DONE with filtering ------------')
                
                %  gBW = imbinarize(uint16(IMs),'adaptive','Sensitivity',threshold);
                gBW = imbinarize(uint16(IMs));
                gBW = bwareaopen(gBW,10000);
                se = strel('disk',20);
                gBW = imclose(gBW,se);
                
                %gBW = imfill(gBW,'holes');
                disp('Extracting Contour')
                
                BWdata = regionprops(gBW,{'Area', 'pixelIDXList'});
                
                [~,maxArea] = max([BWdata.Area]);
                
                finalBW = zeros(size(gBW));
                finalBW(BWdata(maxArea).PixelIdxList) = 1;
                
                %TODO BW CHECK
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.Lamina.mask = finalBW;
                obj.Lamina.contour = contour;
                %
            end
        end
        
        
        function [finalBW] = segmentNUP(obj)
            
            idx = contains({obj.raw.movInfo.fileName},'NUP','IgnoreCase',true);
            if sum(idx)~=1
                warning('no NUP file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'NUP');
                % Gaussian filtering
                S = 5;
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                sigma = [S,S,S/zFactor];
                IMs = imgaussfilt3(data, sigma);
                disp('DONE with filtering ------------')
                
                gBW = imbinarize(uint16(IMs),'adaptive','Sensitivity',0.2);
                %  gBW = imbinarize(uint16(IMs));
                gBW = bwareaopen(gBW,20000);
                
                se = strel('disk',10);
                gBW = imclose(gBW,se);
                
                BWdata = regionprops(gBW,{'Area', 'pixelIDXList'});
                
                [~,maxArea] = max([BWdata.Area]);
                
                finalBW = zeros(size(gBW));
                finalBW(BWdata(maxArea).PixelIdxList) = 1;
                
                se = strel('disk',25);
                finalBW = imclose(finalBW,se);
                
                %gBW = imfill(gBW,'holes');
                disp('Extracting Contour')
                
                %here we obtain the cell contour
                contour = obj.getContour(finalBW);
                
                obj.NUP.mask = finalBW;
                obj.NUP.contour = contour;
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
                
                obj.Lipid.mask = finalBW;
                obj.Lipid.contour = contour;
            end
            
            
        end
        
        function showMembrane(obj,membrane,idx)
            %get membrane data
            [data,~,contour] = getMembrane(obj,membrane);
            
            
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
                title([membrane ' staining'])
                xlabel('Pixels')
                ylabel('Pixels')
                axis image
                
            end
            
            
        end
        
        function [membrane] = getMembranePos(obj,membrane)
            %get membrane data
            [data,~,contour] = getMembrane(obj,membrane);
            
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

                        x = [currPoint(1,2); closestPoint(1,2)];
                        y = [currPoint(1,1); closestPoint(1,1)];

                        center = [mean(x),mean(y)];

                        m = (y(2)-y(1))/(x(2)-x(1));

                        b = y(1)-m*x(1);

                        [x,idx] = sort(x);

                        linEq.m = m;
                        linEq.b = b;
                        linEq.x = x;
                        linEq.y = y(idx);

                        [idx,xVec,yVec] = obj.getLinePixel(linEq,dim);

                        pxLine = currData(idx);

                        if xVec(1)-xVec(2) ==0
                            vec2Use = yVec;
                        else
                            vec2Use = xVec;
                        end

                        [~,id] = max(pxLine);
                        guess.sig = 10;
                        guess.mu  = vec2Use(id);



                        [gPar, Fit] = SimpleFitting.gauss1D(pxLine,vec2Use,guess);

                        %                   figure
                        %                   plot(pxLine)
                        %                   hold on
                        %                   plot(Fit)

                        if xVec(1)-xVec(2) ==0
                            %then we get y position and x constant
                            fContour(j,1) = gPar(2);
                            fContour(j,2) = xVec(1);

                        else
                            %then x position from fit and y from equation
                            fContour(j,2) = gPar(2);
                            fContour(j,1) = m*gPar(2)+b;

                        end

                    end
                    membranePos{i} = fContour;
                end

            else
                membranePos = [];
            end
            
            switch lower(membrane)
                case 'lamina'
                    obj.Lamina.membranePos = membranePos;
                case 'nup'
                    obj.NUP.membranePos = membranePos;
                case 'lipid'
                    obj.Lipid.membranePos = membranePos;
            end
            
        end
        
        
        function [data,mask,contour] = getMembrane(obj,membrane)
            idx = contains({obj.raw.movInfo.fileName},membrane,'IgnoreCase',true);
            
            if sum(idx)~=1
                warning(['no ' membrane ' file was found, operation aborted']);
                data = [];
                mask = [];
                contour = [];
                
            else
                data = obj.getFrame(1,membrane);
                
                switch lower(membrane)
                    case 'lamina'
                        contour = obj.Lamina.contour;
                        mask    = obj.Lamina.mask;
                    case 'nup'
                        contour = obj.NUP.contour;
                        mask    = obj.NUP.mask;
                    case 'lipid'
                        contour = obj.lipid.contour;
                        mask    = obj.lipid.mask;
                    otherwise
                        error('Unknown membrane requested')
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
        
        function makeMovie(obj)
            
            filename = [obj.raw.movInfo.Path filesep 'VisMov.gif'];
            frameRate = 5;
            trailing = 5;
            sizeParticles = 1000;
            frames = obj.raw.movInfo.nFrame;
            file2Test = [obj.raw.movInfo.path2Cell filesep 'mask.mat'];
            assert(isfile(file2Test),'no mask was found, please run getCellMask first')
            %load mask data
            allData = load(file2Test);
            field = fieldnames(allData);
            allData = allData.(field{1});
            
            %get traces and clean
            traces = obj.traces3D;
            allHeight = cellfun(@height,traces);
            traces(allHeight<10) = [];
            
            Fig = figure(1);
            hold on
            xLimit =[1 obj.raw.movInfo.Width*obj.raw.movInfo.pxSize];
            yLimit =[1 obj.raw.movInfo.Length*obj.raw.movInfo.pxSize];
            zLimit =[0 length(obj.raw.movInfo.planePos)*obj.info.dZ];
            
            [x,y,z] = sphere(32);
            x = x*sizeParticles/2;
            y = y*sizeParticles/2;
            z = z*sizeParticles/2;
            
            for i = 1:frames
                data2Render = allData(:,:,:,i);
                smoothISurface = isosurface(data2Render,1/2);
                
                % smoothing using compiled c code
                % smoothISurface = rendering3D.smoothpatch(iSurface,0,1);
                %comnvert to px
                smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1))*obj.raw.movInfo.pxSize;
                smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2))*obj.raw.movInfo.pxSize;
                smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3))*obj.info.dZ;
                
                p2 = patch(smoothISurface);
                p2.FaceColor = [0.4 0.4 0.4];
                p2.EdgeColor = 'none';
                
                title(['Frame ' num2str(i)]);
                axis image
                xlim(xLimit);
                ylim(yLimit);
                zlim(zLimit)
                xlim manual;
                ylim manual;
                zlim manual;
                
                
                for j = 1:length(traces)
                    currTrace = traces{j};
                    idx2Frame = currTrace.t==i;
                    idx = i-trailing:i;
                    idx = ismember(currTrace.t,idx);
                    if ~all(idx==0)
                        
                        data2Plot = currTrace(idx,:);
                        hold on
                        plot3(data2Plot.col,data2Plot.row,data2Plot.z,'color',[0 0 1])
                        try
                            X = x+currTrace.col(idx2Frame);
                            Y = y+currTrace.row(idx2Frame);
                            Z = z+currTrace.z(idx2Frame);
                        catch
                            X = x+data2Plot.col(end,:);
                            Y = y+data2Plot.row(end,:);
                            Z = z+data2Plot.z(end,:);
                        end
                        
                        surf(X,Y,Z,'LineStyle','none','Facecolor',[0.78 0.2 0.2])
                        
                        
                    end
                end
                camlight
                view(3)
                lighting gouraud
                axis image
                xlim(xLimit);
                ylim(yLimit);
                zlim(zLimit)
                xlim manual;
                ylim manual;
                zlim manual;
                
                drawnow;
                frame = getframe(Fig);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                if i == 1
                    
                    imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);
                    
                else
                    
                    imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');
                    
                end
                clf;
            end
            
            
            
            
        end
        
        function plotTraces(obj)
            frames = obj.raw.movInfo.nFrame;
            
            file2Test = [obj.raw.movInfo.path2Cell filesep 'mask.mat'];
            assert(isfile(file2Test),'no mask was found, please run getCellMask first')
            %load mask data
            allData = load(file2Test);
            field = fieldnames(allData);
            allData = allData.(field{1});
            
            
            %average the mask over time
            maskData = squeeze(mean(allData,4));
            %get traces  and clean up
            traces = obj.traces3D;
            allHeight = cellfun(@height,traces);
            traces(allHeight<10) = [];
            
            Fig = figure;
            hold on
            data2Render = maskData;
            smoothISurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
            % smoothISurface = rendering3D.smoothpatch(iSurface,0,1);
            %convert to px
            smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1))*obj.raw.movInfo.pxSize;
            smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2))*obj.raw.movInfo.pxSize;
            smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3))*obj.info.dZ;
            
            p2 = patch(smoothISurface);
            p2.FaceColor = [0.4 0.4 0.4];
            p2.EdgeColor = 'none';
            
            title(['Tracked traces']);
            axis image
            box on
            
            for j = 1:length(traces)
                currTrace = traces{j};
                plot3(currTrace.col,currTrace.row,currTrace.z)
                
            end
            view(2)
            
            
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
            
            factor_distance = 1;
            
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