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
        
        function [gBW] = segmentLamina(obj)
            
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
                se = strel('disk',5);
                gBW = imclose(gBW,se);
                gBW = bwareaopen(gBW,5000);
                gBW = imfill(gBW,'holes');            
                disp('Extracting Contour')

                %here we obtain the cell contour
                contour = obj.getLaminaContour(gBW);

                obj.Lamina.mask = gBW;
                obj.Lamina.contour = contour;
%                 
%                 %% building 3D mask
%                 mask3D = zeros(size(IMs));
%                 for j = 1:size(IMs,3)
%                     mask3D(:,:,j) = Misc.mpoly2mask(contourtmp{1},[size(IMs,1) size(IMs,2)],'style','ij');                    
% %                     if ~isempty(fContour)
% %                         idx = find(fContour(:,3)==j);
% %                         mask3D(:,:,j) = poly2mask(fContour(idx,2),fContour(idx,1),size(IMs,1),size(IMs,2));
% %                     end
%                 end
% 
%                 %% Cleaning mask
% 
%                 test = bwlabeln(mask3D);
%                 data = regionprops3(test,'Volume','VoxelList','VoxelIdxList');
% 
%                 for j = 1:height(data)
%                    currentIdx = data.VoxelList{j,:};
%                    if length(unique(currentIdx(:,3)))>1
% 
%                    else
%                        idx2Delete = data.VoxelIdxList{j,:};
%                        mask3D(idx2Delete) = 0;
%                    end
% 
%                 end
%              

%                 allMask(:,:,:,i) = logical(mask3D);

    %             
    %             filename = [obj.raw.movInfo.path2Cell filesep 'mask.mat'];
    %             save(filename,'allMask','-v7.3');
            end
        end
        
        
        function [gBW] = segmentNUP(obj)
            
            idx = contains({obj.raw.movInfo.fileName},'NUP','IgnoreCase',true);
            if sum(idx)~=1
                warning('no lamina file was found, operation aborted');
                
            else
                data = obj.getFrame(1,'NUP');
                % Gaussian filtering
                S = 2; 
                % size of pixel in z vs x/y
                zFactor = obj.raw.movInfo(1).zSpacing/obj.info.pxSize;
                sigma = [S,S,S/zFactor];
                IMs = imgaussfilt3(data, sigma);
                disp('DONE with filtering ------------')

              %  gBW = imbinarize(uint16(IMs),'adaptive','Sensitivity',threshold);
                gBW = imbinarize(uint16(IMs));
                se = strel('disk',5);
                gBW = imclose(gBW,se);
                %gBW = bwareaopen(gBW,5000);
                gBW = imfill(gBW,'holes');            
                disp('Extracting Contour')

                %here we obtain the cell contour
                contour = obj.getLaminaContour(gBW);

                obj.Lamina.mask = gBW;
                obj.Lamina.contour = contour;
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
        
         function [contour] = getLaminaContour(~,gBW)
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
    end
end