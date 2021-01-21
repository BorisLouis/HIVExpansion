classdef HIVMovie < handle
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        raw
        dataType = {'lamina','NUP','HIV','lipid'}
    end
    
    properties
        info
    end
    
    methods
        %Constructor
        function obj = HIVMovie(ext,info,varargin)
            %MOVIE Construct an instance of this class
            %   We allow the user to create this object with various number
            %   of input allowing therefore to restart the analysis at any
            %   steps in the process.
            
            if isempty(varargin{1})          
                [path] = uigetdir();
            else
                path = varargin{1};
            end
            
            obj.raw = {path; ext};
            obj.info = info;
            
        end
        
        function set.info(obj,inform)
            obj.info = inform;
     
        end
        
        function set.raw(obj,raw)
            
            assert(isfolder(raw{1}),'The path provided does not appear to exist');
            assert(ischar(raw{2}),'Extension needs to be a string');
            obj.checkExtension(raw{2});
            
            
            [files] = obj.getFileInPath(raw{1},raw{2});
            ext = strrep(raw{2},'.','');
            movInfo = struct();
            warning('off')
            
            for i = 1:length(files)
               [movI] = Load.Movie.(ext).getInfo([files(i).folder filesep files(i).name]); 
               
               movInfo(i).path = [files(i).folder filesep files(i).name];
               movInfo(i).fileName = [files(i).name];
               movInfo(i).Width = movI.Width;
               movInfo(i).Length = movI.Length;           
               
               
            end
            warning('on')
            
            obj.raw.ext = raw{2};
            obj.raw.movInfo = movInfo;          
            
            
        end
        
        function getExtraInfo(obj)
            movInfo = obj.raw.movInfo;
            datatype = obj.dataType;
            for i = 1:length(movInfo)
                name = movInfo(i).fileName;                
                
                for j = 1:length(datatype)
                    chk = contains(name,datatype{i},'IgnoreCase',true);
                    if chk
                        obj.raw.movInfo(i).datatype = datatype{i};
                        break;
                                                        
                    end
                end
            end
            
            disp(['We have detected ' num2str(length(movInfo)) ' files, please fill in the additional information']);
            
            prompt = {'FWHM(px)','pxSize(nm)'}';
            dlgTitle = 'Information about experimental parameters';
            numLines = 1;
            defaultVal = {'2','100'};
            answer = inputdlg(prompt,dlgTitle,numLines,defaultVal);
            
            assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')
            
            FWHM_pix = str2double(answer(1));
            assert(~isnan(FWHM_pix),'FWHM should be numerical');
            
            pxSize = str2double(answer(1));
            assert(~isnan(FWHM_pix),'FWHM should be numerical');
            
            obj.info.pxSize = pxSize;
            obj.info.FWHM_pix = FWHM_pix;
            
        end
        
        function giveInfo(obj)
            %Make a prompt asking some question to the user.
            prompt = {'FWHM (px):', 'z interval (nm)'};
            dlgTitle = 'Information about experimental parameters';
            numLines = 1;
            defaultVal = {'2','250'};
            answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
            
            assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')
                        
            FWHM_pix = str2double(answer(1));
            assert(~isnan(FWHM_pix),'FWHM should be numerical');
            
            dZ = str2double(answer(2));
            assert(~isnan(dZ),'Z interval should be numerical, expressed in micrometer');
            
            %Calculate some setup parameters
            %sigma_nm = 0.25 * emW/NA;
            sigmaPix = FWHM_pix/2.355;
            %store info
            obj.info.FWHM_px =  FWHM_pix;
            obj.info.sigma_px = sigmaPix;
            obj.info.dZ = dZ;
            
        end
        
        function loadData(obj,chInfo)
            
            ext = lower(obj.raw.ext);
            path = obj.raw.path;
            [file2Analyze] = Core.HIVMovie.getFileInPath(path,ext);
            nDir = [file2Analyze(1).folder,filesep 'calibrated'];
            if ~isfolder(nDir)
                mkdir(nDir);
            end
            if isempty(file2Analyze)
                warning('Did not find any file of the extension provided in the folder provided');
            end
            extNoDot = strrep(ext,'.','');            
            endFile = file2Analyze(end).name;
            movInfo = Load.Movie.(extNoDot).getInfo([file2Analyze(1).folder filesep endFile]);
            
            path2Files = [file2Analyze.folder filesep file2Analyze.name,'.files'];
            f2Read = dir(path2Files);
            idx2Tif = contains({f2Read.name},'.tif');
            f2Read = f2Read(idx2Tif);
            
            exFile = f2Read(1).name;
                        
            %getIdx
            idx1 = strfind(exFile,'T');
            idx2 = strfind(exFile,'Z');
            idx3 = strfind(exFile,'C');
            idx4 = strfind(exFile,'.tif');
            %check if previous data was found
            try
                test = Core.HIVMovie.getFileInPath([nDir filesep 'part'],'.tif');
            catch ME
                if strcmp(ME.message,'The path given is not a folder')
                    test=0;
                else
                    error(ME.message)
                end
            end
            %if we found data
            if length(test) == movInfo.nPlane
                disp('Previous data found, loading from existing file');
                fileName = [nDir filesep 'cell'];
                movInfo.path2Cell = fileName;
                fileName = [nDir filesep 'part'];
                movInfo.path2Part = fileName;

            else
            
                partData = zeros(movInfo.Length,movInfo.Width,movInfo.nPlane,movInfo.nFrame,'uint16');
                cellData = partData;

                h= waitbar(0,'Separating data channel');
                for i = 1 : movInfo.totFrame
                    currentFile = [f2Read(i).folder filesep f2Read(i).name];
                    currentName = f2Read(i).name;
                    currentIm   = Load.Movie.tif.getFrames(currentFile,1);

                    cFrame = str2double(currentName(idx1+1:idx4-1));
                    cPlane = str2double(currentName(idx2+1:idx1-1));
                    cChan  = str2double(currentName(idx3+1:idx2-1));

                    if strcmp(chInfo.(['ch0' num2str(cChan)]),'particles')
                        partData(:,:,cPlane,cFrame) = currentIm;
                    elseif strcmp(chInfo.(['ch0' num2str(cChan)]),'cells')
                        cellData(:,:,cPlane,cFrame) = currentIm;
                    elseif cChan>2
                        warning('More than two channels, not considering additional channels')
                    else
                        error('Something is wrong, please check data and channel information')
                    end

                    waitbar(i/movInfo.totFrame,h,'Separating data channel')

                end

                fileName = [nDir filesep 'cell'];
                obj.saveChannel(cellData,fileName)
                movInfo.path2Cell = fileName;

                fileName = [nDir filesep 'part'];
                obj.saveChannel(partData,fileName)
                movInfo.path2Part = fileName;
                close(h);
            end
            
            movInfo.dataDir   = nDir;
            dZ = obj.info.dZ;
            movInfo.planePos = 0:dZ:(movInfo.nPlane-1)*dZ;
            
            obj.raw.movInfo = movInfo;
                       
        end     
    end
    methods (Static)
        function [frames]       = checkFrame(frames,maxFrame)
            %Short method that make sure that the frame are making sense.
            testFrame = mod(frames,1);
            
            if all(testFrame<1e-4)
                
            else
                
                frames = round(frames);
                warning('Some frame were not integers, they were rounded');
                
            end
            
            assert(isvector(frames),'Frames should be a vector of integers');
            assert(max(frames) <= maxFrame(1),'Request exceeds max frame');
            assert(min(frames) >0, 'Indexing in matlab start from 1');
            
        end
        
        function [file2Analyze] = getFileInPath(path, ext)
            %Small method to extract the file of a certain extension in a
            %given path
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            assert(isfolder(path),'The path given is not a folder')
            
            folderContent = dir(path);
            index2Images  = contains(lower({folderContent.name}),ext,'IgnoreCase',true);
            file2Analyze  = folderContent(index2Images);
            idx2File      = [file2Analyze.isdir];
            file2Analyze(idx2File) = [];
            
        end
        
        function checkExtension(ext)
            ext = lower(ext);
            extensionList = {'.tif'};
            
            check = contains(extensionList,ext);
            
            assert(sum(check)==1,'Error, unknown extension provided');
            
        end
        
        function [ext] = getFileExt(path)
            
            if isfolder(path)
                ext = [];
            elseif isfile(path)
                
                [~,name,ext1] = fileparts(path);
                [~,~,ext2] = fileparts(name);
                
                ext = [lower(ext2),lower(ext1)];
            end
        end
    end
    
    methods (Access = private)
       function saveChannel(~,data,calDir)
            
            mkdir(calDir);         
                       
            for i = 1:size(data,3)
                data2Store = squeeze(data(:,:,i,:));
                fName = sprintf('calibratedPlane%d.tif',i);
              
                fPathTiff = [calDir filesep fName];
                
                t = Tiff(fPathTiff, 'a');
                
                if isa(data2Store,'uint32')
                    t = dataStorage.writeTiff(t,data2Store,32);
                else
                    t = dataStorage.writeTiff(t,data2Store,16);
                end
                
                t.close;
                
            end
            
        end 
    
    end
end