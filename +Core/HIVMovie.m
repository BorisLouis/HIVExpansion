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
            
            if iscell(raw)
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
                   movInfo(i).nFrames = movI.nFrames;
                   movInfo(i).nPlanes = movI.nPlanes;
                   movInfo(i).zSpacing = movI.zSpacing;
                   movInfo(i).zPosition = [0:movI.nFrames-1]*movI.zSpacing;
                    
                end
                warning('on')

                obj.raw.ext = raw{2};
                obj.raw.movInfo = movInfo;
                obj.raw.path = raw{1};
            else
                assert(isstruct(raw),'raw info is supposed to be a struct');
                obj.raw =raw;
            end

        end
        
        function getExtraInfo(obj)
            movInfo = obj.raw.movInfo;
            datatype = obj.dataType;
            for i = 1:length(movInfo)
                name = movInfo(i).fileName;                
                
                for j = 1:length(datatype)
                    chk = contains(name,datatype{j},'IgnoreCase',true);
                    if chk
                        obj.raw.movInfo(i).datatype = datatype{j};
                        break;
                                                        
                    end
                end
            end
            % here we check that there is only one file per datatype
            dataT = {obj.raw.movInfo.datatype};
            
            assert(isequal(sort(dataT), unique(dataT)),['Found more than one '...
                'file per datatype, please check directory and filenames, '...
                'Expect only one or zero of each of these files: NUP,Lamina,HIV,Lipid'])
                        
            disp(['We have detected ' num2str(length(movInfo)) ' files, please fill in the additional information']);
            
            prompt = {'FWHM(px)','pxSize(nm)'}';
            dlgTitle = 'Information about experimental parameters';
            numLines = 1;
            defaultVal = {'4','100'};
            answer = inputdlg(prompt,dlgTitle,numLines,defaultVal);
            
            assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')
            
            FWHM_px = str2double(answer(1));
            assert(~isnan(FWHM_px),'FWHM should be numerical');
            
            pxSize = str2double(answer(1));
            assert(~isnan(pxSize),'pxSize should be numerical');
            
            obj.info.pxSize = pxSize;
            obj.info.FWHM_px = FWHM_px;
            
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
end