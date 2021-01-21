function [ movieInfo ] = getInfo( path2File )
%GETINFO Summary of this function goes here
%   Detailed explanation goes here

[path,file,ext] = fileparts(path2File);

fID = fopen(path2File);
%read oif main file
tLine = fgetl(fID);
rawMetaData = {};
while ischar(tLine)
    tLine = fgetl(fID);
    rawMetaData = [rawMetaData; {tLine}];
end
fclose(fID);

%clean file from numbers:
idx2Char0 = cellfun(@ischar,rawMetaData);
rawMetaData = rawMetaData(idx2Char0);
%remove char 0 from file
cleanMetaData = strrep(rawMetaData,char(0),'');

%kill empty cells
idx2Empty = cellfun(@isempty,cleanMetaData);
metaData  = cleanMetaData(~idx2Empty);

%Extract data from file
%1 width
width = metaData(contains(metaData,'ImageWidth'));
imWidth = str2double(width{1}(strfind(width{1},'=')+1:end));

%height
height = metaData(contains(metaData,'ImageHeight'));
imHeight = str2double(height{1}(strfind(height{1},'=')+1:end));

%channel
channel = metaData(contains(metaData,'Ch='));
nCh = str2double(channel{1}(strfind(channel{1},'=')+1:end));

%frames
frame = metaData(contains(metaData,'Frames='));
nFrame = str2double(frame{1}(strfind(frame{1},'=')+1:end));

%pxSize
pxSize= metaData(contains(metaData,'HeightConvertValue='));
px = str2double(pxSize{1}(strfind(pxSize{1},'=')+1:end))*1000;

%extract number of time point and zSlice from file name
endFile = metaData(contains(metaData,'PtyFileNameE='));
tmpIdx = strfind(endFile{1},'\');
endFileN = endFile{1}(tmpIdx:end);

idxT = strfind(endFileN,'T');
idxZ = strfind(endFileN,'Z');

zSlice = str2double(endFileN(idxZ+1:idxT-1));
time   = str2double(endFileN(idxT+1:idxT+4));

movieInfo.Width  = imWidth;
movieInfo.Length = imHeight;
movieInfo.nFrame = time;
movieInfo.nPlane = zSlice;
movieInfo.nChan  = nCh;
movieInfo.totFrame = nFrame;
movieInfo.pxSize   = px;
movieInfo.Path   = path;

end

