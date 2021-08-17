function [ movieInfo ] = getInfo( path2File )
%GETINFO Summary of this function goes here
%   Detailed explanation goes here
tObj = Tiff(path2File,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
movieInfo.Path   = fileparts(path2File);

desc = tObj.getTag(270);

idx1 = strfind(desc,'images=');
idx2 = strfind(desc,'slices=');

nImages = str2double(desc(idx1+7:idx2-1));
if isnan(nImages)
    idx3 = strfind(desc,'frames=');
    nImages = str2double(desc(idx2+7:idx3-1));
end
if isnan(nImages)
    nImages =1;
    warning('Could not determine the number of image in file, assuming 1');
%assert(~isnan(nImages),'Something went wrong when reading the file');
end

movieInfo.nPlanes = nImages;
movieInfo.nFrames = 1;
idx1 = strfind(desc,'spacing=');
idx2 = strfind(desc,'loop');

if isempty(idx1)
    zSpacing = 0;
else
    zSpacing = str2double(desc(idx1+8:idx2-1));
end

movieInfo.zSpacing = zSpacing;


tObj.close

end

