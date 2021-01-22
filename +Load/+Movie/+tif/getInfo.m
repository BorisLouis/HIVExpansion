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

movieInfo.nFrames = nImages;
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

