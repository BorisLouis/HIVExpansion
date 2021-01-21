function [ movieInfo ] = getInfo( path2File )
%GETINFO Summary of this function goes here
%   Detailed explanation goes here
tObj = Tiff(path2File,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
movieInfo.Path   = fileparts(path2File);
tObj.close

end

