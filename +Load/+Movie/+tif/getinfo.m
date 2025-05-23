function [ movieInfo ] = getInfo( path2file )
%GETINFO Summary of this function goes here
%   Detailed explanation goes here

warning('off','all')
tObj = Tiff(path2file,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
movieInfo.Path   = fileparts(path2file);

tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end;
    tObj.nextDirectory();
end
tObj.setDirectory(1)
warning('on','all')

movieInfo.Frame_n = tfl; 
tObj.close

end

