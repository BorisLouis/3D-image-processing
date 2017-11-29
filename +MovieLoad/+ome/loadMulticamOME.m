function [ movC1, movC2, idx ] = loadMulticamOME( frameInfo, movieInfo, frames )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%TODO: redefine the number of frame so it changes if one or two camera are
%used, change the corresponding assert

assert(min(size(frames))==1,'Frames must be a vector')
assert(isstruct(frameInfo),'info must be a structure')
assert(max(frames)<= length(frameInfo) / 2, 'your request exceeds the total number of frames')
ImL = movieInfo.Length;
ImW = movieInfo.Width;
path2omes = movieInfo.Path;


nFramres = length(frames);
Cinfo    = {frameInfo.C};
Finfo    = {frameInfo.T};

cam =str2double(unique(Cinfo));% Extract indexes of the camera, the size of
%cam matches the number of camera used.

switch size(cam,2)
    case 1
        idx = zeros(nFramres,1);
    case 2 
        idx = zeros(nFramres,2);
    otherwise
        error('Unexpected number of cam loaded');
end
        
for iCam = 1:size(cam,2)
    C = num2str(cam(iCam));
    for fi = 1:nFramres
        frame = frames(fi);
        F = num2str(frame-1);
        tmp = and(strcmp(Cinfo,{C}),strcmp(Finfo,{F}));
%       test = find(tmp);
%             if length(test)>1
%                 disp('HAAAAA')
%             elseif isempty(test)
%                 disp('NOOOO')
%             end
        idx(fi,cam+1) = find(tmp);
    end
end

movC1 = uint16(zeros(ImL,ImW,nFramres));
movC2 = uint16(zeros(ImL,ImW,nFramres));

h = waitbar(0,'Please wait...');
steps = nFramres*2;

oldPath = '';
warning('off','all')

switch size(idx,2)
    
    case 1 %Only one Cam
        for fi = 1:nFramres
            i      = idx(fi,1);
            f2load = frameInfo(i).File;
            p2file = [path2omes filesep f2load];

            if ~strcmp(oldPath,p2file)

                if exist( 'tObj', 'var' )
                    tObj.close
                end
                tObj   = Tiff(p2file,'r');
                oldPath = p2file;

            end

            TifDir = str2double( frameInfo(i).IFD ) + 1;
        %     disp(TifDir)
            tObj.setDirectory(TifDir)
            movC1(:,:,fi) = tObj.read;
            waitbar(fi / steps)
        end
        
        movC2 = movC1; %We give the same movie for both camera (or leave
        %empty?)
        
    case 2 %when multiCam is used
        for fi = 1:nFramres
            i      = idx(fi,1);
            f2load = frameInfo(i).File;
            p2file = [path2omes filesep f2load];

            if ~strcmp(oldPath,p2file)

                if exist( 'tObj', 'var' )
                    tObj.close
                end
                tObj   = Tiff(p2file,'r');
                oldPath = p2file;

            end

            TifDir = str2double( frameInfo(i).IFD ) + 1;
        %     disp(TifDir)
            tObj.setDirectory(TifDir)
            movC1(:,:,fi) = tObj.read;
            waitbar(fi / steps)
        end
        
        for fi = 1:nFramres
            i      = idx(fi,2);
            f2load = frameInfo(i).File;
            p2file = [path2omes filesep f2load];

            if ~strcmp(oldPath,p2file)

                if exist( 'tObj', 'var' )
                    tObj.close
                end
                tObj   = Tiff(p2file,'r');
                oldPath = p2file;

            end

            TifDir = str2double( frameInfo(i).IFD ) + 1;
        %     disp(TifDir)
            tObj.setDirectory(TifDir)
            movC2(:,:,fi) = tObj.read; 
            waitbar((fi + nFramres) / steps)
        end
        
        
    otherwise
        error('Idx has an unexpected size');
        
end

tObj.close
close(h) 
warning('on','all')

end

