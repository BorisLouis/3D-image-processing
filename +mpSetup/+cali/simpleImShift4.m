function [ imShifts ] = simpleImShift( inFocus, cam1, cam2)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    nPlanes = length(inFocus);
    imShifts = zeros(length(inFocus), 3);
    focus = inFocus(1).frame;
    
    imCh1 = cam1(:,:,1,focus);
    
    for chIdx = 2:nPlanes
    
        focus     = inFocus(chIdx).frame;
        if chIdx<nPlanes/2+1
            % I look at first cam
            imChi = cam1(:,:,chIdx,focus);
        else
            % I look at second cam
            imChi = cam2(:,:,chIdx-nPlanes/2,focus);
        end
        

        figure()
        subplot(1,2,1)
        imshowpair(imCh1,imChi);
        title("before correction");
        hold on 

        config = "monomodal";
        transf = "similarity";
        [optimizer,metric] = imregconfig(config);
        tform = imregcorr(imChi,imCh1,transf);
        movingRegistered = imwarp(imChi,tform, "OutputView",imref2d(size(imCh1)));

        subplot(1,2,2)
        imshowpair(imCh1,movingRegistered);
        title("after correction");
        sgtitle(append("Plane ", num2str(chIdx), " x Plane ", num2str(chIdx+8)));

        % storing the shifts
        imShifts(chIdx,1:2) = [tform.Translation(1), tform.Translation(2)];
        imShifts(chIdx,3) = tform.Scale;
        imShifts(chIdx,4) = multissim(imChi, imCh1);
        
    end
    
end

