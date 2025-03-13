function [ imShifts ] = simpleImShift( inFocus, cam1, cam2)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    nPlanes = length(inFocus);
    % imShifts = zeros(length(inFocus), 3);
    imShifts = {};
    focus = inFocus(1).frame;
    
    imCh1 = cam1(:,:,1,focus);
    se = strel('disk', 12);
    bgCh1 = imopen(imCh1, se);
    imCh1 = imCh1 - bgCh1;

    for chIdx = 2:nPlanes
    
        focus     = inFocus(chIdx).frame;
        if chIdx<nPlanes/2+1
            % I look at first cam
            imChi = cam1(:,:,chIdx,focus);
        else
            % I look at second cam
            imChi = cam2(:,:,chIdx-nPlanes/2,focus);
        end
        bgChi = imopen(imChi, se);
        imChi = imChi - bgChi;
        
        % 
        % figure()
        % subplot(1,2,1)
        % imshowpair(imCh1,imChi);
        % title("before correction");
        % hold on 

        config = "monomodal";
        transf = "similarity";
        [optimizer,metric] = imregconfig(config);
        tform = imregcorr(imChi,imCh1,transf);
        tformChanged = tform;
        tformChanged.RotationAngle = 0;
        tformChanged.R = [1, 0; 0, 1];
        tformChanged.A = [tformChanged.Scale, 0, tformChanged.Translation(1); 0, tformChanged.Scale, tformChanged.Translation(2); 0, 0, 1];
        movingRegistered = imwarp(imChi,tformChanged, "OutputView",imref2d(size(imCh1)));

        % subplot(1,2,2)
        % imshowpair(imCh1,movingRegistered);
        % title("after correction");
        % sgtitle(append("Plane ", num2str(chIdx), " x Plane ", num2str(chIdx+8)));

        % storing the shifts
        % imShifts(chIdx,1:2) = [tform.Translation(1), tform.Translation(2)];
        % imShifts(chIdx,3) = tform.Scale;
        % imShifts(chIdx,4) = tform.RotationAngle;   
        imShifts{chIdx,2} = {multissim(imChi, imCh1)};
        imShifts{chIdx,1} = tform;        
    end
    
end

