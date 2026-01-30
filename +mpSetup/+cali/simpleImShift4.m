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


        figure()
        subplot(1,2,1)
        imshowpair(imCh1,imChi);
        title("before correction");
        hold on 

        config = "monomodal";
        transf = "rigid";
        [optimizer,metric] = imregconfig(config);
        tform = imregcorr(imChi,imCh1,transf);
        % if any(or(tform.Translation > 5, tform.Translation < -5))
        %     transf = "rigid";
        %     tform2 = imregcorr(imChi,imCh1,transf);
        %     if any(or(tform.Translation > 5, tform.Translation < -5))
        %         transf = "translation";
        %         tform3 = imregcorr(imChi,imCh1,transf);
        %         tform.Translation = tform3.Translation;
        %         if any(or(tform.Translation > 5, tform.Translation < -5))
        %             config = "multimodal";
        %             [optimizer,metric] = imregconfig(config);
        %             tform4 = imregcorr(imChi,imCh1,transf);
        %             tform.Translation = tform4.Translation;
        %             if any(or(tform.Translation > 5, tform.Translation < -5))
        %                 tform4 = tform;
        %                 tform.Translation = [0 0];
        %             end
        %         end
        %     end
        % end
        tformChanged = tform;

        tformChanged.RotationAngle = 0;
        tformChanged.R = [1, 0; 0, 1];
        tformChanged.A = [tformChanged.Scale, 0, tformChanged.Translation(1); 0, tformChanged.Scale, tformChanged.Translation(2); 0, 0, 1];

        movingRegistered = imwarp(imChi,tformChanged, "OutputView",imref2d(size(imCh1)));

        subplot(1,2,2)
        imshowpair(imCh1,movingRegistered);
        title("after correction");
        hold on 

        imShifts{chIdx,2} = {multissim(movingRegistered, imCh1)};
        imShifts{chIdx,1} = tform;        
        imShifts{chIdx,3} = transf;
    end  
end

