function [ transformations ] = GetMagnificationScale(cam1, cam2, cam3, cam4, inFocus1, inFocus2)
    SimilarityScore = [];
    
    for chIdx = 1:length(inFocus1)
        Plane9infocus = [];
        Plane1infocus = [];

        idxPlane1 = [inFocus1.globalch]==chIdx;
        idxPlane9 = [inFocus2.globalch]==chIdx;

        focus1     = inFocus1(idxPlane1).frame;
        focus9     = inFocus2(idxPlane9).frame;
    
        camIdx1 = inFocus1(idxPlane1).cam;
        camIdx9 = inFocus2(idxPlane9).cam;

        if camIdx1 == 1
            Plane1infocus = double(cam1(:,:,inFocus1(idxPlane1).ch,focus1));
        else
            Plane1infocus = double(cam2(:,:,inFocus1(idxPlane1).ch,focus1));
        end

        if camIdx9 == 3
            Plane9infocus = double(cam3(:,:,inFocus2(idxPlane9).ch,focus9));
        else
            Plane9infocus = double(cam4(:,:,inFocus2(idxPlane9).ch,focus9));
        end
    
        %Plane9infocus = Plane9infocus./max(Plane9infocus(:));
        %Plane1infocus = Plane1infocus./max(Plane1infocus(:));

        figure()
        subplot(1,2,1)
        imshowpair(Plane1infocus,Plane9infocus)
        title("before correction");
        hold on
        
        config = "multimodal";
        transf = "similarity"; %similarity

        [optimizer,metric] = imregconfig(config);
        tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
        transformations{chIdx, 1} = tform;
        movingRegistered = imwarp(Plane9infocus,tform,"OutputView",imref2d(size(Plane1infocus)));
        SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);
        transformations{chIdx, 2} = SimilarityScore(chIdx, 1);
        transformations{chIdx, 3} = config;
        transformations{chIdx, 4} = transf;

        % if SimilarityScore(chIdx, 1) < 0.10
        %     transf = "translation";
        % 
        %     [optimizer,metric] = imregconfig(config);
        %     tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
        %     movingRegistered = imwarp(Plane9infocus,tform,"OutputView",imref2d(size(Plane1infocus)));
        %     Score = multissim(movingRegistered, Plane1infocus)
        %     if Score > SimilarityScore(chIdx, 1)
        %         transformations{chIdx, 1} = tformstruct;
        %         SimilarityScore(chIdx, 1) = Score;
        %         transformations{chIdx, 2} = SimilarityScore(chIdx, 1);
        %         transformations{chIdx, 3} = config;
        %         transformations{chIdx, 4} = transf;
        %     else
        %     end
        % else
        % end

   
        subplot(1,2,2)
        imshowpair(Plane1infocus,movingRegistered);
        title("after correction");
        sgtitle(append("Plane ", num2str(chIdx), " x Plane ", num2str(chIdx+8)));
    end
end

