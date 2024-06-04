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
        
        config = "monomodal";
        transf = "similarity";

        [optimizer,metric] = imregconfig(config);
        tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
        transformations{chIdx, 1} = tform;
        movingRegistered = imwarp(Plane9infocus,transformations{chIdx, 1},"OutputView",imref2d(size(Plane1infocus)));
        SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);

        if SimilarityScore(chIdx, 1) < 0.2
            config = "multimodal";
    
            [optimizer,metric] = imregconfig(config);
            tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
            transformations{chIdx, 1} = tform;
            movingRegistered = imwarp(Plane9infocus,transformations{chIdx, 1},"OutputView",imref2d(size(Plane1infocus)));
            SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);

            if SimilarityScore(chIdx, 1) < 0.2
                config = "monomodal";
                transf = "affine";
        
                [optimizer,metric] = imregconfig(config);
                tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
                transformations{chIdx, 1} = tform;
                movingRegistered = imwarp(Plane9infocus,transformations{chIdx, 1},"OutputView",imref2d(size(Plane1infocus)));
                SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);

                if SimilarityScore(chIdx, 1) < 0.2
                    config = "multimodal";
            
                    [optimizer,metric] = imregconfig(config);
                    tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
                    transformations{chIdx, 1} = tform;
                    movingRegistered = imwarp(Plane9infocus,transformations{chIdx, 1},"OutputView",imref2d(size(Plane1infocus)));
                    SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);
                else 
                end
            else 
            end
        else
        end

   
        subplot(1,2,2)
        imshowpair(Plane1infocus,movingRegistered);
        title("after correction");
        sgtitle(append("Plane ", num2str(chIdx), " x Plane ", num2str(chIdx+8)));
    end
end

