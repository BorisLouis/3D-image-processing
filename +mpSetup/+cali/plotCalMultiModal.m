function plotCalMultiModal(ch1,ch2, ch3, ch4, inFocus1, inFocus2)
    
    figure()
    for i = 1:4
        channel1 = ch1(:,:,i,inFocus1(i).frame);
        channel9 = ch3(:,:,i,inFocus2(i).frame);
        channel9 = imresize(channel9, size(channel1));
      
        subplot(2,4,i)
        imshowpair(channel1, channel9);
        % subplot(2,8,i)
        % imagesc(channel1)
        % title(append('Plane ', num2str(i)))
        % hold on
        % subplot(2,8,i+8)
        % imagesc(channel9)
        % title(append('Plane ', num2str(i)))
        % hold on
    end

    for i = 1:4
        channel2 = ch2(:,:,i,inFocus1(i+4).frame);
        channel4 = ch4(:,:,i,inFocus2(i+4).frame);
        channel4 = imresize(channel4, size(channel2));
        % Max = max(channel1(:));

        % tform = transformations{i+4,1};
        % channel9 = imwarp(channel9, tform, "OutputView", imref2d(size(channel1)));
        % channel9 = channel9./((max(channel9(:)))*Max);
        subplot(2,4,i+4)
        imshowpair(channel2, channel4);
        % subplot(2,8,i+4)
        % imagesc(channel2)
        % title(append('Plane ', num2str(i+4)))
        % hold on
        % subplot(2,8,i+8+4)
        % imagesc(channel4)
        % title(append('Plane ', num2str(i+4)))
        % hold on
    end

    sgtitle("All frames' ROI's - corrected")
end