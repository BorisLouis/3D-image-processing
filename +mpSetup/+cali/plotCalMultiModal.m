function plotCalMultiModal(ch1,ch2, ch3, ch4)
    
    figure()
    for i = 1:4
        channel1 = mean(ch1(:,:,i,:), 4);
        channel9 = mean(ch3(:,:,i,:), 4);
        Max = max(channel1(:));

        % tform = transformations{i,1};
        % channel9 = imwarp(channel9, tform, "OutputView", imref2d(size(channel1)));
        channel9 = channel9./((max(channel9(:)))*Max);
        subplot(2,8,i)
        imagesc(channel1)
        hold on
        subplot(2,8,i+8)
        imagesc(channel9)
        hold on
    end

    for i = 1:4
        channel1 = mean(ch2(:,:,i,:), 4);
        channel9 = mean(ch4(:,:,i,:), 4);
        Max = max(channel1(:));

        % tform = transformations{i+4,1};
        % channel9 = imwarp(channel9, tform, "OutputView", imref2d(size(channel1)));
        channel9 = channel9./((max(channel9(:)))*Max);
        subplot(2,8,i+4)
        imagesc(channel1)
        hold on
        subplot(2,8,i+8+4)
        imagesc(channel9)
        hold on
    end

    sgtitle("All frames' ROI's - corrected")
end