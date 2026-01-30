function [ new_ROI ] = refineROI2( ROI, im_shifts, multiModal )
%REFINEROI improves the ROI to have the smallest shifts possibles

%   Detailed explanation goes here

    new_ROI(1,:) = ROI(1,:);
    
    for i = 2:size(ROI,1)
        tform = im_shifts{i,1};
        ROIcurrent = ROI(i,:);
        x = ROIcurrent(1);
        y = ROIcurrent(2);
        w = ROIcurrent(3);
        h = ROIcurrent(4);

        %% first apply scaling
        % if strcmp(im_shifts{i,3}, "translation")
        %     wScaled = w;
        %     hScaled = h;
        % else
            wScaled = w / (tform.Scale);
            hScaled = h / (tform.Scale);
        % end
        xScaled = x;
        yScaled = y;

        %% Now do the translation
        % if strcmp(im_shifts{i,3}, "translation")
        %     xTransl = xScaled - (tform.Translation(1));
        %     yTransl = yScaled - (tform.Translation(2));
        % else
            xTransl = xScaled - (tform.Translation(1) / tform.Scale);
            yTransl = yScaled - (tform.Translation(2) / tform.Scale);
        % end

        new_ROI(i,:) = [xTransl, yTransl, wScaled, hScaled];
    end

    MinROI = min(new_ROI(:,2));
    if MinROI < 1
        new_ROI(:,2) = new_ROI(:,2) + ceil(abs(MinROI)); 
    end

    MaxROI = max(new_ROI(:,1) + new_ROI(:,3));
    if MaxROI > 2048
        new_ROI(:,1) = new_ROI(:,1) - ceil(MaxROI - 2048); 
    end
   

end
