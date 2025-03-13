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
        wScaled = w / (tform.Scale);
        hScaled = h / (tform.Scale);
        xScaled = x;
        yScaled = y;

        %% Now do the translation
        xTransl = xScaled - (tform.Translation(1) / tform.Scale);
        yTransl = yScaled - (tform.Translation(2) / tform.Scale);

        new_ROI(i,:) = [xTransl, yTransl, wScaled, hScaled];
    end
end
