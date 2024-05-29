function [ ROI1, ROI2, ROI2FullCam ] = refineROI( ROI1, ROI2, ROI2FullCam, im_shifts)
%REFINEROI improves the ROI to have the smallest shifts possibles
%   Detailed explanation goes here
    
    row_w = ROI1(1,4);
    col_w = ROI1(1,3);
    %get min and max imShift for row
    mars = max(im_shifts(:,1));
    mirs = min(im_shifts(:,1));
    %get in and max imShift for Col
    macs = max(im_shifts(:,2));
    mics = min(im_shifts(:,2));
    %calculate the width 
    row_w = row_w - mars + mirs;
    col_w = col_w - macs + mics;
    
    row_i1 = ROI1(:,2)+(mars-(im_shifts(:,1)));
    row_i2 = ROI2(:,2)+(mars-(im_shifts(:,1)));
    row_i2FullCam = ROI2FullCam(:,2)+(mars-(im_shifts(:,1)));

    col_i = ROI(:,1)+(macs-(im_shifts(:,2)));
    col_i = ROI(:,1)+(macs-(im_shifts(:,2)));
    col_i = ROI(:,1)+(macs-(im_shifts(:,2)));
    col_i = ROI(:,1)+(macs-(im_shifts(:,2)));
    
    ROI(:,4) = row_w;
    ROI(:,2) = row_i;
    ROI(:,3) = col_w;
    ROI(:,1) = col_i;
end

