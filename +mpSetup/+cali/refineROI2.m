function [ new_ROI ] = refineROI2( ROI, im_shifts, multiModal )
%REFINEROI improves the ROI to have the smallest shifts possibles
%   Detailed explanation goes here
    new_ROI(:,1) = round(ROI(:,1) - im_shifts(:,1));
    new_ROI(:,2) = round(ROI(:,2) - im_shifts(:,2));
    new_ROI(:,3) = round(ROI(:,3));
    new_ROI(:,4) = round(ROI(:,4));

    for i = 2:size(new_ROI,1)
        x = new_ROI(i,1);
        y = new_ROI(i,2);
        width = new_ROI(i,3);
        height = new_ROI(i,4);

        x_center = x + width / 2;
        y_center = y + height / 2;

        new_width = width * im_shifts(i,3);
        new_height = height * im_shifts(i,3);

        new_x = x_center - new_width / 2;
        new_y = y_center - new_height / 2;

        new_ROI(i,:) = round([new_x, new_y, new_width, new_height]);
    end
end
