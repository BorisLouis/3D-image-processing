function [ROI, ROIFull] = refineROItransf(ROI,ROIFull, shifts, sizes, idx1, idx2)

    for i = 1:size(shifts, 1)
        Transformation = shifts{i, 1};

        ROI(i,1) = round(ROI(i,1) - Transformation.Translation(1));
        ROI(i,2) = round(ROI(i,2) - Transformation.Translation(2));
        ROI(i,3) = round(ROI(i,3)./Transformation.Scale);
        ROI(i,4) = round(ROI(i,4)./Transformation.Scale);

        if ROI(i, 1) < 0
            ROI(i,1) = ROI(i,1) + ceil(abs(ROI(i, 1)));
        elseif ROI(i, 1) + ROI(i, 3) > sizes(2)
            ROI(i,1) = sizes(2) - ROI(i, 3);
        end

        if ROI(i, 2) < 0
            ROI(i, 2) = 0;
        elseif ROI(i, 2) + ROI(i, 4) > sizes(1)*2
            ROI(i,1) = sizes(1)*2 - ROI(i, 4);
        end

        ROIFull(i,1) = round(ROIFull(i,1) - Transformation.Translation(1));
        ROIFull(i,2) = round(ROIFull(i,2) - Transformation.Translation(2));
        ROIFull(i,3) = round(ROIFull(i,3)./Transformation.Scale);
        ROIFull(i,4) = round(ROIFull(i,4)./Transformation.Scale);

        if ROIFull(i, 1) < 0
            ROIFull(i,1) = ROIFull(i,1) + ceil(abs(ROIFull(i, 1)));
        elseif ROIFull(i, 1) + ROIFull(i, 3) > sizes(2)
            ROIFull(i,1) = sizes(2) - ROIFull(i, 3);
        end

        if ROIFull(i, 2) < 0
            ROIFull(i, 2) = 0;
        elseif ROIFull(i, 2) + ROIFull(i, 4) > sizes(1)*2
            ROIFull(i,1) = sizes(1)*2 - ROIFull(i, 4);
        end
    
    end
end

