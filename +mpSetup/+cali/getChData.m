function [ chData1c, chData2c ] = getChData( data1c, data2c, ROI )
%GETCHDATA get channel data from cam data and ROI, for historical reasons
%it can work on a single cam data, Im not sure if I want to keep this later
%on
 
    if and(~isempty(data1c),~isempty(data2c))
        im_size1 = size(data1c);
        im_size2 = size(data2c);
        assert(all(im_size1 == im_size2), 'camera data are not compatible')
    end
    
    if isempty(data1c)
        chData1c = [];
        warning ('No cam 1 data, weird, returning empty')
    else
        im_size1 = size(data1c);
        c_win = ROI(1,3);
        r_win = ROI(1,4);
        
        if length(im_size1)==3
            N = im_size1(3);
        elseif length(im_size1)==2
            N = 1;
        else
            error('Unexpected problem')
        end
        nChan = 1:size(ROI,1);
        chanIdx = 1:length(nChan)/2;
        if (r_win>=im_size1(1))
            warning('y-size of calibration ROI is different from the data to calibrate');
            disp('re-adapting the ROI size')

            %extract the shift in y
            yShift = ROI(:,2);
            %shift it back related to 0 (minimum shift should be 0)
            yShift = yShift -min(yShift)+1;
            %largest usable ROI is the size of the image - the largest y
            %shift, then we can use the same ROI for everyone.
            maxROI = im_size1(1) - max(yShift);
            
            ROI(chanIdx,2) = yShift(chanIdx);
            ROI(chanIdx,4) = maxROI;
            
            
        end
        
        chData1c = uint16(zeros(ROI(1,4),ROI(1,3),size(ROI,1)/2,N));
        for i=chanIdx
            col1 = ROI(i,1); 
            col2 = ROI(i,1) + ROI(i,3)-1;
            row1 = ROI(i,2);
            row2 = ROI(i,2) + ROI(i,4)-1;
            chData1c(:,:,i,:) = data1c(row1:row2,col1:col2,:);
        end
       
    end
    
    if isempty(data2c)
        chData2c = [];
    else
        im_size2     = size(data2c);
        nChan = 1:size(ROI,1);
        chanIdx = nChan(length(nChan)/2)+1:nChan(end);
        r_win = ROI(chanIdx(1),4);
        
        if length(im_size2)==3
            N = im_size2(3);
        elseif length(im_size2)==2
            N = 1;
        else
            error('Unexpected problem')
        end
        
       
        if (r_win>=im_size2(1))
            warning('y-size of calibration ROI is different from the data to calibrate');
            disp('re-adapting the ROI size')
            %extract the shift in y
%             yShift = ROI(chanIdx,2);
%             %shift it back related to 0 (minimum shift should be 0)
%             yShift = yShift -min(yShift)+1;
%             %largest usable ROI is the size of the image - the largest y
%             %shift, then we can use the same ROI for everyone.
%             maxROI = im_size1(1) - max(yShift);
            
            ROI(chanIdx,2) = yShift(chanIdx);
            ROI(chanIdx,4) = maxROI;
            
            
        end
        chData2c = uint16(zeros(ROI(chanIdx(1),4),ROI(chanIdx(1),3),size(ROI,1)/2,N));
        for i=nChan(length(nChan)/2)+1:nChan(end)
            col1 = ROI(i,1); 
            col2 = ROI(i,1) + ROI(i,3) - 1;
            row1 = ROI(i,2);
            row2 = ROI(i,2) + ROI(i,4) - 1;
            chData2c(:,:,i-size(ROI,1)/2,:) = data2c(row1:row2,col1:col2,:);
        end
    end
    

end

