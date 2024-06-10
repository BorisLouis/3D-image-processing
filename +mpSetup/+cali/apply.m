function [data,isTransmission] = apply( cam1, cam2, cal )
%APPLY rearranges data and corrects for intensity diffecences
%between channels
%   Detailed explanation goes here

if or(isempty(cam1),isempty(cam2))
    error('Calibration expect 2 cameras, update will be performed later');
else
    %test if some of the data is transmission
    S1 = double(median(median(median(cam1))));
    S2 = double(median(median(median(cam2))));
    
    if abs(S1-S2) > 3*min([S1,S2])
        warning('One cam is much brighter than the other, assuming Transmission data');
        nChan = size(cal.ROI1,1)/2;
        if S1< S2
            
            isTransmission1 = [zeros(nChan,1);ones(nChan,1)];
        
        else
            
            isTransmission1 = [ones(nChan,1);zeros(nChan,1)];
            
        end
    else 
        isTransmission1 = zeros(size(cal.ROI1,1),1);
        if cal.multiModal == true
            isTransmission2 = zeros(size(cal.ROI2,1),1);
        else
        end
    end
    
    h = waitbar(0,'Please wait applying calibration');
    % if need be we flip camera 2, this is generally the case
    if cal.flipCam2
        waitbar(.1,h,'Flipping cam')
        cam2 = flip(cam2,2);
    end
    waitbar(.2,h,'Gettingg channel data')
    [ chC1, chC2 ] = mpSetup.cali.getChData( cam1, cam2, cal.ROI1 );
    sTmp = size(chC1);
    sTmp(3) = sTmp(3)*2;
    data1 = ones(sTmp,'uint16');
    
    if cal.correctInt
        C = cal.Icorrf1;
    else
        C = uint16(ones(8,1));
        data1 = uint16(data1);
    end
    
    if cal.reorder
        newor = [cal.neworder1];
        
    else
        newor = 1:sTmp(3);
    end
  
    % correct int
    waitbar(.5,h,'Doing some simple math...')
    for i = 1:size(chC1,3)
        data1(:,:,i,:) = chC1(:,:,i,:)./C(i);
    end
    
    waitbar(.7,h,'Doing some simple math...')
    for i = 1:size(chC1,3)
        data1(:,:,i+size(chC1,3),:) = chC2(:,:,i,:)./C(i+size(chC1,3));
    end
    
    waitbar(.9,h,'Reordering...')
    % reorder planes
    data1 = data1(:,:,newor,:);
    isTransmission1 = isTransmission1(newor);
    close(h)

    data2 = [];

    if cal.multiModal == true
        h = waitbar(0,'Please wait applying calibration');
        waitbar(.2,h,'Gettingg channel data')
        [ chC3_0, chC4_0 ] = mpSetup.cali.getChData( cam1, cam2, cal.ROI2FullCam );
        for k = 1:size(chC3_0, 4)
            for z = 1:size(chC3_0, 3)
                tform = simtform2d(cal.Transformation{z,1}.Scale, cal.Transformation{z,1}.RotationAngle, cal.Transformation{z,1}.Translation);
                chC3(:,:,z,k) = uint16(imwarp(double(chC3_0(:,:,z,k)), tform, "OutputView", imref2d(size(double(chC1(:,:,z,k))))));
                chC4(:,:,z,k) = uint16(imwarp(double(chC4_0(:,:,z,k)), tform, "OutputView", imref2d(size(double(chC2(:,:,z,k))))));
            end
        end

        sTmp2 = size(chC3);
        sTmp2(3) = sTmp2(3)*2;
        data2 = ones(sTmp2,'uint16');
        
        if cal.correctInt
            C2 = cal.Icorrf2;
        else
            C2 = uint16(ones(8,1));
            data2 = uint16(data);
        end
        
        if cal.reorder
            newor2 = [cal.neworder2];
            
        else
            newor2 = 1:sTmp2(3);
        end
      
        % correct int
        waitbar(.5,h,'Doing some simple math...')
        for i = 1:size(chC3,3)
            data2(:,:,i,:) = chC3(:,:,i,:)./C2(i);
        end
        
        waitbar(.7,h,'Doing some simple math...')
        for i = 1:size(chC3,3)
            data2(:,:,i+size(chC3,3),:) = chC4(:,:,i,:)./C2(i+size(chC3,3));
        end
        
        waitbar(.9,h,'Reordering...')
        % reorder planes
        data2 = data2(:,:,newor2,:);
        isTransmission2 = isTransmission2(newor2);
        close(h)
    else 
    end
    data = {data1; data2};
    isTransmission = {isTransmission1; isTransmission2};
end
end
