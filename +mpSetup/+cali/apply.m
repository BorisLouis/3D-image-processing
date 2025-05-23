function [data,isTransmission, ROInew, BackgroundCorr] = apply( cam1, cam2, cal, info, ROI)
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
            nChan = size(cal.ROI,1)/2;
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
        waitbar(.5,h,'Correcting plane intensities...')
        for i = 1:size(chC1,3)
            data1(:,:,i,:) = chC1(:,:,i,:).*C(i);
        end
        
        waitbar(.7,h,'Correcting plane intensities...')
        for i = 1:size(chC1,3)
            data1(:,:,i+size(chC1,3),:) = chC2(:,:,i,:).*C(i+size(chC1,3));
        end
        
        waitbar(.9,h,'Reordering...')
        % reorder planes
        data1 = data1(:,:,newor,:);
        isTransmission1 = isTransmission1(newor);
        close(h)
    
        % if info.rotational == true
            %background substraction
            se = strel('disk', 12);
            h = waitbar(0,'Rotational tracking: background substraction...');
            for i = 1:size(data1,4)
                waitbar(i./size(data1,4),h,'Rotational tracking: background substraction...')
                for k = 1:size(data1,3)
                    frame = data1(:,:,k,i);
                    bg1(:,:,k,i) = imgaussfilt(frame, 15);
                end
            end
            close(h)
        % end

        BackgroundCorr.Ch1 = bg1;
    
        data2 = [];
    
        if cal.multiModal == true
            h = waitbar(0,'Channel2: Please wait applying calibration');
            waitbar(.2,h,'Channel2: Getting channel data')
            [ chC3, chC4 ] = mpSetup.cali.getChData( cam1, cam2, cal.ROI2FullCam );
     
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
            waitbar(.5,h,'Correcting plane intensities...')
            for i = 1:size(chC3,3)
                data2(:,:,i,:) = chC3(:,:,i,:).*C2(i);
            end
            
            waitbar(.7,h,'Correcting plane intensities...')
            for i = 1:size(chC3,3)
                data2(:,:,i+size(chC3,3),:) = chC4(:,:,i,:).*C2(i+size(chC3,3));
            end
            
            waitbar(.9,h,'Channel2: Reordering...')
            % reorder planes
            data2 = data2(:,:,newor2,:);
            isTransmission2 = isTransmission2(newor2);
            close(h)
    
            % if info.rotational == true
                %background substraction
                h = waitbar(0,'Channel2: Rotational tracking: background substraction...');
                se = strel('disk', 12);
                for i = 1:size(data2,4)
                    waitbar(i./size(data2,4),h,'Channel2: Rotational tracking: background substraction...')
                    for k = 1:size(data2,3)
                        frame = data2(:,:,k,i);
                        bg2(:,:,k,i) = imgaussfilt(frame, 15);
                        CorrBgFactor(k,i) = mean((bg1(:,:,k,i)./bg2(:,:,k,i)), 'all');
                    end
                end
                close(h)
            % end
            BackgroundCorr.Ch2 = bg2;
            BackgroundCorr.Ratio = CorrBgFactor;
        else 
        end
        data = {data1; data2};
        isTransmission = {isTransmission1; isTransmission2};
    
        if exist('ROInew','var') == 0
            ROInew = NaN;
        end
    end
end
