function [data,isTransmission, ROInew] = apply( cam1, cam2, cal, info, ROI)
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
    if info.rotational == true
        if exist('ROI','var') == 1
             StartHoriz = ROI(1);
             EndHoriz = ROI(2);
             StartVert = ROI(3);
             EndVert = ROI(4);
        else
            MeanIm = mean(chC1,4);
            MeanIm = mean(MeanIm, 3);
            MeanHoriz = mean(MeanIm, 2);
            MeanHorizDiff = abs(diff(MeanHoriz));
            StartHoriz1 = find(MeanHorizDiff < 5, 1, 'first');
            EndHoriz1 = find(MeanHorizDiff < 5, 1, 'last');
            MeanVert = mean(MeanIm, 1);
            MeanVertDiff = abs(diff(MeanVert));
            StartVert1 = find(MeanVertDiff < 5, 1, 'first');
            EndVert1 = find(MeanVertDiff < 5, 1, 'last');
            
            MeanIm = mean(chC2,4);
            MeanIm = mean(MeanIm, 3);
            MeanHoriz = mean(MeanIm, 2);
            MeanHorizDiff = abs(diff(MeanHoriz));
            StartHoriz2 = find(MeanHorizDiff < 5, 1, 'first');
            EndHoriz2 = find(MeanHorizDiff < 5, 1, 'last');
            MeanVert = mean(MeanIm, 1);
            MeanVertDiff = abs(diff(MeanVert));
            StartVert2 = find(MeanVertDiff < 5, 1, 'first');
            EndVert2 = find(MeanVertDiff < 5, 1, 'last');
            
            StartHoriz = max([StartHoriz1; StartHoriz2]);
            EndHoriz = min([EndHoriz1; EndHoriz2]);
            StartVert = max([StartVert1; StartVert2]);
            EndVert = min([EndVert1; EndVert2]);
            ROInew = [StartHoriz, EndHoriz, StartVert, EndVert];
        end
        chC1 = chC1(StartHoriz:EndHoriz, StartVert:EndVert, :,:);
        chC2 = chC2(StartHoriz:EndHoriz, StartVert:EndVert, :,:);
    end
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

    % if info.rotational == true
    %     %background substraction
    %     h = waitbar(0,'Rotational tracking: background substraction...');
    %     for i = 1:size(data1,4)
    %         waitbar(i./size(data1,4),h,'Rotational tracking: background substraction...')
    %         for k = 1:size(data1,3)
    %             MedianInt = median(data1(:,:,k,i), 'all');
    %             data1(:,:,k,i) = data1(:,:,k,i) - MedianInt;
    %         end
    %     end
    %     close(h)
    % end

    data2 = [];

    if cal.multiModal == true
        h = waitbar(0,'Channel2: Please wait applying calibration');
        waitbar(.2,h,'Channel2: Getting channel data')
        [ chC3, chC4 ] = mpSetup.cali.getChData( cam1, cam2, cal.ROI2FullCam );
        if info.rotational == true
            if exist('ROI','var') == 1
                 StartHorizM = ROI(1);
                 EndHorizM = ROI(2);
                 StartVertM = ROI(3);
                 EndVertM = ROI(4);
                 ROInew = ROI;
            else
                MeanIm = mean(chC3,4);
                MeanIm = mean(MeanIm, 3);
                MeanHoriz = mean(MeanIm, 2);
                MeanHorizDiff = abs(diff(MeanHoriz));
                StartHoriz3 = find(MeanHorizDiff < 5, 1, 'first');
                EndHoriz3 = find(MeanHorizDiff < 5, 1, 'last');
                MeanVert = mean(MeanIm, 1);
                MeanVertDiff = abs(diff(MeanVert));
                StartVert3 = find(MeanVertDiff < 5, 1, 'first');
                EndVert3 = find(MeanVertDiff < 5, 1, 'last');
                
                MeanIm = mean(chC4,4);
                MeanIm = mean(MeanIm, 3);
                MeanHoriz = mean(MeanIm, 2);
                MeanHorizDiff = abs(diff(MeanHoriz));
                StartHoriz4 = find(MeanHorizDiff < 5, 1, 'first');
                EndHoriz4 = find(MeanHorizDiff < 5, 1, 'last');
                MeanVert = mean(MeanIm, 1);
                MeanVertDiff = abs(diff(MeanVert));
                StartVert4 = find(MeanVertDiff < 5, 1, 'first');
                EndVert4 = find(MeanVertDiff < 5, 1, 'last');
                
                StartHorizM = max([StartHoriz3; StartHoriz4; StartHoriz]);
                EndHorizM = min([EndHoriz3; EndHoriz4; EndHoriz]);
                StartVertM = max([StartVert3; StartVert4; StartVert]);
                EndVertM = min([EndVert3; EndVert4; EndVert]);

                ROInew = [StartHorizM, EndHorizM, StartVertM, EndVertM];
                data1 = data1((StartHorizM-StartHoriz)+1:end-(EndHoriz-EndHorizM), (StartVertM-StartVert)+1:end-(EndVert-EndVertM),:,:);
            end
            chC3 = chC3(StartHorizM:EndHorizM, StartVertM:EndVertM, :,:);
            chC4 = chC4(StartHorizM:EndHorizM, StartVertM:EndVertM, :,:);
            
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
        waitbar(.5,h,'Channel2: Doing some simple math...')
        for i = 1:size(chC3,3)
            data2(:,:,i,:) = chC3(:,:,i,:)./C2(i);
        end
        
        waitbar(.7,h,'Channel2: Doing some simple math...')
        for i = 1:size(chC3,3)
            data2(:,:,i+size(chC3,3),:) = chC4(:,:,i,:)./C2(i+size(chC3,3));
        end
        
        waitbar(.9,h,'Channel2: Reordering...')
        % reorder planes
        data2 = data2(:,:,newor2,:);
        isTransmission2 = isTransmission2(newor2);
        close(h)

        if info.rotational == true
            %background substraction
            h = waitbar(0,'Channel2: Rotational tracking: background substraction...');
            for i = 1:size(data2,4)
                waitbar(i./size(data2,4),h,'Channel2: Rotational tracking: background substraction...')
                for k = 1:size(data2,3)
                    MedianInt = median(data2(:,:,k,i), 'all');
                    data2(:,:,k,i) = data2(:,:,k,i) - MedianInt;
                end
            end
            close(h)
        end
    else 
    end
    data = {data1; data2};
    isTransmission = {isTransmission1; isTransmission2};
end
end
