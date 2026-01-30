function [ totCor, Icor ] = findChIntPhase( cam1, cam2, inFocus, doFigure )
    %FINDCHINT finds hte intensity difference between the different channels
    %   Detailed explanation goes here
    
    switch nargin
        case 3
            doFigure = true;
    end
    
    N = size(cam1,4);
    
    for Frame = 1:size(cam1,4)
        for Plane = 1:size(cam1,3)
            I = cam1(:,:,Plane,Frame);
            threshSigma = 2;
            minArea = 6;
            se = strel('disk',5);
            bgImage = imopen(I, se);
            bg = median(bgImage(:));
            I_bgsub = I - bgImage;
            vals = I_bgsub(:);
            mu = median(vals);
            sigma = 1.4826 * median(abs(vals - mu));
            particleMask = I_bgsub > (mu + threshSigma * sigma);
            particleMask = bwareaopen(particleMask, minArea);
            bgCounts1(Frame, Plane) = median(I(~particleMask));
        end
    end
    
    for Frame = 1:size(cam2,4)
        for Plane = 1:size(cam2,3)
            I = cam2(:,:,Plane,Frame);
            threshSigma = 2;
            minArea = 6;
            se = strel('disk',5);
            bgImage = imopen(I, se);
            bg = median(bgImage(:));
            I_bgsub = I - bgImage;
            vals = I_bgsub(:);
            mu = median(vals);
            sigma = 1.4826 * median(abs(vals - mu));
            particleMask = I_bgsub > (mu + threshSigma * sigma);
            particleMask = bwareaopen(particleMask, minArea);
            I(particleMask) = nan;
            cam2(:,:,Plane,Frame) = I;
            bgCounts2(Frame, Plane) = median(I(~particleMask));
        end
    end
    I = double([bgCounts1, bgCounts2]);
       
    % to avoid problems with bleaching I can calculate Int differences between
    % channels by using the frame where both channels are almost in focus. 
    % global channel order
    globCh = cat(1,inFocus.globalch);
    % frame of sharp focus
    focus  = cat(1,inFocus.frame);
    % find the order of the list
    [~, idx] = sort(globCh);
    % get grames of sharp focus
    focus = focus(idx);
    % calculate the common frame to use
    dF = diff(focus);
    commonF = round(focus(1:end-1)+dF./2);
    
    % calculate intensity correction
    % Icor contains the reference channel, the channel to transform, the frame
    % used and the correction factor.
    Icor = ones(length(focus),4);
    Icor(1,1) = idx(1);
    Icor(1,2) = idx(1);
    % for i = 2:length(focus)
    %     % index of the first channel
    %     idx1 = idx(i-1);
    %     % index of the second channel
    %     idx2 = idx(i);
    %     % frame used 
    %     F    = commonF(i-1);
    %     % mean intensity values
    %     I1 = I(F, idx1);
    %     I2 = I(F, idx2);
    %     Icor(i,5) = I1;
    %     Icor(i,6) = I2;
    %     % now we store all values
    %     Icor(i,1) = idx1;
    %     Icor(i,2) = idx2;
    %     Icor(i,3) = F;
    %     % I must multiply by this factor to change int of channel idx2 into the
    %     % same level of idx1;
    %     Icor(i,4) = I1./I2;
    % end
    
    for i = 2:length(focus)
        % index of the first channel
        [I1, idx1] = min(I(:,i-1));
        % index of the second channel
        [I2, idx2] = min(I(:,i));
        % frame used 
        Icor(i,5) = I1;
        Icor(i,6) = I2;
        % now we store all values
        Icor(i,1) = i-1;
        Icor(i,2) = i;
        Icor(i,3) = idx2;
        % I must multiply by this factor to change int of channel idx2 into the
        % same level of idx1;
        Icor(i,4) = I1./I2;
    end

    % so far the corrections are between different channels, I have to find a
    % factor to correct all to the same refence. I do this by its cumulative
    % multiplication
    totCor = cumprod(Icor(:,4));
    % now I just have to order it correctly
    [~,newOrd] = sort(Icor(:,2));
    % sorting
    totCor = totCor(newOrd);
    
    if doFigure
        Itest = I.*repmat(totCor,1,N)';
    
    
        figure(1)
        subplot(1,2,1)
        plot(1:N,I(:,1))
        hold on
        for i=2:8
           plot(1:N,I(:,i)) 
        end
        hold off
        xlim([1, N])
        title('Before correction')
        xlabel('Frame')
        ylabel('Mean intensity of frame')
    
        subplot(1,2,2)
        plot(1:N,Itest(:,1))
        hold on
        for i=2:8
           plot(1:N,Itest(:,i)) 
        end
        hold off
        title('Corrected intensities')
        xlabel('Frame')
        ylabel('Mean intensity of frame')
    
        sgtitle('Corrected to plane 1')
    
    end

end

