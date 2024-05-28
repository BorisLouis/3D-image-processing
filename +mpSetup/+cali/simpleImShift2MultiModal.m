function [ finalShifts ] = simpleImShift2MultiModal( inFocus1, inFocus2, cam1, cam2, cam3, cam4)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    inFocus = table2struct([struct2table(inFocus1); struct2table(inFocus2)]);
    nPlanes = length(inFocus);
    imShifts = zeros(nPlanes-1,2);
    focus = inFocus(1).frame;
    
    imCh1 = cam1(:,:,1,focus);
    % shift must be small thus
    if nPlanes == 2
        maxShift = 300;
    else
        maxShift = 300; %normally 30
    end%[pixels]
    % generating mask, so I only look at a local max around shift 0 to
    % max_shift
    bw = false(size(imCh1).*2 - [1 1]);
    bw_c = round(size(bw)./2);
    bw(bw_c(1),bw_c(2)) = true;
    se = strel('disk',maxShift,8);
    bw = imdilate(bw,se);
    
    for chIdx = 1:nPlanes-1
        %get the idx to the plane in order (1-2, 2-3,3-4...)
        idxPlane1_1 = [inFocus1.globalch]==chIdx;
        idxPlane2_1 = [inFocus2.globalch]==chIdx;
        idxPlane1_2 = [inFocus1.globalch]==chIdx+1;
        idxPlane2_2 = [inFocus2.globalch]==chIdx+1;
        %get the plane where there both equally focus from the mean of
        %their focus and extract the camera
        focus_1     = round((inFocus(idxPlane1_1).frame+inFocus(idxPlane2_1).frame)/2);
        focus_2     = round((inFocus(idxPlane1_2).frame+inFocus(idxPlane2_2).frame)/2);
        camIdx1_1 = inFocus1(idxPlane1_1).cam;
        camIdx2_1 = inFocus1(idxPlane2_1).cam;
        camIdx1_2 = inFocus2(idxPlane1_2).cam;
        camIdx2_2 = inFocus2(idxPlane2_2).cam;

        if camIdx1_1 == 1 && camIdx1_2 == 3;
            imCh1 = cam1(:,:,inFocus1(idxPlane1_1).ch,focus);
            imCh3 = cam1(:,:,inFocus2(idxPlane1_2).ch,focus);
        else
            imCh1 = cam2(:,:,inFocus1(idxPlane1_1).ch,focus);
            imCh3 = cam2(:,:,inFocus2(idxPlane1_2).ch,focus);
        end
         
        if camIdx2_1==1 && camIdx2_2 == 3;
            imCh2 = cam1(:,:,inFocus1(idxPlane2_1).ch,focus);
            imCh4 = cam1(:,:,inFocus2(idxPlane2_2).ch,focus);
        else
            imCh2 = cam2(:,:,inFocus(idxPlane2_1).ch,focus);
            imCh4 = cam2(:,:,inFocus(idxPlane2_2).ch,focus);
        end
        
        % cross correlation
        res    = normxcorr2(imCh1,imCh2);
        % applying the mask
        res    = res.*bw;
        % finding shift
%         [mrow,mcol] = find(res==max(res(:)));
%         mrow = mrow-size(imCh1,1);
%         mcol = mcol-size(imCh1,2);
        [~,m_corr] = max(res(:));
        [mrow,mcol]=ind2sub(size(res),m_corr);
        mrow=-mrow+size(imCh1,1);
        mcol=-mcol+size(imCh1,2);
         % storing the shifts
        imShifts(chIdx,:) = [mrow mcol];
        
       
        
    end
    % calculate shift for a certain reference plane
    refPlane = round(nPlanes/2);
    tmpImShifts = zeros(nPlanes,2);
    idx = 1:nPlanes;
    idx(idx==4) = [];
    tmpImShifts(idx,:) = imShifts;
    for i = 1: size(tmpImShifts,1)
        
        if i<refPlane
            tmpImShifts(i,:) = -sum(imShifts(i:refPlane-1,:),1);
        elseif i>refPlane
            tmpImShifts(i,:) = +sum(imShifts(refPlane:i-1,:),1);
        else
            %refPlane so we do nothing
        end
        
    end
    %reorder shifts
    %imShifts([inFocus.globalch],:) = tmpImShifts;
    finalShifts = tmpImShifts([inFocus.globalch],:);
end

