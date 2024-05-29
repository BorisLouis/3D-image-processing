function [ finalShifts ] = simpleImShift2Test( inFocus1, inFocus2, cam1, cam2, cam3, cam4)
%SIMPLEIMSHIFT fast im shift calculation at the pixel resolution
%   spatial cross-correlation algorithm to determine shift of coordinates
    nPlanes = length(inFocus2) + length(inFocus2);
    imShifts = zeros(nPlanes-1,2);
    focus = inFocus1(1).frame;
    
    imCh1 = cam1(:,:,1,focus);

    % shift must be small thus
    if nPlanes == 2
        maxShift = 300;
    else
        maxShift = 30;
    end%[pixels]
    % generating mask, so I only look at a local max around shift 0 to
    % max_shift
    bw = false(size(imCh1).*2 - [1 1]);
    bw_c = round(size(bw)./2);
    bw(bw_c(1),bw_c(2)) = true;
    se = strel('disk',maxShift,8);
    bw = imdilate(bw,se);
    
    for chIdx = 1:(nPlanes./2)-1
        %get the idx to the plane in order (1-2, 2-3,3-4...)
        idxPlane1 = [inFocus1.globalch]==chIdx;
        idxPlane2 = [inFocus1.globalch]==chIdx+1;
        idxPlane3 = [inFocus2.globalch]==chIdx;
        idxPlane4 = [inFocus2.globalch]==chIdx+1;
        %get the plane where there both equally focus from the mean of
        %their focus and extract the camera
        focus12     = round((inFocus1(idxPlane1).frame+inFocus1(idxPlane2).frame)/2);
        focus34     = round((inFocus2(idxPlane3).frame+inFocus2(idxPlane4).frame)/2);
        camIdx1 = inFocus1(idxPlane1).cam;
        camIdx2 = inFocus1(idxPlane2).cam;
        camIdx3 = inFocus2(idxPlane3).cam;
        camIdx4 = inFocus2(idxPlane4).cam;

        if camIdx1==1
            imCh1 = cam1(:,:,inFocus1(idxPlane1).ch,focus12);
        else
            imCh1 = cam2(:,:,inFocus1(idxPlane1).ch,focus12);
        end
         
        if camIdx2==1
            imCh2 = cam1(:,:,inFocus1(idxPlane2).ch,focus12);
        else
            imCh2 = cam2(:,:,inFocus1(idxPlane2).ch,focus12);
        end

        if camIdx3==1
            imCh3 = cam3(:,:,inFocus2(idxPlane3).ch,focus34);
        else
            imCh3 = cam4(:,:,inFocus2(idxPlane3).ch,focus34);
        end

        if camIdx4==1
            imCh4 = cam3(:,:,inFocus2(idxPlane4).ch,focus34);
        else
            imCh4 = cam4(:,:,inFocus2(idxPlane4).ch,focus34);
        end
        
        % cross correlation
        res12    = normxcorr2(imCh1,imCh2);
        res13    = normxcorr2(imCh1,imCh3);
        res24    = normxcorr2(imCh2,imCh4);
        res34    = normxcorr2(imCh3,imCh4);
        % applying the mask
        res12    = res12.*bw;
        res13    = res13.*bw;
        res24    = res24.*bw;
        res34    = res34.*bw;
        % finding shift
%         [mrow,mcol] = find(res==max(res(:)));
%         mrow = mrow-size(imCh1,1);
%         mcol = mcol-size(imCh1,2);
        [~,m_corr12] = max(res12(:));
        [~,m_corr13] = max(res13(:));
        [~,m_corr24] = max(res24(:));
        [~,m_corr34] = max(res34(:));

        [mrow12,mcol12]=ind2sub(size(res12),m_corr12);
        [mrow13,mcol13]=ind2sub(size(res13),m_corr13);
        [mrow24,mcol24]=ind2sub(size(res24),m_corr24);
        [mrow34,mcol34]=ind2sub(size(res34),m_corr34);

        mrow12=-mrow12+size(imCh1,1);
        mrow13=-mrow13+size(imCh1,1);
        mrow24=-mrow24+size(imCh1,1);
        mrow34=-mrow34+size(imCh1,1);

        mcol12=-mcol12+size(imCh1,2);
        mcol13=-mcol13+size(imCh1,2);
        mcol24=-mcol24+size(imCh1,2);
        mcol34=-mcol34+size(imCh1,2);

        %Check if shifts make sense:
        diff1 = (abs(mrow12)+abs(mrow13)) - (abs(mrow34)+abs(mrow24));
        diff2 = (abs(mcol12)+abs(mcol13)) - (abs(mcol34)+abs(mcol24));
        % if diff1 > 1 || diff2 > 1
        %     error("correlations between planes 1-2 and 3-4 are not equal")
        % else
        % end
         % storing the shifts
        imShifts12(chIdx,:) = [mrow12 mcol12];
        imShifts13(chIdx,:) = [mrow13 mcol13];
        imShifts24(chIdx,:) = [mrow24 mcol24];
        imShifts34(chIdx,:) = [mrow34 mcol34];  
        
    end
    % calculate shift for a certain reference plane
    refPlane = round(nPlanes/4);
    tmpImShifts = zeros(nPlanes,2);
    idx = 1:nPlanes./2;
    idx(idx==4) = [];
    tmpImShifts1(idx,:) = [imShifts12];
    tmpImShifts2(idx,:) = [imShifts34];
    for i = 1:7
        if i<refPlane
            tmpImShifts1(i,:) = -sum(imShifts12(i:refPlane-1,:),1);
            tmpImShifts2(i,:) = -sum(imShifts34(i:refPlane-1,:),1);
        elseif i>refPlane
            tmpImShifts1(i,:) = +sum(imShifts12(refPlane:i-1,:),1);
            tmpImShifts2(i,:) = +sum(imShifts34(refPlane:i-1,:),1);
        else
            %refPlane so we do nothing
        end
        
    end

    %reorder shifts
    %imShifts([inFocus.globalch],:) = tmpImShifts;
    tmpImShifts2 = tmpImShifts2 + imShifts13(refPlane,:);
    
    finalShifts1 = tmpImShifts1([inFocus1.globalch],:);
    finalShifts2 = tmpImShifts2([inFocus2.globalch],:);
    finalShifts = [finalShifts1; finalShifts2];

end