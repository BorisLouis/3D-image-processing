function [cal, movInfo, MagnificationFactors] = calculate(fPath,nChan, DarkFieldPhase, Phase, correctInt, flipCam2)
%CALCULATE calculates calibration for the multiplane setup. NOTE that if
%you choose to correct for intensity differences the data changes form
%uint16 to double because we have to multiply by a correction factor
%(double).

switch nargin
    case 3
        
        nChan = 4;
        flipCam2 = true;
        
    case 4
        % this is the normal way, I cant see how it would be different
        correctInt = true;
        flipCam2 = true;
        
    case 5
        
        flipCam2 = true;
        
end

cal.correctInt = correctInt;


cal.reorder    = true;


% load general information about the calibration movie
[frameInfo, movInfo, ~ ]= Load.Movie.ome.getInfo( fPath );

% load the data 
[ movC1, movC2, idx ] = Load.Movie.ome.load( frameInfo, movInfo, 1:movInfo.maxFrame );

h = waitbar(0,'Please wait...');
% store information about position given by the motors
posInfo = cat(1,frameInfo.Pos);
Z1 = posInfo(idx(:,1),3);
Z2 = posInfo(idx(:,2),3);

% if need be we flip camera 2, this is generally the case
cal.flipCam2 = flipCam2;
if cal.flipCam2
    movC2 = flip(movC2,2);
end
meanIm1 = mean(movC1,3);
meanIm2 = mean(movC2,3);
if size(movC1,1)>1000
    multiModal = true;
    cal.multiModal = multiModal;
    [chan1,chan3, idx1] = mpSetup.cali.splitMultiModalCamera(meanIm1, 650);
    [chan2,chan4, idx2] = mpSetup.cali.splitMultiModalCamera(meanIm2, 650);
else 
    multiModal =false;
    cal.multiModal = multiModal;
    chan1 = meanIm1;
    chan2 = meanIm2;
end

waitbar(.1,h,'Finding channels')
% find channels
[ chCentCam1, ~, commonW1 ] = mpSetup.cali.findChannels( chan1, true, nChan, DarkFieldPhase );
[ chCentCam2, ~, commonW2 ] = mpSetup.cali.findChannels( chan2, true, nChan, DarkFieldPhase );
commonwin = min([commonW1; commonW2]); 
if multiModal == true
    [ chCentCam3, ~, commonW3 ] = mpSetup.cali.findChannels( chan3, true, nChan, DarkFieldPhase );
    [ chCentCam4, ~, commonW4 ] = mpSetup.cali.findChannels( chan4, true, nChan, DarkFieldPhase );
    commonwin = min([commonW1; commonW2; commonW3; commonW4]); 
else 
end

waitbar(.2,h,'getting ROIs')
% get ROI
imS = size(meanIm1);
[ cal.ROI1 ] = mpSetup.cali.defineROI( commonwin, chCentCam1, chCentCam2, imS );

waitbar(.3,h,'getting channel data')
% get data for each channel identified. chData has dim im_size1 im_size2
% 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI1, []);

waitbar(.4,h,'getting focus metric')
% getting the focus metric as used in EPFL we might want to change this to
% a gradient method.
[ cal.focusMet1, cal.inFocus1, cal.fit1 ] = mpSetup.cali.getFocusMetric( chData1c, chData2c , Z1, Z2, 0);
cal.Zpos = Z1;

waitbar(.5,h,'getting new order for channels')
% find the new order for the camera channels
[ cal.neworder1, cal.inFocus1 ] = mpSetup.cali.getNewOrder( cal.inFocus1 );

waitbar(.7,h,'getting image shifts')
% find image shift in order to have the same ROIs to a pixel resoltuon
[ cal.imShifts1 ] = mpSetup.cali.simpleImShift4( cal.inFocus1, chData1c, chData2c );

waitbar(.8,h,'refining ROIs')
% refine the ROIs to consider the shifts
[ cal.ROI1 ] = mpSetup.cali.refineROI2( cal.ROI1, cal.imShifts1 );

if multiModal == true
    [ cal.ROI2 ] = mpSetup.cali.defineROI( commonwin, chCentCam3, chCentCam4, imS );
    
    for i = 1:4
        cal.ROI2FullCam(i,2) = cal.ROI2(i,2)+idx1;
        cal.ROI2FullCam(i,1) = cal.ROI2(i,1);
        cal.ROI2FullCam(i,3) = cal.ROI2(i,3);
        cal.ROI2FullCam(i,4) = cal.ROI2(i,4);
        cal.ROI2FullCam(i+4,2) = cal.ROI2(i+4,2)+idx2;
        cal.ROI2FullCam(i+4,1) = cal.ROI2(i+4,1);
        cal.ROI2FullCam(i+4,3) = cal.ROI2(i+4,3);
        cal.ROI2FullCam(i+4,4) = cal.ROI2(i+4,4);
    end
    cal.cutCamerasMultiModal = [idx1, idx2];
    [ chData3c, chData4c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI2FullCam, []);
    
    [ cal.focusMet2, cal.inFocus2, cal.fit2 ] = mpSetup.cali.getFocusMetric( chData3c, chData4c , Z1, Z2, 1);
    [ cal.neworder2, cal.inFocus2 ] = mpSetup.cali.getNewOrder( cal.inFocus2 );

    [ MagnificationFactors ] = mpSetup.cali.GetMagnificationScale(chData1c, chData2c, chData3c, chData4c, cal.inFocus1, cal.inFocus2);
    [ cal.ROI2, cal.ROI2FullCam ] = mpSetup.cali.refineROItransf(cal.ROI2, cal.ROI2FullCam, MagnificationFactors, size(meanIm1), idx1, idx2);
    [ chData3c, chData4c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI2FullCam, []);

    [ cal.imShifts2 ] = mpSetup.cali.simpleImShift4( cal.inFocus2, chData3c, chData4c );  
    [ cal.ROI2 ] = mpSetup.cali.refineROI2( cal.ROI2, cal.imShifts2 );
    [ cal.ROI2FullCam ] = mpSetup.cali.refineROI2( cal.ROI2FullCam, cal.imShifts2 );
    
end

figure()
subplot(2,1,1)
imagesc(meanIm1)
axis image
for i = 1:nChan
    rectangle('Position', cal.ROI1(i,:))
end
if multiModal == true
    hold on
    for i = 1:nChan
    rectangle('Position', cal.ROI2FullCam(i,:))
    end
else
end
title('Camera 1 with ROIs')
hold on 

subplot(2,1,2)
imagesc(meanIm2)
axis image
for i = nChan+1:2*nChan
    rectangle('Position', cal.ROI1(i,:))
end
if multiModal == true
    hold on
    for i = nChan+1:2*nChan
    rectangle('Position', cal.ROI2FullCam(i,:))
    end
else
end
title('Camera 2 with ROIs')

mpSetup.cali.plotCal(meanIm1,meanIm2, cal.ROI1);
sgtitle('Planes 1-8')
if multiModal == true
    mpSetup.cali.plotCal(meanIm1, meanIm2, cal.ROI2FullCam);
    sgtitle('Planes 9-15')
else 
end

if cal.correctInt
    waitbar(.9,h,'Correcting intensity')
    % update the channel data
    [ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI1, []);
    % calculate intensity correction    
    if Phase == 1
        [ cal.Icorrf1, IntCh1] = mpSetup.cali.findChIntPhase( chData1c, chData2c, cal.inFocus1 );
    else
        [ cal.Icorrf1, IntCh1] = mpSetup.cali.findChInt( chData1c, chData2c, cal.inFocus1 );
    end
    cal.Icorrf1 = cal.Icorrf1./(cal.Icorrf1(4));
    % maxInt1 = max(cal.fit1(:,2:2:end),[],1);
    % cal.Icorrf1 = maxInt1./max(maxInt1);
    if multiModal == true
        [ chData3c, chData4c] = mpSetup.cali.getChData( movC1, movC2, cal.ROI2FullCam, []);
        if Phase == 1
            [cal.Icorrf2, IntCh2] = mpSetup.cali.findChIntPhase( chData3c, chData4c, cal.inFocus2);
            cal.Icorrf2 = cal.Icorrf2/(cal.Icorrf2(4));
            [cal.Icorrf2] = mpSetup.cali.findChIntChannels(cal.Icorrf1, cal.Icorrf2, IntCh1, IntCh2);
        else
            [cal.Icorrf2, IntCh2] = mpSetup.cali.findChInt( chData3c, chData4c, cal.inFocus2);
            cal.Icorrf2 = cal.Icorrf2/(cal.Icorrf2(4));
            [cal.Icorrf2] = mpSetup.cali.findChIntChannels(cal.Icorrf1, cal.Icorrf2, IntCh1, IntCh2);
        end
        mpSetup.cali.plotCalMultiModal(chData1c, chData2c, chData3c, chData4c, cal.inFocus1, cal.inFocus2)
        sgtitle('All planes corrected')
    else
        [ MagnificationFactors ] = [];
    end
else
    cal.Icorrf = ones(8,1);
end

close(h)
end