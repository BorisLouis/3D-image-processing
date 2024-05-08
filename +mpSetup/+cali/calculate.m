function [cal, movInfo] = calculate(fPath,nChan, correctInt, flipCam2)
%CALCULATE calculates calibration for the multiplane setup. NOTE that if
%you choose to correct for intensity differences the data changes form
%uint16 to double because we have to multiply by a correction factor
%(double).

switch nargin
    case 1
        
        nChan = 4;
        flipCam2 = true;
        
    case 2
        % this is the normal way, I cant see how it would be different
        correctInt = true;
        flipCam2 = true;
        
    case 3
        
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
    [chan1,chan3, idx1] = mpSetup.cali.splitMultiModalCamera(meanIm1);
    [chan2,chan4, idx2] = mpSetup.cali.splitMultiModalCamera(meanIm2);
else 
    multiModal =false;
end

waitbar(.1,h,'Finding channels')
% find channels
% [ chCentCam1, ~, commonW1 ] = mpSetup.cali.findChannels( meanIm1, false, nChan );
% [ chCentCam2, ~, commonW2 ] = mpSetup.cali.findChannels( meanIm2, false, nChan );
[ chCentCam1, ~, commonW1 ] = mpSetup.cali.findChannels( chan1, false, nChan );
[ chCentCam2, ~, commonW2 ] = mpSetup.cali.findChannels( chan2, false, nChan );
if multiModal == true
    [ chCentCam3, ~, commonW3 ] = mpSetup.cali.findChannels( chan3, false, nChan );
    [ chCentCam4, ~, commonW4 ] = mpSetup.cali.findChannels( chan4, false, nChan );
else
end

waitbar(.2,h,'getting ROIs')
% get ROI
commonwin = min([commonW1; commonW2; commonW3; commonW4]);
imS = size(meanIm1);
if multiModal == false
[ cal.ROI ] = mpSetup.cali.defineROI( commonwin, chCentCam1, chCentCam2, imS );
elseif multiModal == true
    [ cal.ROI ] = mpSetup.cali.defineROIMultiModal( commonwin, chCentCam1, chCentCam2, chCentCam3, chCentCam4, imS );
else
end

waitbar(.3,h,'getting channel data')
% get data for each channel identified. chData has dim im_size1 im_size2
% 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );
if multiModal == true
    [ chData3c, chData4c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );
else
end

waitbar(.4,h,'getting focus metric')
% getting the focus metric as used in EPFL we might want to change this to
% a gradient method.
if multiModal == false
[ cal.focusMet, cal.inFocus, cal.fit ] = mpSetup.cali.getFocusMetric( chData1c, chData2c , Z1, Z2 );
elseif multiModal == true
    [ cal.focusMet, cal.inFocus, cal.fit ] = mpSetup.cali.getFocusMetricMultiModal( chData1c, chData2c, chData3c, chData4c, Z1, Z2 );
end
cal.Zpos = Z1;

% figure(2)
% for i = 1:4
%     plot(Z1,focusMet(:,i))
%     hold on
% end
% for i = 5:8
%     plot(Z2,focusMet(:,i))
% end
% hold off

waitbar(.5,h,'getting new order for channels')
% find the new order for the camera channels
if multiModal == false
[ cal.neworder, cal.inFocus ] = mpSetup.cali.getNewOrder( cal.inFocus );
elseif multiModal == true
[ cal.neworder, cal.inFocus ] = mpSetup.cali.getNewOrderMultiModal( cal.inFocus );
end

waitbar(.7,h,'getting image shifts')
% find image shift in order to have the same ROIs to a pixel resoltuon
if multiModal == false
[ imShifts ] = mpSetup.cali.simpleImShift( cal.inFocus, chData1c, chData2c );
elseif multiModal == true
[ imShifts ] = mpSetup.cali.simpleImShiftMultiModal( cal.inFocus, chData1c, chData2c, chData3c, chData4c );
end

waitbar(.8,h,'refining ROIs')
% refine the ROIs to consider the shifts
[ cal.ROI ] = mpSetup.cali.refineROI( cal.ROI, imShifts );

if multiModal == true
    for i = 9:12
    cal.ROI(i,2) = cal.ROI(i,2)+idx1
    cal.ROI(i+4,2) = cal.ROI(i+4,2)+idx2
    end
end


figure()
subplot(2,1,1)
imagesc(meanIm1)
axis image
for i = 1:nChan
    rectangle('Position', cal.ROI(i,:))
end
if multiModal == true
    hold on
    for i = (nChan*2+1):(nChan*3)
    rectangle('Position', cal.ROI(i,:))
    end
else
end
title('Camera 1 with ROIs')
hold on 

subplot(2,1,2)
imagesc(meanIm2)
axis image
for i = nChan+1:2*nChan
    rectangle('Position', cal.ROI(i,:))
end
if multiModal == true
    hold on
    for i = (nChan*3+1):(nChan*4)
    rectangle('Position', cal.ROI(i,:))
    end
else
end
title('Camera 2 with ROIs')

mpSetup.cali.plotCal(meanIm1,meanIm2, cal.ROI);

if cal.correctInt
    waitbar(.9,h,'Correcting intensity')
    % update the channel data
    [ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );
    % calculate intensity correction
    [ cal.Icorrf ] = mpSetup.cali.findChInt( chData1c, chData2c, cal.inFocus );
    maxInt = max(cal.fit(:,2:2:end),[],1);
    cal.Icorrf = maxInt./max(maxInt);
else
    cal.Icorrf = ones(8,1);
end

close(h)
end


