function [cal] = calculate(fPath, correctInt, flipCam2)
%CALCULATE calculates calibration for the multiplane setup. NOTE that if
%you choose to correct for intensity differences the data changes form
%uint16 to double because we have to multiply by a correction factor
%(double).

switch nargin
    case 1
        correctInt = true;
        flipCam2 = true;
    case 2
        % this is the normal way, I cant see how it would be different
        flipCam2 = true;
end

cal.correctInt = correctInt;


cal.reorder    = true;


% load general information about the calibration movie
[frameInfo, movInfo, ~ ]= loadMovie.ome.getInfo( fPath );

% load the data 
[ movC1, movC2, idx ] = loadMovie.ome.load( frameInfo, movInfo, 1:movInfo.maxFrame );

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

% max projection image
max_im1 = max(movC1,[],3);
max_im2 = max(movC2,[],3);

waitbar(.1,h,'Finding channels')
% find channels
[ chCentCam1, ~, commonW1 ] = mpSetup.cali.findChannels( max_im1 );
[ chCentCam2, ~, commonW2 ] = mpSetup.cali.findChannels( max_im2 );

waitbar(.2,h,'getting ROIs')
% get ROI
commonwin = min([commonW1; commonW2]);
imS = size(max_im1);
[ cal.ROI ] = mpSetup.cali.defineROI( commonwin, chCentCam1, chCentCam2, imS );

waitbar(.3,h,'getting channel data')
% get data for each channel identified. chData has dim im_size1 im_size2
% 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );

waitbar(.4,h,'getting focus metric')
% getting the focus metric as used in EPFL we might want to change this to
% a gradient method.
[ focusMet, cal.inFocus ] = mpSetup.cali.getFocusMetric( chData1c, chData2c , Z1, Z2 );

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
[ cal.neworder, cal.inFocus ] = mpSetup.cali.getNewOrder( cal.inFocus );

waitbar(.6,h,'getting image shifts')
% find image shift in order to have the same ROIs to a pixel resoltuon
[ imShifts ] = mpSetup.cali.simpleImShift( cal.inFocus, chData1c, chData2c );

waitbar(.7,h,'refining ROIs')
% refine the ROIs to consider the shifts
[ cal.ROI ] = mpSetup.cali.refineROI( cal.ROI, imShifts );

figure()
subplot(2,1,1)
imagesc(max_im1)
axis image
for i = 1:4
    rectangle('Position', cal.ROI(i,:))
end
title('Camera 1 with ROIs')
subplot(2,1,2)
imagesc(max_im2)
axis image
for i = 5:8
    rectangle('Position', cal.ROI(i,:))
end
title('Camera 2 with ROIs')


if cal.correctInt
    waitbar(.8,h,'Correcting intensity')
    % update the channel data
    [ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );
    % calculate intensity correction
    [ cal.Icorrf ] = mpSetup.cali.findChInt( chData1c, chData2c, cal.inFocus );
else
    cal.Icorrf = ones(8,1);
end

close(h)
end

% % %% TODO I still have to fix this part, however I do not think is needed for our purposes
% % 
% % % init output
% %     tf=cell(sys.nplanes-1,1);
% % 
% %     for m=(1:sys.nplanes-1)-1
% %         
% %         ch = sys.nplanes-m
% % 
% %         % get image sequence for 2 channels to compare
% %         img_seq_ch_fix=squeeze(data12(:,:,ch,:));
% %         img_seq_ch_moving=squeeze(data12(:,:,ch-1,:));
% % 
% %         % see when each channel is in focus % TODO this I can do better
% %         focus_metric_fix = squeeze(mean(max(img_seq_ch_fix)));
% %         focus_metric_moving = squeeze(mean(max(img_seq_ch_moving)));
% %         
% %         % get a frame when both cameras are simultaneously in best focus
% %         [~,indcom]=max(focus_metric_fix.*focus_metric_moving);
% %         im_fix=medfilt2(img_seq_ch_fix(:,:,indcom),[3 3]);
% %         im_mov=medfilt2(img_seq_ch_moving(:,:,indcom),[3 3]);
% %         
% %         % identify corresponding beads in the two channels
% %         [my, mx] = ccrShiftEstimation(im_fix,im_mov,10); % ! Stefan switched mx and my
% %         
% % 
% %         % extracting the center of gravities of the beads
% %         out=struct;
% %         out.ru=5;
% %         sys.bg=1;
% %         out.nh=0;
% % 
% %         out = hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
% %         out = hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
% %         cog_fix = [out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
% %         
% %         out = hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
% %         out = hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
% %         cog_mov = [out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
% % 
% %         %identify corresponding center of gravities
% %         
% %         pixel_tolerance = 2;
% %         mainSet = cog_fix - repmat([mx my],size(cog_fix,1),1);
% %         testSet = cog_mov;
% %         
% %         [ cog_common_fix, cog_common_mov ] = ...
% %              ImageProc.consolidatePos( mainSet, testSet, pixel_tolerance );
% %          
% %         d = cog_common_mov - cog_common_fix;
% %         
% %         disp(['avg error of COG coordinates using pure displacement at pixel level: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);
% %         
% %         cog_common_fix     = cog_common_fix + ...
% %                                  repmat([mx my],size(cog_common_fix,1),1);
% %         
% %         
% %         % affine transformation
% %         ccm=fliplr(cog_common_mov);
% %         ccf=fliplr(cog_common_fix);
% %         % Infer spatial transformation from control point pairs
% %         tf{m+1}=cp2tform(ccm,ccf,'similarity');
% %         
% % %         % I think this is here for example purposes
% % %         % Apply forward spatial transformation.
% % %         [x,y] = tformfwd(tf{m+1},cog_common_mov(:,2),cog_common_mov(:,1));
% %     end
% %     
% %     % tf{1} transforms ch7 into ch8
% %     % tf{2} transforms ch6 into ch7
% %     % tf{3} transforms ch5 into ch6
% %     % tf{4} transforms ch4 into ch5
% %     % tf{5} transforms ch3 into ch4
% %     % tf{6} transforms ch2 into ch3
% %     % tf{7} transforms ch1 into ch2
% %        
% %     cal.tf = tf;



