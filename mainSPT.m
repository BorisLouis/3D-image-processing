% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

% ad general paths that can be usefull
addpath(genpath('Ext'));

% path to the callibration
filePath = '..\data\Multiplane\PlaneCalib\BeadsCalibrationZStack_1';
fName = 'BeadsCalibrationZStack_1_MMStack_Pos0.ome.tif';

fullfilePath = [filePath filesep fName];

% Calculate calibration
[cal, info] = mpSetup.cali.calculate(fullfilePath, false);
calib.path = filePath;
calib.info = info;
calib.file = cal;
disp('Done with calibration')
%%
% load and calibrate, when applied to the calibration data then we should
% be able to demonstrate that it works

folderPath = '..\data\Multiplane\Data\TL-OD2-200msExposure_1';
%fName = 'TL-OD2-200msExposure_1_MMStack_Pos0.ome.tif';
%filePath = [folderPath filesep fName];

% % load general information about the multi-plane movie
% [~, movInfo, ~ ]= Load.Movie.ome.getInfo( filePath );
% 
% raw.path = filePath;
% raw.info = movInfo;
raw = folderPath;
testMov = Core.Movie(raw,[],calib.path);
%%
frame = 1:movInfo.maxFrame(1);
[data, frameInfo, movInfo] = mpSetup.loadAndCal( testMov.raw.Path, testMov.cal.file, frame);
step = 100;
fid = fopen([folderPath filesep 'CalibratedInfo.txt'],'w');
fprintf(fid,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
for i = 1:size(data,3)
    
    fName = sprintf('calibratedPlane%d.tif',i);
    fPathTiff = [folderPath filesep fName];
    t = Tiff(fPathTiff, 'w');
    for j = 1:step:size(data,4)
        range = j:j+step-1;
        if max(range)>= size(data,4)
            range = j:size(data,4);
        end
    t = dataStorage.writeTiff(t,squeeze(data(:,:,i,range)),16);
    
    end
    t.close;
    fprintf(fid,...
        'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
        i,testMov.cal.file.inFocus(testMov.cal.file.neworder==i).cam,...
        testMov.cal.file.inFocus(testMov.cal.file.neworder==i).ch,...
        testMov.cal.file.ROI(testMov.cal.file.neworder==i,1),...
        testMov.cal.file.ROI(testMov.cal.file.neworder==i,1)+...
        testMov.cal.file.ROI(testMov.cal.file.neworder==i,3),...
        testMov.cal.file.inFocus(testMov.cal.file.neworder==i).zpos-...
        testMov.cal.file.inFocus(testMov.cal.file.neworder==1).zpos);
        
end
fclose(fid);


%% example of a frame list I will grow this into the frame object
frameList = mcodekit.list.dl_list();
for i = 1:8
    tmp = data(:,:,i,1);
    imP = Core.imPlane(tmp);
    imP.setPixSizeNm(100);
    imP.setTime(uint16(1)); 
    frameList.append_key(imP);   
end

%% detect particles
% for a water immersion obj the best-fit gasuss to the PSF has 
% sigma = 0.25 wavelength / NA
objNA  = 1.2;
emWave = 600;
sigma_nm = 0.25 * emWave/objNA;
FWHMnm = sigma_nm * sqrt(8*log(2));         

GLRTprops.delta  = 6;
GLRTprops.pxSnm  = 100;
GLRTprops.FWHMnm = FWHMnm;
GLRTprops.chi2   = 80;

rTh = 5; % in pixels
ROIrad = 10;

pList = [];
p = [];
partList = mcodekit.list.dl_list();

for fIdx = frame(1:10)
    
    sfData = data(:,:,:,fIdx);
    imSize = size(sfData);
    % detect particles
    [consLoc,totLoc] = mpSetup.localize(sfData, rTh, GLRTprops);
    % build ROIs
    [ROIs] = Misc.getROIs(consLoc,ROIrad,imSize);
    
    for i = 1:size(consLoc,1)
        tmpLoc = consLoc(i,:);
        tmpROI = ROIs(i,:);
        p = Core.particle(tmpLoc,fIdx,tmpROI,sfData);

        pList = [pList, p];
        partList.append_key(p);
        
    end
    disp(['done for frame ' num2str(fIdx)])

end

%%
objNA    = 1.2;
emWave   = 600;
pxSizeNm = 100;
tic
pList.setPSFprops(objNA, emWave, pxSizeNm);
toc

tic
iterator = partList.get_iterator(); 
            
while (iterator.has_next())

    pTmp = iterator.next();
    pTmp.setPSFprops(objNA, emWave, pxSizeNm);
%     idx = idx + 1;

end
toc

%%
ptest = partList.get_key(10);
%%

tic
pList.superResolve;
disp('Done with SR-loc')
toc

tic
iterator = partList.get_iterator(); 
            
while (iterator.has_next())

    pTmp = iterator.next();
    pTmp.superResolve();
%     idx = idx + 1;

end
toc

ptest1 = pList(2);
ptest2 = partList.get_key(2);

%%
tic
test = findobj(pList,'frame',1);
toc

tic
iterator = partList.get_iterator(); 
ftotal = zeros(partList.size_,1);
idx = 0;
while (iterator.has_next())

    idx = idx + 1;
    pTmp = iterator.next();
    ftotal(idx) = pTmp.frame;

end

idxList = find(ftotal==1);
toc
bla = [];
for idx = idxList'
    idx
    tmp = partList.get_key(idx);
    bla = [bla tmp];
end
%%
% test = findobj(pList,'frame',1);
% tmpVal = cat(1,test.superResLoc);
% tmpX = tmpVal(:,1);
% tmpY = tmpVal(:,2);
% tmpZ = tmpVal(:,3);
% % scatter3 (tmpX,tmpY, tmpZ)
% scatter (tmpX,tmpY,'kx')
% 
% % axis image
% shg

cols = {'k','r','g','b','y'};

for ii = 1:250
    
    test = findobj(pList,'frame',ii);
    cidx = ceil(ii/50);
    tmpVal = cat(1,test.superResLoc);
    
    tmpX = tmpVal(:,1);
    tmpY = tmpVal(:,2);
    tmpZ = tmpVal(:,3);
%     scatter3 (tmpX,tmpY, tmpZ)
    subplot(1,2,1)
    scatter (tmpX,tmpY,[cols{cidx} 'x'])
    hold on
    subplot(1,2,2)
    scatter (tmpX,tmpY,[cols{cidx} 'x'])
    hold on
end
subplot(1,2,1)
hold off
axis image
d = .6;
c = [270.5,243.6];
xlim([c(1)-d c(1)+d])
ylim([c(2)-d c(2)+d])
% ylim([242 243])

subplot(1,2,2)
hold off
axis image
d = .6;
c = [430.3,277.2];
xlim([c(1)-d c(1)+d])
ylim([c(2)-d c(2)+d])
shg
%%
% subplot(1,2,1)
% xlim([269.5,271.5])
% ylim([243, 244])

% make a common list of sm detections? should I have a test for seeing a
% molecule in at least 3 (or X) planes? once I have the common list I have
% to cropt the ROIs [dx, dy, 8] and do the fine localization. We will for
% the moment pick the best as I do not have a super-resolved registration
% matrix between all channels. This should be generated in order to use
% info from multiple planes to increase fit accuracy. it is from z tacks of
% beads that we can create such a registration.
