%%
clear
clc
close all;
%% get path to SRCalibration
file.path = 'E:\Data Steven\testdata yuzu\2d cal';
file.ext  = '.ome.tif';
path2Cal  = 'E:\Data Steven\testdata yuzu\2d cal';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod = 'Phasor';
info.zMethod   = 'Intensity';
info.multiModal = 1;
info.detectionMethod = 'MaxLR'; 
info.rotational = 0;
info.PxSize = 95;
info.rotationalCalib = 0;
info.FWHM = 3;
info.Dimension = '3D';
info.IntCorr = 'off';
info.Channel1 = 'Translational Tracking';
info.Channel2 = 'Translational Tracking';
info.AlignChannels = 'off';
info.IntCorr = 'on';

testSRCal = Core.SRCalibrationMultiModal(file,path2Cal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 50;
detectParam{1}.consThresh = 4;
detectParam{1}.fitting = 'on';
detectParam{2}.delta = 6;
detectParam{2}.chi2  = 40;
detectParam{2}.consThresh = 15;
detectParam{2}.fitting = 'on';

trackParam.commonPlanes = 1; 
trackParam.euDistPx = 20;

testSRCal.SRAnalysis(detectParam, trackParam);