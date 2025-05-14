%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'E:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal';
file.ext  = '.ome.tif';
path2Cal  = 'E:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal';

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

testSRCal = Core.SRCalibration(file,path2Cal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 50;
detectParam{1}.consThresh = 6;
detectParam{2}.delta = 6;
detectParam{2}.chi2  = 50;
detectParam{2}.consThresh = 6;

trackParam.commonPlanes = 1; 
trackParam.euDistPx = 6;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibratewwwnBB
%% calc translation
refPlane = 4;
testSRCal.corrTranslation(refPlane);
testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);
testSRCal.checkAccuracy(refPlane);

%% calc channel translations
testSRCal.CalcAccuracyChannels(refPlane);