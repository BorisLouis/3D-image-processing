%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'G:\multicolor_polarization\polarisation\20241223_AuBPs_193x90_calib_polaris\2D_AuNPs';
file.ext  = '.ome.tif';
path2Cal  = 'G:\multicolor_polarization\polarisation\20241223_AuBPs_193x90_calib_polaris\2D_AuNPs';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod = 'Phasor';
info.zMethod   = 'Intensity';
info.multiModal = 1;
info.detectionMethod = 'MaxLR'; 
info.rotational = 0;

testSRCal = Core.SRCalibration(file,path2Cal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 40;
detectParam.consThresh = 4;

trackParam.commonPlanes = 1; 
trackParam.euDistPx = 3;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibratewwwnBB
%% calc translation
refPlane = 4;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);