%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'S:\Dual Color\20250121_dualcolor\2DCal';
file.ext  = '.ome.tif';
path2Cal  = 'S:\Dual Color\20250121_dualcolor\2DCal';

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

testSRCal = Core.SRCalibrationMultiModal(file,path2Cal,info);

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
trackParam.euDistPx = 10;

testSRCal.SRAnalysis(detectParam, trackParam);