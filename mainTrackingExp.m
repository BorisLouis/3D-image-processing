clc 
clear 
close all;


path2ZCal = 'D:\Documents\Unif\PhD\2019-Data\Roger\02-Feb\26\3D Cal';
path2SRCal = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\2DCal';

path2File = 'D:\Documents\Unif\PhD\2019-Data\Roger\02-Feb\26\Tracking';
path2Cal = 'D:\Documents\Unif\PhD\2019-Data\Roger\02-Feb\26\2DCal\Crop';

detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,calib.getCal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData
detectParam.delta = 6;
detectParam.chi2 = 60;
trackParam.euDistXY = 400;
trackParam.euDistZ  = 400;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% save Data

trackingExp.saveData;


