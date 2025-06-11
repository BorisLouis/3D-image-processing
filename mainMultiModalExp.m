clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';

%file info
file.path  = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\Glycerol\glycerol 80\sample1';
file.ext   = '.ome.tif';
path2Cal = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';
dimension = '3D';

%% Storing info for both channels
info.PxSize = 95;
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.Channel1 = 'Rotational Tracking'; %Translational tracking, rotational tracking, segmentation, phase
info.Channel2 = 'Rotational Tracking'; %Translational tracking, rotational tracking, segmentation, phase
info.FWHM = 3;

%% Storing info for channel 1
%%% for tracking
info1.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info1.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info1.detectionMethod = 'MaxLR';%'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info1.calibrate = false; %true to recalibrate;
info1.euDist = 1000; %Error distance between particles in different channels
info1.multiTracking = 'MultiColor'; %MultiColor or Rotation
info1.rotational = 1; %Rotational tracking 
info1.rotationalCalib = 0;
info1.testFrame = 95;
info1.detectParam.delta = 6;
info1.detectParam.chi2  = 50;
info1.detectParam.consThresh = 4;

info1.trackParam.radius  = 2500;%nm
info1.trackParam.memory  = 3;

%% Storing info for channel 2
%%% for tracking
info2.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info2.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info2.detectionMethod = 'MaxLR';%'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info2.calibrate = false; %true to recalibrate;
info2.euDist = 1000; %Error distance between particles in different channels
info2.multiTracking = 'MultiColor'; %MultiColor or Rotation
info2.rotational = 1; %Rotational tracking 
info2.rotationalCalib = 0;
info2.testFrame = 95;
info2.detectParam.delta = 6;
info2.detectParam.chi2  = 65;
info2.detectParam.consThresh = 4;
info2.trackParam.radius  = 2500;%nm
info2.trackParam.memory  = 3;

%% create experiments
trackingExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
trackingExp.RetrieveMovies;
trackingExp.RunAnalysis;


%% get TrackingData
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);
traces = trackingExp.getTraces3D;
trackingExp.ConsolidateChannels;

%% Get Intensity
[int,SNR] = trackingExp.getAvgIntensity;

%% show traces
trackingExp.showTraces;

%% save Data
trackingExp.saveData;