clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
<<<<<<< HEAD
file.path  = 'E:\DDM_TestData\PS_100nm\sample1';
file.ext   = '.his';
path2Cal = [];
dimension = '2D';
=======
file.path  = 'D:\Dual Color\20250121\PS_200_green_PS_100_red\sample2\15_min';
file.ext   = '.ome.tif';
path2Cal = 'D:\Dual Color\20250121\2DCal';
dimension = '3D';
>>>>>>> 26cb929d6524975666fa866e33b73b324ffeb515

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 50;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 2500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'%'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.euDist = 1000; %Error distance between particles in different channels
<<<<<<< HEAD
info.multiModal = 0;
info.multiTracking = 0; %MultiColor or Rotation
=======
info.multiTracking = 'MultiColor'; %MultiColor or Rotation
info.rotational = 0; %Rotational tracking 
info.rotationalCalib = 0;
>>>>>>> 26cb929d6524975666fa866e33b73b324ffeb515

%% create experiments
trackingExp = Core.TrackingExperimentMultiModal(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies
trackingExp.retrieveMovies;

%% test detection parameters
frame = 23;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);


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