clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'D:\Documents\2024 - Data\11 - November\LAIA\2D Cal';

%file info
file.path  = 'D:\Documents\2024 - Data\11 - November\LAIA\mov1';
file.ext   = '.ome.tif';
path2Cal = 'D:\Documents\2024 - Data\11 - November\LAIA\2D Cal';
dimension = '3D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 25;
detectParam.consThresh = 6;
%tracking parameter
trackParam.radius  = 2000;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame =10;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData


val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

%[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


