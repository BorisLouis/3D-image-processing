clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'E:\DDM_TestData\PS_100nm\sample1';
file.ext   = '.his';
path2Cal = [];
dimension = '2D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 40;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 3500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 0; %multiModal (1) or not (0)
info.rotational = 0;
info.PxSize = 81; %in nm

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame = 50;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData


val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

%[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

%[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

%trackingExp.showTraces(1);
%% save Data
trackingExp.saveData;


