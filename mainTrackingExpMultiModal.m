clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'F:\multicolor_polarization\20240527_PAA_NPs\PS NPs 300 nm fluo\10 min';
file.ext   = '.ome.tif';
path2Cal = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\20240513_spheric_PS_NPs_2Dcal_fluo';
dimension = '3D';

%detection parameter
detectParam.delta = 4;
detectParam.chi2  = 100;
detectParam.consThresh = 1;
%tracking parameter
trackParam.radius  = 500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 0; %first check the 'normal' planes no multiModal (0)

%% create experiments
trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies
trackingExp.retrieveMovies;

%% test detection parameters
frame =850;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);
traces = trackingExp.getTraces3D;

%% Get Intensity
[int,SNR] = trackingExp.getAvgIntensity;

%% show traces
trackingExp.showTraces(1);

%% save Data
trackingExp.saveData;




%%% Repeat the whole procedure for the multimodal images
info.multiModal = 1; %Now check the 'multiModal' planes (1)

%% create experiments
trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies
trackingExp.retrieveMovies;

%% test detection parameters
frame =60;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);
traces = trackingExp.getTraces3D;

%% Get Intensity
[int,SNR] = trackingExp.getAvgIntensity;

%% show traces
trackingExp.showTraces(1);

%% save Data
trackingExp.saveData;
