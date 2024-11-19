clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\Testdata_rotational\20241115_AuBPs_2DCal\phaseplate_phaseimaging\sample_2';
file.ext   = '.ome.tif';
path2Cal = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\Testdata_rotational\20241115_AuBPs_2DCal\2DCal_200nm_PS';
dimension = '3D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 40;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 1; %multiModal (1) or not (0)
info.rotational = 1; %Rotational tracking 1 => bg substraction,...
info.euDist = 5000;

%% create experiments
trackingExp = Core.TrackingExperimentMultiModal(file,path2Cal,info,path2SRCal,path2ZCal);

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
trackingExp.ConsolidateChannels2;
trackingExp.saveData;



