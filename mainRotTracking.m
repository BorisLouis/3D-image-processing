clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'G:\multicolor_polarization\polarisation\20241223_AuBPs_193x90_calib_polaris\2D_AuBPs_polarization\sample_1';
file.ext   = '.ome.tif';
path2Cal = 'G:\multicolor_polarization\polarisation\20241223_AuBPs_193x90_calib_polaris\2D_AuNPs';
dimension = '3D';

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 55;
detectParam{1}.consThresh = 4;

detectParam{2}.delta = 6;
detectParam{2}.chi2  = 40;
detectParam{2}.consThresh = 4;
%tracking parameter
trackParam.radius  = 1500;%nm
trackParam.memory  = 50;

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
info.euDist = 500;

%% create experiments
trackingExp = Core.TrackingExperimentRotational(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies
trackingExp.retrieveMovies;

%% test detection parameters
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam);
testMov.getROIs;
testMov.showCandidate(1);

%% get TrackingData
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);
traces = trackingExp.getTraces3D;
trackingExp.ConsolidateChannels2;
trackingExp.saveData;



