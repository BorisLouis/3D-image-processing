clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'G:\multicolor_polarization\polarisation\20241202_AuBPs\2D_Cal\150x50_1-10_spincoating\sample_4';
file.ext   = '.ome.tif';
path2Cal = 'G:\multicolor_polarization\polarisation\20241202_AuBPs\2D_Cal\planes';
dimension = '3D';

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 55;
detectParam{1}.consThresh = 4;

detectParam{2}.delta = 6;
detectParam{2}.chi2  = 60;
detectParam{2}.consThresh = 4;
%tracking parameter
trackParam.radius  = 5000;%nm
trackParam.memory  = 50;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 1; %multiModal (1) or not (0)
info.rotational = 1; %Rotational tracking 1 => bg substraction,...
info.euDist = 1000;

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



