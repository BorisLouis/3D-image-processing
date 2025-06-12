clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';

%file info
file.path  = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\Glycerol\glycerol 80\sample1';
file.ext   = '.ome.tif';
path2Cal = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';
dimension = '3D';

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 30;
detectParam{1}.consThresh = 4;

detectParam{2}.delta = 6;
detectParam{2}.chi2  = 30;
detectParam{2}.consThresh = 4;
%tracking parameter
trackParam.radius  = 3500;%nm
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
info.rotational = 1; %Rotational tracking 
info.rotationalCalib = 0;
info.euDist = 1500;
info.expTime = 0.010; %in sec
info.RadTime = []; %in degrees per second (speed of rotating waveplate)
info.PxSize = 95;

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
trackingExp.ConsolidateChannels3;
trackingExp.saveData;



