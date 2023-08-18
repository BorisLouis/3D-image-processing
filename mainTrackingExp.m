clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'D:\Documents\Unif\PhD\2021-Data\08 - August\Gold Particle code\1P data\22 nm silica';
file.ext   = '.ome.tif';
path2Cal = 'D:\Documents\Unif\PhD\2021-Data\08 - August\Gold Particle code\2DCal';
dimension = '3D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 60;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 1000;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.calibrate = false; %true to recalibrate;

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame =50;
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


