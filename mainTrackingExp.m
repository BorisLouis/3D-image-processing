clc 
clear 
close all;

path2ZCal = 'D:\Dropbox\4Apr Data\ZCal CS';
path2SRCal = 'D:\Dropbox\4Apr Data\2DCal';

file.path  = 'D:\Dropbox\4Apr Data\XYZ - OD35\X';
file.ext   = '.ome.tif';
path2Cal = 'D:\Dropbox\4Apr Data\2DCal';
dimension = '3D';
detectParam.delta = 6;
detectParam.chi2 = 60;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'PSFE'; %Intensity, 3DFit or PSFE
info.calibrate = false; %true to recalibrate;

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame = 20;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;

val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


