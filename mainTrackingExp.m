clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\20240513_spheric_PS_NPs_2Dcal_fluo\tracking_test';
file.ext   = '.ome.tif';
path2Cal = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\20240513_spheric_PS_NPs_2Dcal_fluo';
dimension = '3D';
multiModal = 'on'; %multiModal on or off

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 20;
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
info.detectionMethod = 'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame =60;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,1,frame);
testMov.showCandidate(frame, 1);
if multiModal == 'on'
    testMov.findCandidatePos(detectParam,2,frame);
    testMov.showCandidate(frame, 2);
end

%% get TrackingData


val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,1);
if multiModal == 'on'
    trackingExp.retrieveTrackData(detectParam,trackParam,2);
end

traces = trackingExp.getTraces3D(1);
if multiModal == 'on'
    traces2 = trackingExp.getTraces3D(2);
end



%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity(1);
if multiModal == 'on'
    [int2,SNR2] = trackingExp.getAvgIntensity(2);
end

%% Get MSD

%[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


