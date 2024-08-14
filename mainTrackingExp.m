clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'F:\multicolor_polarization\20240813_oil_MC_NPs_200nm\2D_Cal';

%file info
file.path  = 'F:\multicolor_polarization\20240813_oil_MC_NPs_200nm\200nm_PS_MC_NPs_in_water_1-25_dilution\sample1_meaurement__1';
file.ext   = '.ome.tif';
path2Cal = 'F:\multicolor_polarization\20240813_oil_MC_NPs_200nm\2D_Cal';
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

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

%[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data
trackingExp.saveData;


