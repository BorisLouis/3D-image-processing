clc;
clear all;
close all;

%% Calibration info
path2ZCal = [];
path2SRCal = [];

%% Pathinfo
file.path = 'E:\MultiColor - lysosome tracking\20260311-new analysis\test';

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = false; %true to recalibrate;
file.ext   = '.his'; %extenstion of video
dimension = '2D';
info.PxSize = 81; %in nm %% ADJUST
info.FWHM = 3; %Half with of gauss
info.multiModal = 0; 
info.detectionMethod = 'Intensity'; %always take Intensity
info.ExpTime = 10;
MakeMovie = 1; %Make movie from traces (takes long)

%% Detection parameters
detectParam.size = [5 100]; % min and mix size of FAs in pixels %% ADJUST % VinCT: [5 120], VinCCM: [5 100], PaxCT&CCM: [5 100]
trackParam.radius = 3*info.PxSize; %nm % VinCT: 6, VinCCM: 4, PaxCT: 3
trackParam.memory = 5; %If trace is lost, how many frames to keep
info.file = file;

path2Cal = [];
path2SRCal = [];
path2ZCal = [];

%% get TrackingData
trackingExp = Core.TrackingExperimentNoFit(file,path2Cal,info,path2SRCal,path2ZCal);
trackingExp.retrieveMovies;
trackingExp.retrieveTrackData(detectParam,trackParam);
traces = trackingExp.getTraces3D;

if MakeMovie == 1
    %% get an amazing movie with the traces
    trackingExp.MakeMovie(10, diffInfo.minSize, 500, 1000);
end