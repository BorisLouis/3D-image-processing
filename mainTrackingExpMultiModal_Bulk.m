clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
%file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\sample1_meaurement__2';
file.ext   = '.ome.tif';
path2Cal = 'G:\multicolor_polarization\multicolor\20241105_polymerisation_dual_color_PS_air_obj\2D_cal';
dimension = '3D';

MainFolder = 'G:\multicolor_polarization\multicolor\20241105_polymerisation_dual_color_PS_air_obj';
SubFolders = {'sample1', 'sample2', 'sample3'};
SubsubFolders = {'3 min', '5 min', '7 min', '9 min', '11 min', '13 min', '15 min'};

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 40;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 2500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'%'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.euDist = 1000; %Error distance between particles in different channels
info.multiTracking = 'MultiColor'; %MultiColor or Rotation


for m = 1:numel(SubFolders)
    SubFolder = SubFolders{m};
    for u = 1:numel(SubsubFolders)
        SubsubFolder = SubsubFolders{u};
        file.path = append(MainFolder, filesep, SubFolder, filesep, SubsubFolder);
        try
            %% create experiments
            trackingExp = Core.TrackingExperimentMultiModal(file,path2Cal,info,path2SRCal,path2ZCal);
            
            %% get Movies
            trackingExp.retrieveMovies;
            
            %% test detection parameters
            frame = 23;
            testMov = trackingExp.trackMovies.mov1;
            testMov.findCandidatePos(detectParam,frame);
            testMov.showCandidate(frame);
            
            
            %% get TrackingData
            val2Use = 'bestFocus';
            trackingExp.retrieveTrackData(detectParam,trackParam);
            traces = trackingExp.getTraces3D;
            %trackingExp.ConsolidateChannels;
            
            %% Get Intensity
            [int,SNR] = trackingExp.getAvgIntensity;
            
            %% save Data
            trackingExp.saveData;


        catch
            continue
        end
    end
end