clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
MainFolder = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';
subFolders = {'sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5', 'sample_6', 'sample_7', 'sample_8', 'sample_9', 'sample_10'};
file.ext   = '.ome.tif';
path2Cal = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal';
dimension = '3D';

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 40;
detectParam{1}.consThresh = 4;
detectParam{2}.delta = 6;
detectParam{2}.chi2  = 40;
detectParam{2}.consThresh = 4;

%tracking parameter
trackParam.radius  = 800;%nm
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
info.rotational = 1; %Rotational tracking 
info.rotationalCalib = 1;
info.euDist = 1500;
info.expTime = 0.010; %in sec
info.RadTime = 25; %in degrees per second (speed of rotating waveplate)
info.PxSize = 95;

for i = 1:size(subFolders, 2)
    % try
        file.path = append(MainFolder, filesep, subFolders{i});
    
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
        trackingExp.RotationalCalibration;
        trackingExp.saveData;
    % catch
    % end
end

