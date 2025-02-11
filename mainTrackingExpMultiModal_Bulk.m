clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
%file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\sample1_meaurement__2';
file.ext   = '.ome.tif';
path2Cal = 'D:\Rotational Tracking\20250211\2D_cal';
dimension = '3D';

MainFolder = 'D:';
SubFolders = {'Rotational Tracking'};
SubsubFolders = { '20250211'};
SubsubsubFolders = {'2DCal_AuBPs_193x90_rotation'};
SubsubsubsubFolders = {'sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8', 'sample9', 'sample10'}; %, 

%detection parameter
detectParam{1}.chi2  = 35; 
detectParam{1}.delta = 6;
detectParam{1}.consThresh = 4;
detectParam{2}.chi2  = 30;
detectParam{2}.delta = 6;
detectParam{2}.consThresh = 4;


%tracking parameter
trackParam.radius  = 1500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.euDist = 1000; %Error distance between particles in different channels
info.multiTracking = 'Rotation'; %MultiColor or Rotation
info.rotational = 1; %Rotational tracking 
info.rotationalCalib = 1;


for t = 1:numel(SubFolders)
    for r = 1:numel(SubsubFolders)
        for a = 1:numel(SubsubsubFolders)
            for c = 1:numel(SubsubsubsubFolders)
                file.path = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r},...
                    filesep, SubsubsubFolders{a}, filesep, SubsubsubsubFolders{c});

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
                    % if r == 1
                    %     trackingExp.ConsolidateChannels3;
                    % end
                    
                    %% Get Intensity
                    [int,SNR] = trackingExp.getAvgIntensity;
                    
                    %% save Data
                    trackingExp.saveData;
        
        
                catch
                    continue
                end
            end
        end
    end
end






        
