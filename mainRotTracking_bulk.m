clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';

%file info
file.ext   = '.ome.tif';
path2Cal = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';
dimension = '3D';

%path info
MainFolder = 'E:\Rotational Tracking';
SubFolders = {'20250407_AuBPs_184s92_glycerol'};
SubsubFolders = {'Glycerol'};
SubsubsubFolders = {'glycerol 95', 'glycerol 100'}; %, 'glycerol_85', , 
SubsubsubsubFolders = {'sample1','sample2', 'sample3', 'sample4','sample5'}; %, 'sample1', 'sample2', 'sample3','sample5'

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.chi2  = 40;
detectParam{1}.consThresh = 4;
detectParam{2}.delta = 6; %High for rotational tracking
detectParam{2}.chi2  = 40;
detectParam{2}.consThresh = 4;

%tracking parameter
trackParam.radius  = 3500;%nm
trackParam.memory  = 15;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %'Intensity'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 1; %multiModal (1) or not (0)
info.rotational = 1; %Rotational tracking 
info.rotationalCalib = 0;
info.euDist = 1000;
info.expTime = 0.010; %in sec
info.RadTime = []; %in degrees per second (speed of rotating waveplate)
info.PxSize = 95;

for t = 1:numel(SubFolders)
    for r = 1:numel(SubsubFolders)
        for a = 1:numel(SubsubsubFolders)
            for c = 1:numel(SubsubsubsubFolders)
                file.path = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r},...
                    filesep, SubsubsubFolders{a}, filesep, SubsubsubsubFolders{c});
                try
                    %% create experiments
                    trackingExp = Core.TrackingExperimentRotational(file,path2Cal,info,path2SRCal,path2ZCal);
                    
                    %% get Movies
                    trackingExp.retrieveMovies;
                    
                    %% test detection parameters
                    frame = 211;
                    testMov = trackingExp.trackMovies.mov1;
                    testMov.findCandidatePos(detectParam,frame);
                    testMov.getROIs;
                    testMov.showCandidate(frame);
                    
                    %% get TrackingData
                    val2Use = 'bestFocus';
                    trackingExp.retrieveTrackData(detectParam,trackParam);
                    traces = trackingExp.getTraces3D;
                    trackingExp.ConsolidateChannels3;
                    
                    %% save Data
                    trackingExp.saveData;
                catch
                end
            end
        end
    end
end