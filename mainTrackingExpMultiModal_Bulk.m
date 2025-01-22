clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
%file.path  = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Multicolor Project\sample1_meaurement__2';
file.ext   = '.ome.tif';
path2Cal = 'D:\Dual Color\20250121\2DCal';
dimension = '3D';

MainFolder = 'D:\Dual Color';
SubFolders = {'20250121', '20250122'};
SubsubFolders = {'Multicolor_particles', 'PS_200_green_PS_100_red', 'PS_300_green_PS_100_red'};
SubsubsubFolders = {'0_min_measurements','In_water', 'sample_1', 'sample_2', 'sample_3'};
SubsubsubsubFolders = {'0_min1', '0_min2', '0_min3', '0_min4', '0_min5',...
    '3 min', '5 min', '7 min', '9 min', '11 min', '13 min', '15 min', '17 min'};

%detection parameter
detectParam{1}.delta = 6;
detectParam{1}.consThresh = 4;
detectParam{2}.delta = 6;
detectParam{2}.consThresh = 4;


%tracking parameter
trackParam.radius  = 2000;%nm
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
info.rotational = 0; %Rotational tracking 
info.rotationalCalib = 0;


for t = 1:numel(SubFolders)
    for r = 1:numel(SubsubFolders)
        for a = 1:numel(SubsubsubFolders)
            for c = 1:numel(SubsubsubsubFolders)
                file.path = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r},...
                    filesep, SubsubsubFolders{a}, filesep, SubsubsubsubFolders{c});

                if r == 1
                    detectParam{1}.chi2  = 25;
                    detectParam{2}.chi2  = 50;
                else
                    if c == 12
                        detectParam{1}.chi2  = 50;
                        detectParam{2}.chi2  = 30;
                    elseif c == 13
                        detectParam{1}.chi2  = 50;
                        detectParam{2}.chi2  = 30;
                    else
                        detectParam{1}.chi2  = 50;
                        detectParam{2}.chi2  = 50;
                    end
                end

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
                    trackingExp.ConsolidateChannels;
                    
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






        
