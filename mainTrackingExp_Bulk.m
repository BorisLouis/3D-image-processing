clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
MainFolder = 'S:\Dual Color\20250220\Test_stickyness_2D\PS100 green';
SizeFolder = {'sample1', 'sample2', 'sample3'};
SampleFolder = {'3_min', '5_min', '4_min', '6_min', '7_min', '8_min', '9_min', '10_min',....
                '11_min',  '12_min', '13_min', '14_min', '15_min', '16_min','17_min', '18_min'}; % 

file.ext   = '.his';
path2Cal = [];
dimension = '2D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 20;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 3500;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.detectionMethod = 'MaxLR'; %MaxLR (for maximum likehood ratio) %Intensity
info.calibrate = false; %true to recalibrate;
info.multiModal = 0; %multiModal (1) or not (0)
info.rotational = 0;
info.rotationalCalib = 0;
info.PxSize = 81; %in nm

for t = 1:numel(SizeFolder)
    for r = 1:numel(SampleFolder)
        %try
            file.path = append(MainFolder, filesep, SizeFolder{t}, filesep, SampleFolder{r});
    
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
            
            %% save Data
            trackingExp.saveData;
        % catch
        % end
    end
end

