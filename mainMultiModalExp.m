clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'E:\Data Steven\testdata yuzu\2d cal';

%file info
file.path  = 'E:\Data Steven\testdata yuzu\data';
path2Cal = 'E:\Data Steven\testdata yuzu\2d cal';

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;
disp('=== Analysis done ===')