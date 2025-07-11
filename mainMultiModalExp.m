clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
file.path  = 'S:\Indra\20250710\Cam1\All data';
path2Cal = 'S:\Indra\20250710\Cam1\2DCal';
[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;