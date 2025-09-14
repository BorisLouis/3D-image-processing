clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'D:\Polymer Dynamics\20250911\2DCal';

%file info

file.path  = 'D:\Polymer Dynamics\20250911\AllData\PAA_93kPa_normal';
path2Cal = 'D:\Polymer Dynamics\20250911\2DCal';

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;