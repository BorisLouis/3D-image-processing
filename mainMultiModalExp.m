clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'S:\Dual Color\20250121_dualcolor\2DCal';

%file info

file.path  = 'S:\Dual Color\20250121_dualcolor\Multicolor_particles\In_water';
path2Cal = 'S:\Dual Color\20250121_dualcolor\2DCal';

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;