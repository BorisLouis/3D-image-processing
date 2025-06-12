clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';

%file info
file.path  = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\Glycerol\glycerol 80\sample1';
file.ext   = '.ome.tif';
path2Cal = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\2DCal';
dimension = '3D';

[info, info1, info2] = UserInput.infoGUI();

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;