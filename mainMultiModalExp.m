clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info

file.path  = 'D:\Data Steven - GEMs\data paper\KM12C';
path2Cal = [];

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;