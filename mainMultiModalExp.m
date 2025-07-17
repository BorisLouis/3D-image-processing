clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info

file.path  = 'C:\Users\steve\OneDrive\Documenten\TestData Indra\mSiPEI';
path2Cal = 'C:\Users\steve\OneDrive\Documenten\TestData Indra\2DCal';
[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;