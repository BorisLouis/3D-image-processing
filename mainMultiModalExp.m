clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'E:\Data Steven\test_3D_multiplane_2channel\2DCal';

%file info
file.path  = 'E:\Data Steven\test_3D_multiplane_2channel\data';
path2Cal = 'E:\Data Steven\test_3D_multiplane_2channel\2DCal';

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;
disp('=== Analysis done ===')