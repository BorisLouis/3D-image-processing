clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
<<<<<<< HEAD

file.path  = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\Brightfield\500nmPS';
path2Cal = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\2D-cal';
=======
<<<<<<< HEAD
file.path  = 'D:\Data Hannah\20260203\2DCal_multimodal_before';
path2Cal = 'D:\Data Hannah\20260203\2DCal_multimodal_before';
=======
<<<<<<< HEAD
file.path  = 'D:\Data Steven - GEMs\test_laser_fluctuations';
path2Cal = [];
=======
file.path  = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\Brightfield\500nmPS';
path2Cal = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\2D-cal';
>>>>>>> 1927bfea7804a1814a3277c5899c3d031208d162
>>>>>>> 06f3aa7c780662434d34d2c94a4eb0296e7ae53f
>>>>>>> 09be1aa3dac6260aa0e3e35225f5da6883189da5

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;
disp('=== Analysis done ===')