clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
<<<<<<< HEAD
file.path  = 'D:\Data Steven - GEMs\test_laser_fluctuations';
path2Cal = [];
=======
file.path  = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\Brightfield\500nmPS';
path2Cal = 'C:\Users\steve\data_no_onedrive\Data Hannah\202512117\2D-cal';
>>>>>>> 1927bfea7804a1814a3277c5899c3d031208d162

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;
disp('=== Analysis done ===')