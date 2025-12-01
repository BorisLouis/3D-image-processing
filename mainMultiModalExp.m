clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
<<<<<<< HEAD

file.path  = 'E:\Data Hannah\20251105\test';
path2Cal = 'E:\Data Hannah\20251113\2DCal';
=======
file.path  = 'D:\Data Steven - GEMs\test';
path2Cal = [];
>>>>>>> 04400a85e93e248320c57d7a9eb1d753095b1273

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;
disp('=== Analysis done ===')