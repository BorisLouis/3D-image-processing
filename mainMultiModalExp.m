clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info

<<<<<<< HEAD
file.path  = 'D:\Data Steven - GEMs\20250725_GEMScarlet_37kPa_48h_KM12\KM12SM';
=======
file.path  = 'C:\Users\steve\OneDrive\Documenten\data multicolor\testdata Hannah\500 nm';
>>>>>>> db0cfca37ae0e23ff0cf171a44caeb95b904c125
path2Cal = [];

[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;