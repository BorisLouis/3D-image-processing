clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info

<<<<<<< HEAD
file.path  = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\alldata';
path2Cal = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\2DCal';
=======
file.path  = 'D:\Multimodal tracking\20250724\alldata';
path2Cal = 'D:\Multimodal tracking\20250724\2DCal';

>>>>>>> 9ec929933c8979e268d62c2493e4ffa72a6e8177
[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;