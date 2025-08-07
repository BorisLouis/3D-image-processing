clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\2DCal';

%file info
<<<<<<< HEAD
file.path  = 'D:\Multimodal tracking\20250724\alldata';
path2Cal = 'D:\Multimodal tracking\20250724\2DCal';
=======
file.path  = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\alldata';
path2Cal = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\2DCal';
>>>>>>> 789a2706a2455919fa94798e7f7a1f6f4c444a9a
[info, info1, info2, file] = UserInput.infoGUI(file);

%% create experiments
MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);

%% get Movies
MultiModalExp.RetrieveMovies;
MultiModalExp.RunAnalysis;