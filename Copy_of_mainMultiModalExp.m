clc 
clear 
close all;
%file info
Paths = {'D:\Polymer Dynamics\20260121'};
Paths2DCal = {};
<<<<<<< HEAD
TimePaths = {'1x_bAA', '2x_AA', '3x_bAA', '3x_AA'};
=======
TimePaths = {'1x_bAA', '3x_bAA', '2x_AA', '3x_AA'};
>>>>>>> 963d6598c839f63c99ffb9d48fee704268afaea0
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        file.path  = append(Paths{j}, filesep, TimePaths{i});   
<<<<<<< HEAD
        path2Cal = 'D:\Polymer Dynamics\20260121\2DCal';
=======
        path2Cal = 'D:\Polymer Dynamics\20260121\2DCal_after';
>>>>>>> 963d6598c839f63c99ffb9d48fee704268afaea0
        path2ZCal = [];
        path2SRCal = [];
    
        %% create experiments
        MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
        
        %% get Movies
        MultiModalExp.RetrieveMovies;
        MultiModalExp.RunAnalysis;
    end
end
%