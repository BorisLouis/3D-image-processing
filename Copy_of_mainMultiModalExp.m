clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'S:\Dual Color\20251015_localisation_error\2DCal';

%file info
Paths = {'S:\Dual Color\20251015_localisation_error\10ms_expTime', 'S:\Dual Color\20251015_localisation_error\500ms_expTime'};
path2Cal = 'S:\Dual Color\20251015_localisation_error\2DCal';
file.path = Paths{1};
[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    file.path  = Paths{j};   
    
    %% create experiments
    MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
    
    %% get Movies
    MultiModalExp.RetrieveMovies;
end
%MultiModalExp.RunAnalysis;