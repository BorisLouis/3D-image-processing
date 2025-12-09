clc 
clear 
close all;
%file info
Paths = {'F:\Polymer Dynamics\PAA_1x_bAA', 'F:\Polymer Dynamics\PAA_2x_bAA',...
        'F:\Polymer Dynamics\PAA_3x_bAA', 'F:\Polymer Dynamics\PAA_4x_bAA',...
        'F:\Polymer Dynamics\PAA_2x_AA', 'F:\Polymer Dynamics\PAA_3x_AA'};
Paths2DCal = {'F:\Polymer Dynamics\20251004\2DCal', 'F:\Polymer Dynamics\20251004\2DCal',...
    'F:\Polymer Dynamics\20251008\2DCal', 'F:\Polymer Dynamics\20251008\2DCal',...
    'F:\Polymer Dynamics\20251009\2DCal', 'F:\Polymer Dynamics\20251009\2DCal'};
TimePaths = {'time3', 'time6', 'time12'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        file.path  = append(Paths{j}, filesep, TimePaths{i});   
        path2Cal = Paths2DCal{j};
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