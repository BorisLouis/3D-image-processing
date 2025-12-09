clc 
clear 
close all;
%file info
Paths = {'F:\Polymer Dynamics\20250924\PAA_1x_bAA', 'F:\Polymer Dynamics\20250924\PAA_2x_bAA',...
        'F:\Polymer Dynamics\20250925\PAA_2x_bAA', 'F:\Polymer Dynamics\20250925\PAA_3x_bAA'};
Paths2DCal = {'F:\Polymer Dynamics\20250924\2DCal', 'F:\Polymer Dynamics\20250924\2DCal',...
    'F:\Polymer Dynamics\20250925\2DCal', 'F:\Polymer Dynamics\20250925\2DCal'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    file.path  = Paths{j};   
    path2Cal = Paths2DCal{j};
    path2ZCal = [];
    path2SRCal = Paths2DCal{j};

    %% create experiments
    MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
    
    %% get Movies
    MultiModalExp.RetrieveMovies;
    MultiModalExp.RunAnalysis;
end
%