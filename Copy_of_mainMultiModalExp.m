clc 
clear 
close all;
%file info
Paths = {'D:\Polymer Dynamics\20251016'};
Paths2DCal = {'D:\Data Hannah\202512117\2D-cal', 'D:\Data Hannah\202512117\2D-cal'};
TimePaths = {'PAA_2x_APS_TEMED', 'PAA_3x_AA'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        file.path  = append(Paths{j}, filesep, TimePaths{i});   
        path2Cal = 'D:\Polymer Dynamics\20251016\2DCal';
        path2ZCal = [];
        path2SRCal = 'D:\Polymer Dynamics\20251016\2DCal';
    
        %% create experiments
        MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
        
        %% get Movies
        MultiModalExp.RetrieveMovies;
        MultiModalExp.RunAnalysis;
    end
end
%