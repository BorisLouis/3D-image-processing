clc 
clear 
close all;
%file info
Paths = {'D:\Polymer Dynamics\20260121', 'D:\Polymer Dynamics\20260120'};
Paths2DCal = {};
TimePaths = {'1x_bAA', '2x_bAA', '3x_bAA', '2x_AA', '3x_AA'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
            file.path  = append(Paths{j}, filesep, TimePaths{i});   
            path2Cal = 'D:\Polymer Dynamics\20260121\2DCal_after';
            path2ZCal = [];
            path2SRCal = [];
        
            %% create experiments
            MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
            
            %% get Movies
            MultiModalExp.RetrieveMovies;
            MultiModalExp.RunAnalysis;
        catch
        end
    end
end
%