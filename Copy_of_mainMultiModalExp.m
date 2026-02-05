clc 
clear 
close all;
%file info
Paths = {'D:\Data Hannah\20260203\200nm_in_PAA', 'D:\Data Hannah\20260203\200nm_in_PAA'};
Paths2DCal = {};
TimePaths = {'brightfield', 'darkfield'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
            file.path  = append(Paths{j}, filesep, TimePaths{i});   
            path2Cal = 'D:\Data Hannah\20260203\2DCal_multimodal_before';
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