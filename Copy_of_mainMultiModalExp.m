clc 
clear 
close all;
%file info
Paths = {'E:\Data Steven - GEMs\data paper'};
Paths2DCal = {};
TimePaths = {'KM12C', 'KM12L4a', 'KM12SM', 'SW480', 'SW620'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
            file.path  = append(Paths{j}, filesep, TimePaths{i});   
            path2Cal = [];
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