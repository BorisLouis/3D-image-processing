clc 
clear 
close all;
%file info
Paths = {'D:\Data Hannah\202512117\Darkfield'};
Paths2DCal = {};
TimePaths = {'500nmPS', '1000nmPS', '2000nmPS'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        file.path  = append(Paths{j}, filesep, TimePaths{i});   
        path2Cal = 'D:\Data Hannah\202512117\2D-cal';
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