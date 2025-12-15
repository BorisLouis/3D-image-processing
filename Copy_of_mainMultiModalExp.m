clc 
clear 
close all;
%file info
Paths = {'D:\Data Hannah\20251211\200 nm PS beads'};
Paths2DCal = {''};
TimePaths = {'200nm_PS_phase_filter_480_60'};
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        file.path  = append(Paths{j}, filesep, TimePaths{i});   
        path2Cal = 'D:\Data Hannah\20251211\2D Cal';
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