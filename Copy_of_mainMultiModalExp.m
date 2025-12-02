clc 
clear 
close all;
%file info
Paths = {'E:\Data Hannah\20251113\500 nm polystyrene beads\500nm_closed_video', 'E:\Data Hannah\20251113\500 nm polystyrene beads\500nm_closed_z-stack',...
    'E:\Data Hannah\20251113\1000 nm polystyrene beads\1000nm_closed_video', 'E:\Data Hannah\20251113\1000 nm polystyrene beads\1000nm_closed_z-stack',...
    'E:\Data Hannah\20251113\2000 nm polysterene beads\2000nm_closed_video', 'E:\Data Hannah\20251113\2000 nm polysterene beads\2000nm_closed_z-stack'};
file.path = Paths{1};
[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    file.path  = Paths{j};   
    path2Cal = 'E:\Data Hannah\20251113\2DCal';
    path2ZCal = [];
    path2SRCal = [];

    %% create experiments
    MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
    
    %% get Movies
    MultiModalExp.RetrieveMovies;
    MultiModalExp.RunAnalysis;
end
%