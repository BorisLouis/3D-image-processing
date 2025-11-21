clc 
clear 
close all;
%calibration info
path2ZCal = [];
% path2SRCal = 'S:\Dual Color\20251015_localisation_error\2DCal';

%file info
Paths = {'D:\DDM_TestData\PS_500nm_highconc', 'D:\DDM_TestData\PS_1000nm_highconc'};
file.path = Paths{1};
[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    file.path  = Paths{j};   
    path2Cal = [];
    path2SRCal = [];
    
    %% create experiments
    MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
    
    %% get Movies
    MultiModalExp.RetrieveMovies;
    MultiModalExp.RunAnalysis;
end
%