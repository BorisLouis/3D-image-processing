clc 
clear 
close all;
%file info
<<<<<<< HEAD
Paths = {'E:\Data Steven - GEMs\data paper'};
Paths2DCal = {};
TimePaths = {'KM12C', 'KM12L4a', 'KM12SM', 'SW480', 'SW620'};
=======

Paths = {'E:\DDM_TestData\PS_100nm_highconc', 'E:\DDM_TestData\PS_200nm_highconc', 'E:\DDM_TestData\PS_500nm_highconc', 'E:\DDM_TestData\PS_1000nm_highconc',...
    'E:\DDM_TestData\PS_100nm_lowconc', 'E:\DDM_TestData\PS_200nm_lowconc', 'E:\DDM_TestData\PS_500nm_lowconc', 'E:\DDM_TestData\PS_1000nm_lowconc'};
Paths2DCal = {};
TimePaths = {'sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'};
>>>>>>> 06cf564fac4eed7ba0d7e886222a6806c64289cf
file.path = Paths{1};

[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
            file.path  = append(Paths{j}, filesep, TimePaths{i});   
<<<<<<< HEAD
            path2Cal = [];
=======
            path2Cal = []; %'D:\Data Hannah\20260203\2DCal_multimodal_before';
>>>>>>> 06cf564fac4eed7ba0d7e886222a6806c64289cf
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