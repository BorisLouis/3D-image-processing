clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = [];

%file info
Paths = {'D:\Polymer Dynamics\20251009\PAA_2x_AA', 'D:\Polymer Dynamics\20251009\PAA_3x_AA',...
    'D:\Polymer Dynamics\20251008\PAA_3x_bAA', 'D:\Polymer Dynamics\20251008\PAA_4x_bAA',...
    'D:\Polymer Dynamics\20251004\PAA_1x_bAA', 'D:\Polymer Dynamics\20251004\PAA_2x_bAA',...
    'D:\Polymer Dynamics\20250925\PAA_2x_bAA', 'D:\Polymer Dynamics\20250925\PAA_3x_bAA',...
    'D:\Polymer Dynamics\20250924\PAA_1x_bAA', 'D:\Polymer Dynamics\20250924\PAA_2x_bAA'};
Cals2D = {'D:\Polymer Dynamics\20251009\2DCal', 'D:\Polymer Dynamics\20251009\2DCal',...
    'D:\Polymer Dynamics\20251008\2DCal', 'D:\Polymer Dynamics\20251008\2DCal',...
    'D:\Polymer Dynamics\20251004\2DCal', 'D:\Polymer Dynamics\20251004\2DCal',...
    'D:\Polymer Dynamics\20250925\2DCal', 'D:\Polymer Dynamics\20250925\2DCal',...
    'D:\Polymer Dynamics\20250924\2DCal', 'D:\Polymer Dynamics\20250924\2DCal'};
file.path = Paths{1};
[info, info1, info2, file] = UserInput.infoGUI(file);
for j = 1:numel(Paths)
    file.path  = Paths{j};
    path2Cal = Cals2D{j};
    path2SRCal = Cals2D{j};    
    
    %% create experiments
    MultiModalExp = Core.MultiModalExperiment(file,path2Cal, info, info1, info2,path2SRCal,path2ZCal);
    
    %% get Movies
    MultiModalExp.RetrieveMovies;
end
%MultiModalExp.RunAnalysis;