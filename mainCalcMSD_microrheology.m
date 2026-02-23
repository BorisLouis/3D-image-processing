clc ;
clear ;
close all;

[raw.FilePath, info.Experiment, info.FilenameRaw, info.Dimension, info.expTime, info.Temp, info.Radius1, info.Radius2, info.DiffFit, info.MinSize, info.Ext, info.ParticleType, info.path2RotCal, info.CutTraces, info.ExpModel, info.StepsizeAnalysis] = UserInput.CalcMSDinfoGUI;

MainFolder = 'E:\Polymer Dynamics';
SubFolders = {'20260121'}; %
SubSubFolders = {'2x_initiator'}; 
ExpTimes = [0.033];
% CutTraces = [50];

for i = 1:numel(SubFolders)
    for j = 1:numel(SubSubFolders)
        try
            raw.FilePath = append(MainFolder, filesep, SubFolders{i}, filesep, SubSubFolders{j});
            % info.exptime = 0.05;
            % info.CutTraces = 50;
            info.FilenameRaw = "traces3D_";
    
            Microrheology = MicrorheologyAnalysis.Microrheology(raw, info);
            Microrheology.setMovies;
            Microrheology.RunAnalysis;
        catch
        end
    end
end
