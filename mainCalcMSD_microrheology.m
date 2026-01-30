clc ;
clear ;
close all;

[raw.FilePath, info.Experiment, info.FilenameRaw, info.Dimension, info.expTime, info.Temp, info.Radius1, info.Radius2, info.DiffFit, info.MinSize, info.Ext, info.ParticleType, info.path2RotCal, info.CutTraces, info.ExpModel, info.StepsizeAnalysis] = UserInput.CalcMSDinfoGUI;

MainFolder = 'D:\Polymer Dynamics';
<<<<<<< HEAD
SubFolders = {'PAA_1x_bAA'}; %, 'PAA_2x_bAA', 'PAA_3x_bAA', 'PAA_4x_bAA', 'PAA_2x_AA', 'PAA_3x_AA'
SubSubFolders = {'time12\test'}; %, , 'time6', 'time6', 'time3'
ExpTimes = [0.05, 0.01, 0.01];
CutTraces = [50, NaN, NaN];
=======
SubFolders = {'20260121'};
SubSubFolders = {'2x_AA', '3x_AA'};
ExpTimes = [0.05];
CutTraces = [50];
>>>>>>> 1927bfea7804a1814a3277c5899c3d031208d162

for i = 1:numel(SubFolders)
    for j = 1:numel(SubSubFolders)
        raw.FilePath = append(MainFolder, filesep, SubFolders{i}, filesep, SubSubFolders{j});
        info.exptime = 0.05;
        info.CutTraces = 50;
        info.FilenameRaw = "traces3D_";

        Microrheology = MicrorheologyAnalysis.Microrheology(raw, info);
        Microrheology.setMovies;
        Microrheology.RunAnalysis;
    end
end
