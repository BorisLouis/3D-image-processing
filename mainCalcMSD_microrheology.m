clc ;
clear ;
close all;

[raw.FilePath, info.Experiment, info.FilenameRaw, info.Dimension, info.expTime, info.Temp, info.Radius1, info.Radius2, info.DiffFit, info.MinSize, info.Ext, info.ParticleType, info.path2RotCal, info.CutTraces, info.ExpModel, info.StepsizeAnalysis] = UserInput.CalcMSDinfoGUI;

MainFolder = 'F:\Polymer Dynamics';
SubFolders = {'PAA_1x_bAA', 'PAA_2x_bAA', 'PAA_3x_bAA', 'PAA_4x_bAA', 'PAA_2x_AA', 'PAA_3x_AA'};
SubSubFolders = {'time12'};
ExpTimes = [0.05];
CutTraces = [50];

for i = 1:numel(SubFolders)
    for j = 1:numel(SubSubFolders)
        raw.FilePath = append(MainFolder, filesep, SubFolders{i}, filesep, SubSubFolders{j});
        info.exptime = ExpTimes(j);
        info.CutTraces = CutTraces(j);
        info.FilenameRaw = "traces3D_";

        Microrheology = MicrorheologyAnalysis.Microrheology(raw, info);
        Microrheology.setMovies;
        Microrheology.RunAnalysis;
    end
end
