clc ;
clear ;
close all;

[raw.FilePath, info.Experiment, info.FilenameRaw, info.Dimension, info.expTime, info.Temp, info.Radius1, info.Radius2, info.DiffFit, info.MinSize, info.Ext, info.ParticleType, info.path2RotCal, info.CutTraces, info.ExpModel, info.StepsizeAnalysis] = UserInput.CalcMSDinfoGUI;


Microrheology = MicrorheologyAnalysis.Microrheology(raw, info);
Microrheology.setMovies;
Microrheology.RunAnalysis;

