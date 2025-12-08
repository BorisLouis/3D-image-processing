clc ;
clear ;
close all;

MainMainFolders = {'D:\Data Steven - GEMs\20250725_GEMScarlet_37kPa_48h_KM12\KM12SM'};

[raw.FilePath, info.Experiment, info.FilenameRaw, info.Dimension, info.expTime, info.Temp, info.Radius1, info.Radius2, info.DiffFit, info.MinSize, info.Ext, info.ParticleType, info.path2RotCal, info.CutTraces, info.ExpModel] = UserInput.CalcMSDinfoGUI;

Microrheology = MicrorheologyAnalysis.Microrheology(raw, info);
Microrheology.setMovies;
Microrheology.RunAnalysis;