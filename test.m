MainFolder = 'S:\Dual Color';
DateFolders = {'20250121', '20250122'};
ParticleFolders = {'Multicolor_particles'};
ExperimentFolders = {'In_water'};
SampleFolders = {'0_min1', '0_min2', '0_min3', '0_min4', '0_min4'};



for a = 1:numel(DateFolders)
    ExcelName = append(MainFolder, filesep, 'Results_MC_InWater_', DateFolders{a}, '.xlsx');
    for z = 1:numel(ParticleFolders)
        for e = 1:numel(ExperimentFolders)
            for r = 1:numel(SampleFolders)
                filename = append(MainFolder, filesep, DateFolders{a}, filesep, ParticleFolders{z},...
                    filesep, ExperimentFolders{e}, filesep, SampleFolders{r});

                %%% Channel 1
                MSDRes = load(append(filename, filesep, 'msdRes1.mat'));
                MSDRes = MSDRes.allRes;
                Diff = [];
                Visc = [];
                AnExp = [];
                Rg = [];
                for t = 1:size(MSDRes, 2)
                    Diff(t,1) = MSDRes(t).DR;
                    Visc(t,1) = MSDRes(t).nR;
                    AnExp(t,1) = MSDRes(t).aR;
                    Rg(t,1) = MSDRes(t).Rg;
                end
                
                if r == 1
                    Range = "A5:A5000";
                elseif r == 2
                    Range = "B5:B5000";
                elseif r == 3
                    Range = "C5:C5000";
                elseif r == 4
                    Range = "D5:D5000";
                elseif r == 5
                    Range = "E5:E5000";
                end
                writematrix(Diff, ExcelName, 'Sheet', 'Diff', 'Range', Range);
                writematrix(AnExp, ExcelName, 'Sheet', 'AnExp', 'Range', Range);
                writematrix(Visc, ExcelName, 'Sheet', 'Visc', 'Range', Range);
                writematrix(Rg, ExcelName, 'Sheet', 'Rg', 'Range', Range);


                %%% Channel 2
                MSDRes = load(append(filename, filesep, 'msdRes2.mat'));
                MSDRes = MSDRes.allRes;
                Diff = [];
                Visc = [];
                AnExp = [];
                Rg = [];
                for t = 1:size(MSDRes, 2)
                    Diff(t,1) = MSDRes(t).DR;
                    Visc(t,1) = MSDRes(t).nR;
                    AnExp(t,1) = MSDRes(t).aR;
                    Rg(t,1) = MSDRes(t).Rg;
                end
                
                if r == 1
                    Range = "F5:F5000";
                elseif r == 2
                    Range = "G5:G5000";
                elseif r == 3
                    Range = "H5:H5000";
                elseif r == 4
                    Range = "I5:I5000";
                elseif r == 5
                    Range = "J5:J5000";
                end
                writematrix(Diff, ExcelName, 'Sheet', 'Diff', 'Range', Range);
                writematrix(AnExp, ExcelName, 'Sheet', 'AnExp', 'Range', Range);
                writematrix(Visc, ExcelName, 'Sheet', 'Visc', 'Range', Range);
                writematrix(Rg, ExcelName, 'Sheet', 'Rg', 'Range', Range);

            end
        end
    end
end