classdef Microrheology < handle
    %CALCMSD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        raw
        info
        Movies
    end
    
    methods
        function obj = Microrheology(raw, info)
            obj.raw  = raw;
            obj.info = info;
        end
        
        function setMovies(obj)
             MainFolder = dir(obj.raw.FilePath);
             disp('Unpacking movie folders')
             n = 0;
             for j = 3 : size(MainFolder,1)
                 if MainFolder(j).isdir == 1
                     try
                         n = n+1;
                         Path = append(MainFolder(j).folder, filesep, MainFolder(j).name);
                         if contains(obj.info.Experiment, 'Rotational')
                             Movie = MicrorheologyAnalysis.RotationalTracking(Path, obj.info);
                         else
                             Movie = MicrorheologyAnalysis.TranslationalTracking(Path, obj.info);
                         end
                         obj.Movies.(append("Movie", num2str(n))) = Movie;
                     catch
                     end
                 end
             end
        end

        function RunAnalysis(obj)
             Movies = fieldnames(obj.Movies);
             Results1 = [];
             Results2 = [];
             for j = 1 : size(Movies, 1)
                 try
                     CurrentMovie = obj.Movies.(Movies{j});
                     if strcmp(obj.info.Experiment, "Rotational Tracking")
                        CurrentMovie.LoadTraces(obj.info.FilenameRaw);
                        CurrentMovie.Analysis;
                     elseif strcmp(obj.info.Experiment, "Dual color tracking")
                        for loop = 1:2
                            Radius = obj.info.(append("Radius", num2str(loop)));
                            FileName = append(obj.info.FilenameRaw, num2str(loop));
                            CurrentMovie.LoadTraces(FileName);
                            if strcmp(obj.info.StepsizeAnalysis, 'on')
                                CurrentMovie.CalculateStepsizes(loop);
                                CurrentMovie.FitDiffPopulations(loop);
                                if loop == 1
                                    Results1(end+1,1) = CurrentMovie.PopulationFractions.D(1,1);
                                elseif loop == 2
                                    Results2(end+1,1) = CurrentMovie.PopulationFractions.D(1,1);
                                end
                                if ~isnan(obj.info.CutTraces)
                                    [Results] = CurrentMovie.FitPopulationFractions(loop);
                                    % if loop == 1
                                    %     Results1 = [Results1; Results];
                                    % elseif loop == 2
                                    %     Results2 = [Results2; Results];
                                    % end
                                end
                            else
                                CurrentMovie.TracesAnalysis(Radius, loop);
                            end
                        end
                     else
                         Radius = obj.info.Radius1;
                         FileName = append(obj.info.FilenameRaw, "1");
                         CurrentMovie.LoadTraces(FileName);
                         CurrentMovie.TracesAnalysis(Radius, 1);
                         CurrentMovie.PlotTrends;
                         if ~isnan(obj.info.CutTraces)
                             CurrentMovie.PlotDistributions(1);
                         end
                     end
                 catch
                     disp("fail");
                 end
                 close all
             end


             Results1 = array2table(Results1, 'VariableNames', {'Base','Height','slope1', 'slope2','Inflection point'});
             FileNameSave = append(obj.raw.FilePath, filesep, 'FitResultsCh1.mat');
             save(FileNameSave, "Results1")

             Results2 = array2table(Results2, 'VariableNames', {'Base','Height','slope1', 'slope2','Inflection point'});
             FileNameSave = append(obj.raw.FilePath, filesep, 'FitResultsCh2.mat');
             save(FileNameSave, "Results2")
        end
    end
end
