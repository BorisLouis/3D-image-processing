clc 
clear 
close all;
%file info

Paths = {'D:\Data Hannah\20260203'};
Paths2DCal = {};
TimePaths = {'200nm_in_PAA\brightfield', '300nm_in_PAA\brightfield', '500nm_in_PAA\brightfield',...
    '200nm_in_PAA\darkfield', '300nm_in_PAA\darkfield', '500nm_in_PAA\darkfield'};

for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
            MinMatrix = [];
            MaxMatrix = [];
            Folder = dir(append(Paths{j}, filesep, TimePaths{i}));  
            for k = 3:size(Folder, 1)
                try
                file = append(Folder(k).folder, filesep, Folder(k).name, filesep, 'PhaseTrends', filesep, 'ResultsPhaseCalibration.mat');
                load(file);
                MinMatrix(:, end+1) = [tracks.MinRangeMedian]';
                MaxMatrix(:, end+1) = [tracks.MaxRangeMedian]';
                catch
                end
            end
            MinAv = mean(MinMatrix, 2);
            MaxAv = mean(MaxMatrix, 2);
        catch
        end
    end
end
%