clc 
clear 
close all;
%file info

Paths = {'D:\Polymer Dynamics'};
Paths2DCal = {};
TimePaths = {'20260121\1x_bAA', '20260120\2x_bAA', '20260121\2x_AA', '20260121\3x_AA'};

for j = 1:numel(Paths)
    for i = 1:numel(TimePaths)
        try
           load(append(Paths{j}, filesep, TimePaths{i}, filesep, 'ResCalcMSD1'));
           ResCh1(i, :) = mean(table2array(FinalResults.Fitting(:,1:end-1)),1);

           load(append(Paths{j}, filesep, TimePaths{i}, filesep, 'ResCalcMSD2'));
           ResCh2(i, :) = mean(table2array(FinalResults.Fitting(:,1:end-1)),1);
        catch
        end
    end
end
%