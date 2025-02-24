Table = readmatrix('S:\Dual Color\NewPartSizesResultsPAA_no_negatives.xlsx', 'Sheet', 'Viscosity');
Table = Table(3:end, :);

PS200 = Table(:, 1:9);
PS100 = Table(:, 11:19);
PS300 = Table(:, 22:30);

timepoints = [0 3 5 7 9 11 13 15 17];
numTimepoints = length(timepoints);

groupedData = [];
groupLabels = [];
timeLabels = [];

for i = 1:numTimepoints
    groupedData = [groupedData; PS100(:, i); PS200(:, i); PS300(:, i)];
    groupLabels = [groupLabels; repmat({'100'}, size(PS100,1), 1); ...
                              repmat({'200'}, size(PS200,1), 1); ...
                              repmat({'300'}, size(PS300,1), 1)];
    timeLabels = [timeLabels; repmat(timepoints(i), size(PS100,1) + size(PS200,1) + size(PS300,1), 1)];
end

figure;
boxplot(groupedData, {timeLabels, groupLabels}, 'FactorGap', [5 1], 'LabelVerbosity', 'minor', 'Symbol', '');

xlabel('Tijdstip');
ylabel('Viscosity (cP)');
ylim([0 10])
title('Gegroepeerde Boxplot van PS100, PS200 en PS300');
xtickangle(-45)


hold on;
positions = unique(get(gca, 'XTick')); % Posities van de boxplots op de x-as
means = []; % Opslaan van gemiddelde waarden

for i = 1:length(positions)
    idx = find(timeLabels == timepoints(floor((i-1)/3) + 1)); 
    if mod(i,3) == 1
        group_idx = strcmp(groupLabels(idx), '100');
    elseif mod(i,3) == 2
        group_idx = strcmp(groupLabels(idx), '200');
    else
        group_idx = strcmp(groupLabels(idx), '300');
    end
    mean_val = mean(groupedData(idx(group_idx)));
    means = [means; mean_val];
    plot(positions(i), mean_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

hold off;