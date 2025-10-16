% Define paths
rotPath = 'S:\Rotational Tracking\20250708_AuBPs_184x92_PAA\data\';
transPath = 'S:\Rotational Tracking\20250708_AuBPs_184x92_PAA\Translational analysis\translational_results.mat';

% Define filenames (rotational data)
rotFiles = {'1s_3min', '2s_5min', '3s_7min', '4s_9min', '5s_11min', ...
            '6s_13min', '7s_15min', '8s_17min', '9s_19min'};

% Define timepoints for x-axis
timepoints = [3 5 7 9 11 13 15 17 19]; % minutes

% Load translational data
load(transPath, 'results');

% Preallocate cell arrays
eta_rot = cell(1, numel(rotFiles));
eta_trans = cell(1, numel(rotFiles));

% Loop through timepoints
for i = 1:numel(rotFiles)
    % --- Rotational ---
    data = load(fullfile(rotPath, rotFiles{i}));
    allRes = data.AllMovieResults;
    eta_rot{i} = [allRes.nr]';  % Force column vector

    % --- Translational ---
    fieldName = ['t_' rotFiles{i}];
    eta_trans{i} = results.(fieldName).eta(:);  % Force column vector
end

%% Prepare figure
figure('Color', 'w'); hold on;
set(gca, 'YScale', 'log');

% Box positions
positions_trans = timepoints - 0.15;
positions_rot = timepoints + 0.15;

boxplotData = [];
group = [];
positions = [];

for i = 1:numel(timepoints)
    % Add translational (green)
    boxplotData = [boxplotData; eta_trans{i}];
    group = [group; repmat(i*2-1, numel(eta_trans{i}), 1)];

    % Add rotational (yellow)
    boxplotData = [boxplotData; eta_rot{i}];
    group = [group; repmat(i*2, numel(eta_rot{i}), 1)];
end

positions = reshape([positions_trans; positions_rot], 1, []); % interleave

% Create boxplots
h = boxplot(boxplotData, group, 'Positions', positions, 'Colors', 'k', ...
    'Symbol', '', 'Widths', 0.3, 'Whisker', 1.5);

% Color the boxes (green = translational, yellow = rotational)
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    if mod(j,2)==0
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), [1 1 0], 'FaceAlpha', 0.6, 'EdgeColor', 'k'); % yellow
    else
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), [0 1 0], 'FaceAlpha', 0.6, 'EdgeColor', 'k'); % green
    end
end

% Axis formatting
xlim([min(timepoints)-1, max(timepoints)+1]);
ylim([1e-1, 1e3]);
xticks(timepoints);
xlabel('Polymerization time (min)');
ylabel('Viscosity (cP)');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on; box on;

%% Secondary axis
yyaxis right
secondary_vals = [1.388e25, 1.055e25, 1.5066e25, 3.42e23, 5.28e23, ...
                  6.918e23, 3.616e24, 5.371e24, 1.168e25];
plot(timepoints, secondary_vals, '--ok', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6);
ylabel('Secondary quantity (a.u.)');
set(gca, 'YScale', 'log');

yyaxis left
title('Viscosity evolution during polymerization');

legend({'Translational (green)', 'Rotational (yellow)', 'Secondary data'}, ...
       'Location', 'northwest');
