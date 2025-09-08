%% Script for Rotational Tracking Data Analysis and Visualization
% This script loads data from a specified directory structure,
% generates a boxplot of rotational speeds, and plots individual trace
% intensities.

%% User Configuration
% Define the base path to your data.
basePath = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational';

% Define the output folder for figures.
figureOutputPath = fullfile(basePath, 'Figures');
tracePlotsPath = fullfile(figureOutputPath, 'TracePlots');

% Create the directories if they do not exist.
if ~exist(figureOutputPath, 'dir')
    mkdir(figureOutputPath);
end
if ~exist(tracePlotsPath, 'dir')
    mkdir(tracePlotsPath);
end

%% Part 1: Boxplot of Rotational Speed (vTheta) with Interrupted Y-Axis
disp('Generating vTheta boxplot with interrupted y-axis...');

% Load the 10 ms data.
data10ms = load(fullfile(basePath, '10ms_exp.mat'));
vTheta10ms = [data10ms.AllMovieResults.vTheta];

% Load the 100 ms data.
data100ms = load(fullfile(basePath, '100ms_exp.mat'));
vTheta100ms = [data100ms.AllMovieResults.vTheta];

% Combine data for plotting.
combinedData = [vTheta10ms, vTheta100ms];
groupLabels = [repmat({'10 ms'}, size(vTheta10ms)), repmat({'100 ms'}, size(vTheta100ms))];

% Create a new figure.
fig1 = figure('Name', 'vTheta Boxplot', 'Color', 'w');

% Define subplot positions to stack them and remove space.
ax1 = axes('Position', [0.13 0.53 0.775 0.39]); % Top part of the plot
ax2 = axes('Position', [0.13 0.11 0.775 0.39]); % Bottom part of the plot

% Define colors for the boxplots.
colors = [0.1, 0.4, 0.7; 0.9, 0.3, 0.1];

% --- Top Plot (Y-axis from 20 to 35) ---
hold(ax1, 'on');
h1 = boxplot(ax1, combinedData, groupLabels, 'Labels', {'10 ms Exposure', '100 ms Exposure'});
set(h1, 'LineWidth', 2);
title(ax1, 'Rotational Speed ($v_{\theta}$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel(ax1, 'Rotational Speed (rad/s)', 'FontSize', 16);
set(ax1, 'YLim', [20, 35], 'XTick', [], 'box', 'on', 'LineWidth', 1.5, 'FontSize', 14);
set(ax1, 'YTick', 20:5:35);

% Customize top boxplot colors.
box_patches1 = findobj(h1, 'Tag', 'Box');
for i = 1:numel(box_patches1)
    set(box_patches1(i), 'Color', colors(i,:), 'LineWidth', 2);
    patch(get(box_patches1(i), 'XData'), get(box_patches1(i), 'YData'), colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% --- Bottom Plot (Y-axis from 0 to 10) ---
hold(ax2, 'on');
h2 = boxplot(ax2, combinedData, groupLabels, 'Labels', {'10 ms Exposure', '100 ms Exposure'});
set(h2, 'LineWidth', 2);
xlabel(ax2, 'Exposure Time', 'FontSize', 16);
set(ax2, 'YLim', [4, 6], 'box', 'on', 'LineWidth', 1.5, 'FontSize', 14);
set(ax2, 'YTick', 0:2:10);

% Customize bottom boxplot colors.
box_patches2 = findobj(h2, 'Tag', 'Box');
for i = 1:numel(box_patches2)
    set(box_patches2(i), 'Color', colors(i,:), 'LineWidth', 2);
    patch(get(box_patches2(i), 'XData'), get(box_patches2(i), 'YData'), colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add break lines to show the axis discontinuity.
break_line_top = plot(ax1, ax1.XLim, [20 20], '--k', 'LineWidth', 1.5);
break_line_bottom = plot(ax2, ax2.XLim, [10 10], '--k', 'LineWidth', 1.5);

% Clean up the labels and ticks to create the illusion of a single plot.
set(ax1, 'XTickLabel', []);
set(ax2, 'XTickLabel', {'10 ms Exposure', '100 ms Exposure'});
grid(ax1, 'on');
grid(ax2, 'on');
set(findall(fig1, 'type', 'text'), 'Interpreter', 'latex');

% Save the figure in high resolution.
disp('Saving vTheta boxplot with interrupted y-axis...');
saveas(fig1, fullfile(figureOutputPath, 'vTheta_boxplot_interrupted.png'));
saveas(fig1, fullfile(figureOutputPath, 'vTheta_boxplot_interrupted.svg'));
close(fig1);
disp('vTheta boxplot saved successfully.');

%% Part 2: Individual Trace Plots
disp('Generating individual trace plots...');

% Exposure times to loop through.
expTimes = {'10ms_exp', '100ms_exp'};

for i = 1:numel(expTimes)
    currentExp = expTimes{i};
    disp(['Processing traces for ' currentExp '...']);
    
    % Loop through sample folders from 1 to 10.
    for sampleNum = 1:10
        sampleFolderName = ['sample_' num2str(sampleNum)];
        
        % Construct the path to the Traces3DCommon.mat file.
        filePath = fullfile(basePath, currentExp, sampleFolderName, 'Traces3DCommon.mat');
        
        if exist(filePath, 'file')
            
            % Load the file.
            data = load(filePath);
            
            % Loop through each trace in the table.
            for traceNum = 1:size(data.CommonTraces, 1)
                
                % Extract the data for the current trace.
                time = data.CommonTraces.Time{traceNum};
                int1 = data.CommonTraces.Int1{traceNum};
                int2 = data.CommonTraces.Int2{traceNum};
                totInt = data.CommonTraces.TotInt{traceNum};
                
                % Create a new figure for each trace.
                fig2 = figure('Name', ['Trace ' num2str(traceNum) ' - ' currentExp], 'Color', 'w');
                ax = gca;
                hold on;
                
                % Plot the three intensity curves with different colors and line styles.
                plot(ax, time, int1, 'Color', [0.1 0.4 0.7], 'LineWidth', 2, 'DisplayName', 'Channel 1');
                plot(ax, time, int2, 'Color', [0.9 0.3 0.1], 'LineWidth', 2, 'DisplayName', 'Channel 2');
                plot(ax, time, totInt, 'Color', [0.2 0.8 0.2], 'LineWidth', 2, 'DisplayName', 'Total Intensity');
                
                hold off;
                
                % Pimp the plot for scientific publication.
                title(['Trace ' num2str(traceNum) ' - ' strrep(currentExp, '_', ' ')], 'FontSize', 18);
                xlabel('Time (s)', 'FontSize', 14);
                ylabel('Intensity (A.U.)', 'FontSize', 14);
                legend('show', 'Location', 'best');
                ax.FontSize = 12;
                ax.LineWidth = 1.5;
                box on;
                
                % Save the figure.
                fileName = ['trace_' num2str(traceNum) '_' currentExp];
                disp(['Saving ' fileName '...']);
                saveas(fig2, fullfile(tracePlotsPath, [fileName '.png']));
                saveas(fig2, fullfile(tracePlotsPath, [fileName '.svg']));
                
                close(fig2);
            end
        else
            warning('File not found: %s', filePath);
        end
    end
end

disp('All trace plots saved successfully.');
disp('Script finished.');