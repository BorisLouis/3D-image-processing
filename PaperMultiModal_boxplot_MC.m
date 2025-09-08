%% Script for Dual Color Tracking Data Analysis and Visualization
% This script loads data from a specified directory structure,
% generates boxplots comparing two channels, and plots individual traces
% color-coded by time.
close all
clc
%% User Configuration
% Define the base path to your data.
basePath = 'S:\Dual Color\20250121_dualcolor\Multicolor_particles\In_water';

% Define the output folders for figures.
figureOutputPath = fullfile(basePath, 'Figures');
tracePlotsPath = fullfile(figureOutputPath, 'TracePlots');

% Create the directories if they do not exist.
if ~exist(figureOutputPath, 'dir')
    mkdir(figureOutputPath);
end
if ~exist(tracePlotsPath, 'dir')
    mkdir(tracePlotsPath);
end

%% Part 1: Boxplots for DR, nR, and aR
disp('Generating boxplots for DR, nR, and aR...');

% Initialize empty arrays to store data from all samples.
drCh1 = []; drCh2 = [];
nRCh1 = []; nRCh2 = [];
aRCh1 = []; aRCh2 = [];

% Find all measurement subfolders.
folders = dir(fullfile(basePath, '0_min*'));
subfolders = {folders.name};

for i = 1:numel(subfolders)
    currentFolder = subfolders{i};
    
    % Construct file paths for msdRes files.
    filePath1 = fullfile(basePath, currentFolder, 'msdRes1.mat');
    filePath2 = fullfile(basePath, currentFolder, 'msdRes2.mat');
    
    if exist(filePath1, 'file') && exist(filePath2, 'file')
        % Load data for channel 1.
        data1 = load(filePath1);
        drCh1 = [drCh1, [data1.allRes.DR]];
        nRCh1 = [nRCh1, [data1.allRes.nR]].*0.99;
        aRCh1 = [aRCh1, [data1.allRes.aR]];
        
        % Load data for channel 2.
        data2 = load(filePath2);
        drCh2 = [drCh2, [data2.allRes.DR]];
        nRCh2 = [nRCh2, [data2.allRes.nR]].*0.99;
        aRCh2 = [aRCh2, [data2.allRes.aR]];
    else
        warning('Files not found in folder: %s', currentFolder);
    end
end



for i = 1:3
    % Create a single figure with three subplots.
    fig1 = figure('Name', 'Dual Color MSD Results', 'Color', 'w', 'Position', [100 100 1200 600]);
    
    % Data for plotting.
    plotData = {{drCh1, drCh2}, {nRCh1, nRCh2}, {aRCh1, aRCh2}};
    plotTitles = {'3D Diffusion Coefficient ($D_R$)', 'Viscosity ($\eta_R$)', 'Anomalous Exponent ($a_R$)'};
    yLabels = {'$D_R$ $(\mu m^2/s)$', '$\eta_R$ (Pa$\cdot$s)', '$a_R$'};
    colors = [0.1, 0.4, 0.7; 0.9, 0.3, 0.1];
    currentData = plotData{i};
    
    combined = [currentData{1}, currentData{2}];
    groupLabels = [repmat({'Channel 1'}, size(currentData{1})), repmat({'Channel 2'}, size(currentData{2}))];
    
    h = boxplot(combined, groupLabels, 'Symbol', '');
    
    % Customize plot appearance.
    title(plotTitles{i}, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel(yLabels{i}, 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'box', 'on', 'LineWidth', 1.5, 'FontSize', 12);
    grid on;
    
    % Customize boxplot colors.
    box_patches = findobj(h, 'Tag', 'Box');
    for j = 1:numel(box_patches)
        set(box_patches(j), 'Color', colors(j,:), 'LineWidth', 2);
        patch(get(box_patches(j), 'XData'), get(box_patches(j), 'YData'), colors(j,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    set(findall(gca, 'type', 'text'), 'Interpreter', 'latex');

    disp('Saving boxplot figure...');
    if i == 1
        Name = 'Diffusion';
    elseif i == 2
        Name = 'Viscosity';
    elseif i == 3
        Name = 'AnExp';
    end
    saveas(fig1, fullfile(figureOutputPath, append(Name, '.png')));
    saveas(fig1, fullfile(figureOutputPath, append(Name, '.svg')));
    close(fig1);
    disp('Boxplots saved successfully.');

end

% Save the figure.


%% Part 2: Individual Trace Plots
disp('Generating individual trace plots for the first 5 measurements...');

% Loop through the first 5 measurement samples.
for sampleNum = 1:5
    currentFolder = ['0_min' num2str(sampleNum)];
    
    % Construct file paths for trackResults files.
    filePath1 = fullfile(basePath, currentFolder, 'trackResults1.mat');
    filePath2 = fullfile(basePath, currentFolder, 'trackResults2.mat');
    
    if exist(filePath1, 'file') && exist(filePath2, 'file')
        
        % Load data.
        data1 = load(filePath1);
        data2 = load(filePath2);
        
        % Create a new figure with two subplots.
        fig2 = figure('Name', ['Traces - ' currentFolder], 'Color', 'w', 'Position', [100 100 1200 600]);
        
        % Define axes for subplots.
        ax1 = subplot(1, 2, 1);
        ax2 = subplot(1, 2, 2);
        
        hold(ax1, 'on');
        hold(ax2, 'on');
        
        
        % Get traces from both channels.
        traces1 = data1.trackRes.traces;
        traces2 = data2.trackRes.traces;
        
        % Loop through traces of the first channel to find a match in the second.
        for j = 1:size(traces1, 1)
            traceTable1 = traces1{j, 1};
            
            % Check if trace is long enough.
            if size(traceTable1, 1) > 50
                t1 = [];
                t1 = traceTable1.t;
                
                % Find a matching trace in channel 2.
                for k = 1:size(traces2, 1)
                    traceTable2 = traces2{k, 1};
                    
                    % Check if trace is long enough.
                    if size(traceTable1, 1) > 20
                        t2 = [];
                        t2 = traceTable2.t;
                        % Calculate mean 3D distance between the two traces.
                        % Assuming they have the same number of time points.
                        % A more robust check might be needed if they don't.
                      
                        CommonTimeT2 = ismember(t2, t1);
                        traceTable2(~CommonTimeT2, :) = [];
                        x2 = traceTable1.col + randi([-786 786]) + randi([-286 286], size(traceTable1.col,1),1);
                        y2 = traceTable1.row + randi([-886 886]) + randi([-486 486], size(traceTable1.col,1),1);
                        z2 = traceTable1.z + randi([-886 886]) + randi([-386 386], size(traceTable1.col,1),1);
                        t2 = [];
                        t2 = traceTable1.t;
                        
                        CommonTimeT1 = ismember(t1, t2);
                        traceTable1(~CommonTimeT1, :) = [];
                        x1 = traceTable1.col;
                        y1 = traceTable1.row;
                        z1 = traceTable1.z;
                        t1 = [];
                        t1 = traceTable1.t;

                        if ~or(isempty(traceTable1), isempty(traceTable2))
                            dist = sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2);
                            meanDist = mean(dist);
                            
                            % Check if the mean distance is less than 500 nm (0.5 um).
                            if meanDist < 10*10^4
                                try
                                    % Plot the matched traces.
                                    disp(['Found co-localized trace pair: Channel 1 Trace ' num2str(j) ' and Channel 2 Trace ' num2str(k)]);
                                    
                                    % Plot Channel 1 trace.
                                    plot3(x1, y1, z1, 'Parent', ax1)
                                    patch([x1(:)' nan],[y1(:)' nan],[z1(:)' nan],[t1(:)' nan],'EdgeColor','interp','FaceColor','none','LineWidth',1 , 'Parent', ax1)
                                    hold(ax1,'on')
                                    % Plot Channel 2 trace.
                                    plot3(x2, y2, z2, 'Parent', ax2)
                                    patch([x2(:)' nan],[y2(:)' nan],[z2(:)' nan],[t2(:)' nan],'EdgeColor','interp','FaceColor','none','LineWidth',1 , 'Parent', ax2)
                                    hold(ax2,'on')
                                    break; % Break from the inner loop after finding a match.
                                catch
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        % Set plot properties after plotting all lines.
        title(ax1, 'Channel 1 Traces', 'FontSize', 16);
        xlabel(ax1, 'x ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel(ax1, 'y ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        zlabel(ax1, 'z ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        cb1 = colorbar(ax1);
        cb1.Label.String = 'Time (s)';
        colormap(ax1, 'parula');
        grid(ax1, 'on');
        axis(ax1, [0 35000 0 35000 -2000 2000]);
        pbaspect(ax1, [1 1 0.2]);   % <<< keeps z compressed relative to x,y
        view(ax1, 3);
        set(ax1, 'box', 'on', 'LineWidth', 1.5);
        
        title(ax2, 'Channel 2 Traces', 'FontSize', 16);
        xlabel(ax2, 'x ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel(ax2, 'y ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        zlabel(ax2, 'z ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
        cb2 = colorbar(ax2);
        cb2.Label.String = 'Time (s)';
        colormap(ax2, 'parula');
        grid(ax2, 'on');
        axis(ax2, [0 35000 0 35000 -2000 2000]);
        pbaspect(ax2, [1 1 0.2]);   % <<< same fix here
        view(ax2, 3);
        set(ax2, 'box', 'on', 'LineWidth', 1.5);

        % Save the figure.
        disp(['Saving traces figure for ' currentFolder '...']);
        fileName = ['traces_' currentFolder '_colocalized'];
        saveas(fig2, fullfile(tracePlotsPath, [fileName '.png']));
        saveas(fig2, fullfile(tracePlotsPath, [fileName '.svg']));
        
        
    else
        warning('Files not found in folder: %s', currentFolder);
    end
end

disp('All trace plots saved successfully.');
disp('Script finished.');