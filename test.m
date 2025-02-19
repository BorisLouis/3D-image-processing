clc
clear all

DiffusionMatrix = readmatrix('D:\Dual Color\NewPartSizesResultsPAA_no_negatives.xlsx', 'Sheet', 'Viscosity');
DiffusionMatrix = DiffusionMatrix(4:end, :);

% Read matrix
DataMatrix = DiffusionMatrix;

% Define time points
time_points = [0, 3, 5, 7, 9, 11, 13, 15, 17];
num_timepoints = length(time_points);
num_samples = 9;  % Each set has 9 samples

% Extract data for each set
PS200 = DataMatrix(:, 1:9);   % Columns 1 to 9
PS100 = DataMatrix(:, 11:19); % Columns 11 to 19
PS300 = DataMatrix(:, 22:30); % Columns 22 to 30

% Reshape data into column vectors for boxplot
for i = 1:size(PS100, 2)
    data(:,i*3-2) = PS100(:, i);
    data(:,i*3-1) = PS200(:, i);
    data(:,i*3) = PS300(:, i);
end
%data = [PS100, PS200, PS300];  

% Define x positions for boxplots
offset = 0.2; % Adjust spacing between boxplots within the same timepoint
x_positions = reshape(repmat(time_points, 3, 1) + [-offset;0; offset], 1, []);

% Define colors
color_PS100 = [119, 172, 48] / 255;  % Soft green
color_PS200 = [217, 83, 25] / 255;   % Soft red
color_PS300 = [0, 114, 189] / 255;   % Soft blue
box_colors = [color_PS100; color_PS200; color_PS300]; 

% Compute mean values for each timepoint

mean_PS100 = nanmean(PS100, 1);
mean_PS200 = nanmean(PS200, 1);
mean_PS300 = nanmean(PS300, 1);
% Create figure
figure(); % Large figure
hold on;

% % Create boxplot
boxplot(data, 'Positions', x_positions, 'Colors', 'k', 'Widths', 0.2, 'Symbol', '*');  

% Fill boxplots with color
h = findobj(gca, 'Tag', 'Box');
for i = 1:length(h)
    if mod(i,3) == 2
        c = 1;
    elseif mod(i,3) == 1
        c = 2;
    elseif mod(i,3) == 0
        c = 3;
    end
    patch(get(h(i), 'XData'), get(h(i), 'YData'), box_colors(c, :), ...
        'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.2);
end

% Connect mean values with lines
plot(time_points, mean_PS100, '-o', 'Color', color_PS100, 'MarkerFaceColor', color_PS100, 'LineWidth', 2);
plot(time_points, mean_PS200, '-o', 'Color', color_PS200, 'MarkerFaceColor', color_PS200, 'LineWidth', 2);
plot(time_points, mean_PS300, '-o', 'Color', color_PS300, 'MarkerFaceColor', color_PS300, 'LineWidth', 2);

% Customize plot
xlabel('Time (min)');
ylabel('Viscosity (cP)');
ylim([0 1500])
title('Viscosity: outliers');
xticks(time_points);
legend({'PS100', 'PS200', 'PS300'}, 'Location', 'Best');
set(gca, 'FontSize', 24);
grid off;
hold off;