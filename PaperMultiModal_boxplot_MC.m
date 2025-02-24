% Read dataanE
DiffusionMatrix = readmatrix('D:\Dual Color\Results_MC_InWater.xlsx', 'Sheet', 'Diff');
DiffusionMatrix = DiffusionMatrix(2:end, :);
M = DiffusionMatrix;

% Define dimensions
num_measurements = 5;  % 5 measurements per sample
num_channels = 2;      % 2 channels (ch1, ch2)
num_dates = 2;         % 2 measurement dates
gap_size = 1.5; 

% Split data correctly by ROWS
M1 = M(:, 1:10);   % First date (first 10 rows)
M2 = M(:, 12:21);  % Second date (next 10 rows)

% Separate into channels
M1_ch1 = M1(:, 1:num_measurements);       % Rows 1-5 → Channel 1
M1_ch2 = M1(:, num_measurements+1:end);   % Rows 6-10 → Channel 2
M2_ch1 = M2(:, 1:num_measurements);       % Rows 11-15 → Channel 1
M2_ch2 = M2(:, num_measurements+1:end);   % Rows 16-20 → Channel 2

% Convert data to columns for boxplot
data = [M1_ch1, M1_ch2, M2_ch1, M2_ch2];

% === Group labels for x-axis ===
num_samples = num_measurements;  % 5 samples

% X positions: Keep red & green boxplots closer together
x_positions_1 = repmat(1:num_measurements, num_channels, 1) + [-0.2; 0.2];
x_positions_2 = x_positions_1 + num_measurements + gap_size; % Shift second date
x_positions = [x_positions_1(:); x_positions_2(:)]; % Combine both sets

% Custom colors
soft_green = [119, 172, 48] / 255;  % #77AC30
soft_red = [217, 83, 25] / 255;     % #D95319
box_colors = [soft_green; soft_red]; % Assign colors in order

% ========== PLOT FOR FIRST DATE ==========
figure('Position', [100, -100, 1200, 800]);
hold on;
boxplot(data, 'Positions', x_positions, 'Colors', 'k', 'Widths', 0.3, 'Symbol', '');  

% Fill colors
h = findobj(gca, 'Tag', 'Box');
for i = 1:length(h)
    patch(get(h(i), 'XData'), get(h(i), 'YData'), box_colors(mod(i-1, 2) + 1, :), ...
        'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.2);
end

% Customize plot
xlabel('Samples');
ylabel('Diffusion coefficient (µm^2s^-^1');
ylim([0 7])
title('Diffusion Measurements');
xticks([1:num_measurements, num_measurements + gap_size + (1:num_measurements)]);
xticklabels(["Sample1", "Sample2", "Sample3", "Sample4", "Sample5", ...
             "Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]);

set(gca, 'FontSize', 22);
grid off;


text(mean(1:num_measurements), 5.75, '2025 01 21', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(mean(num_measurements + gap_size + (1:num_measurements)), 5.75, '2025 01 22', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
legend({'Channel 1 (Green)', 'Channel 2 (Red)'}, 'Location', 'Best');
hold off;

