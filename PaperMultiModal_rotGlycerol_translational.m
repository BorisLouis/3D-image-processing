% Load data
load('S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\translational analysis\translational_results.mat');

% Define glycerol contents and ground truth viscosities
glycerol_contents = [80, 85, 90, 95];
ground_truth = [75, 130, 241, 485];

% Extract viscosity data
eta_80 = results.glycerol_80.eta;
eta_85 = results.glycerol_85.eta;
eta_90 = results.glycerol_90.eta;
eta_95 = results.glycerol_95.eta;

% Combine for boxplot
all_eta = [eta_80; eta_85; eta_90; eta_95];
group = [repmat(80, length(eta_80), 1);
         repmat(85, length(eta_85), 1);
         repmat(90, length(eta_90), 1);
         repmat(95, length(eta_95), 1)];

% Create boxplot
figure;
h = boxplot(all_eta, group, 'Colors', 'k', 'Symbol', '', 'Widths', 0.6);
xlabel('Glycerol content (%)');
ylabel('Viscosity (mPa·s)');
title('Viscosity vs. Glycerol Content');
grid on;
set(gca, 'FontSize', 12);
ylim([0 500])
box on;
hold on;

% --- Fix: use x positions corresponding to boxplot groups ---
% boxplot places boxes at positions 1:N, not at 80,85,90,95
numGroups = numel(glycerol_contents);
xPositions = 1:numGroups;

% Overlay dashed line for ground truth viscosities
plot(xPositions, ground_truth, 'k--', 'LineWidth', 1.5);
plot(xPositions, ground_truth, 'ko', 'MarkerFaceColor', 'w');

% Fix x-axis tick labels to show actual glycerol content
set(gca, 'XTick', xPositions, 'XTickLabel', string(glycerol_contents));

