% MATLAB script to plot diffusion data from Excel
% File and folder paths
dataFile = 'S:\Dual Color\20250220_stickyness\Test_stickyness_2D\Results_stickyness_particles.xlsx';
saveFolder = 'S:\Dual Color\20250220_stickyness\Test_stickyness_2D\Figures';

% Timepoints (minutes) corresponding to columns E-L
timepoints = [4 6 8 10 12 14 16 18];

% Read data ranges for each sample
sample1 = readmatrix(dataFile,'Range','E4:L25');
sample2 = readmatrix(dataFile,'Range','E26:L71');
sample3 = readmatrix(dataFile,'Range','E72:L94');

% Combine all samples
allData = [sample1; sample2; sample3];

% Reshape into long format for boxplotting
timeVec = [];
diffusionVec = [];
for i = 1:numel(timepoints)
    tp = timepoints(i);
    validValues = allData(:,i);
    validValues = validValues(~isnan(validValues)); % remove NaN
    diffusionVec = [diffusionVec; validValues];
    timeVec = [timeVec; repmat(tp, numel(validValues), 1)];
end

% Create figure
figure('Color','w');
boxplot(diffusionVec, timeVec, 'Symbol','k+');
hold on;

% Calculate and plot mean values
meanVals = grpstats(diffusionVec, timeVec, 'mean');
plot(timepoints, meanVals, '-o','LineWidth',2,'Color','r');

xlabel('Time (minutes)');
ylabel('Diffusion');
title('Diffusion over Time (all samples combined)');
set(gca,'FontSize',12);

% Save outputs
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end
saveas(gcf, fullfile(saveFolder,'Diffusion_over_time.svg'));
saveas(gcf, fullfile(saveFolder,'Diffusion_over_time.png'));