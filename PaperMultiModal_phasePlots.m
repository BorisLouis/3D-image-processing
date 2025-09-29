close all
clc

FolderName = 'E:\Multimodal tracking\20250724\alldata';
OutputFolder = 'E:\Multimodal tracking\20250724\Figures\Gradient';
Folder = dir(FolderName);
Parameter = 'GradientMagnitude';

Diff0 = [];
a0 = [];
n0 = [];
phase0 = [];
Diff2 = [];
a2 = [];
n2 = [];
phase2 = [];
Diff4 = [];
a4 = [];
n4 = [];
phase4 = [];
Diff6 = [];
a6 = [];
n6 = [];
phase6 = [];
Diff8 = [];
a8 = [];
n8 = [];
phase8 = [];
Diff10 = [];
a10 = [];
n10 = [];
phase10 = [];
Diff13 = [];
a13 = [];
n13 = [];
phase13 = [];

for i = 3:size(Folder, 1)
    try
        File = append(FolderName, filesep, Folder(i).name);
        if Folder(i).isdir == 1
            load(append(File, filesep, 'msdResPhase.mat'));
            if strcmp(File(end-4:end), 'n0__1')
                GradMag = [allRes.GradientMagnitude];
                allRes = allRes(GradMag < 0.5);
                Diff0 = [Diff0; [allRes.DR]'];
                a0 = [a0; [allRes.aR]'];
                n0 = [n0; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase0 = [phase0; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase0 = [phase0; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase0 = [phase0; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase0 = [phase0; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase0 = [phase0; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), 'n2__1')
                Diff2 = [Diff2; [allRes.DR]'];
                a2 = [a2; [allRes.aR]'];
                n2 = [n2; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase2 = [phase2; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase2 = [phase2; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase2 = [phase2; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase2 = [phase2; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase2 = [phase2; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), 'n4__1')
                allDR = [allRes.DR];
                cutoff = prctile(allDR, 10);
                allRes = allRes(allDR > cutoff);

                Diff4 = [Diff4; [allRes.DR]'];
                a4 = [a4; [allRes.aR]'];
                n4 = [n4; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase4 = [phase4; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase4 = [phase4; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase4 = [phase4; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase4 = [phase4; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase4 = [phase4; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), 'n6__1')
                allDR = [allRes.DR];
                cutoff = prctile(allDR, 20);
                allRes = allRes(allDR > cutoff);

                Diff6 = [Diff6; [allRes.DR]'];
                a6 = [a6; [allRes.aR]'];
                n6 = [n6; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase6 = [phase6; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase6 = [phase6; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase6 = [phase6; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase6 = [phase6; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase6 = [phase6; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), 'n8__1')
                allDR = [allRes.DR];
                cutoff = prctile(allDR, 25);
                allRes = allRes(allDR > cutoff);
                Diff8 = [Diff8; [allRes.DR]'];
                a8 = [a8; [allRes.aR]'];
                n8 = [n8; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase8 = [phase8; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase8 = [phase8; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase8 = [phase8; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase8 = [phase8; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase8 = [phase8; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), '10__1')
                allDR = [allRes.DR];
                cutoff = prctile(allDR, 50);
                allRes = allRes(allDR > cutoff);

                Diff10 = [Diff10; [allRes.DR]'];
                a10 = [a10; [allRes.aR]'];
                n10 = [n10; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase10 = [phase10; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase10 = [phase10; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase10 = [phase10; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase10 = [phase10; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase10 = [phase10; [allRes.SharpnessLaplacian]'];
                end
            elseif strcmp(File(end-4:end), '13__1')
                Diff13 = [Diff13; [allRes.DR]'];
                a13 = [a13; [allRes.aR]'];
                n13 = [n13; [allRes.nR]'];
                if strcmp(Parameter, 'Phase')
                    phase13 = [phase13; [allRes.Phase]'];
                elseif strcmp(Parameter, 'IntPhaseCh')
                    phase13 = [phase13; [allRes.IntPhaseCh]'];
                elseif strcmp(Parameter, 'GradientMagnitude')
                    phase13 = [phase13; [allRes.GradientMagnitude]'];
                elseif strcmp(Parameter, 'LocalVariance')
                    phase13 = [phase13; [allRes.LocalVariance]'];
                elseif strcmp(Parameter, 'SharpnessLaplacian')
                    phase13 = [phase13; [allRes.SharpnessLaplacian]'];
                end
            else
                warning('Filename not found')
            end
        end
    catch
    end

    % load(append(File, '\SegmentMovie\SegmentMovie5'));
    % Frame = Load.Movie.tif.getframes(append(File, '\calibrated1\calibratedPlane1.tif'), 500);
    % Mask = Mask{500,1};
    % Fig = figure();
    % 
    % subplot(1,3,1)
    % imagesc(Frame);
    % colormap("gray")
    % axis image
    % title('Raw data')
    % subplot(1,3,2)
    % imagesc(Mask);
    % axis image
    % title('Mask')
    % subplot(1,3,3)
    % imshowpair(Frame, Mask);
    % axis image
    % title("overlay")
    % 
    % SaveName = append(File, filesep, 'SegmentTestframe500.png');
    % saveas(Fig, SaveName)
end



%% Trends over time

timepoints = [0 2 4 6 8 10 13];
% Define variable base names and labels

if strcmp(Parameter, 'Phase')
    PhaseLabel = 'Phase';
elseif strcmp(Parameter, 'IntPhaseCh')
    PhaseLabel = 'IntPhaseCh';
elseif strcmp(Parameter, 'GradientMagnitude')
    PhaseLabel = 'GradientMagnitude';
elseif strcmp(Parameter, 'LocalVariance')
    PhaseLabel = 'LocalVariance';
elseif strcmp(Parameter, 'SharpnessLaplacian')
    PhaseLabel = 'SharpnessLaplacian';
end
varNames = {'Diff', 'a', 'n', 'phase'};
yLabels = {'Diffusion', 'Anomalous exponent', 'Viscosity', PhaseLabel};

for v = 1:length(varNames)
    Fig = figure;
    hold on;
    means = zeros(size(timepoints));

    for i = 1:length(timepoints)
        % Dynamically access variables
        var = eval([varNames{v} num2str(timepoints(i))]);

        % Boxplot at x = timepoint
        boxchart(ones(size(var)) * timepoints(i), var, 'BoxFaceColor', [0.2 0.6 0.8]);

        % Store mean
        means(i) = mean(var, 'omitnan');
    end

    % Plot mean line
    plot(timepoints, means, '-o', 'LineWidth', 2, 'Color', 'k');

    title([yLabels{v} ' over time']);
    xlabel('Time (minutes)');
    ylim([0.05 0.25])
    ylabel(yLabels{v});
    % if v == 4
    %     ylim([-0.01 0.01])
    % end
    grid on;
    hold off;

    saveas(Fig, append(OutputFolder, filesep, varNames{v}, '.png'))
    saveas(Fig, append(OutputFolder, filesep, varNames{v}, '.svg'))
end


timepoints = [13 10 8  6  4  2  0];
allDiff = [];
allPhase = [];

for i = 1:length(timepoints)
    d = eval(['Diff' num2str(timepoints(i))]);
    p = eval(['phase' num2str(timepoints(i))]);

    allDiff = [allDiff; d];
    allPhase = [allPhase; p];
end

Fig = figure;
scatter(allPhase, allDiff, 20, 'filled');
if strcmp(Parameter, 'Phase')
    xlabel('Phase');
elseif strcmp(Parameter, 'IntPhaseCh')
    xlabel('IntPhaseCh');
elseif strcmp(Parameter, 'GradientMagnitude')
    xlabel('GradientMagnitude');
elseif strcmp(Parameter, 'LocalVariance')
    xlabel('LocalVariance');
elseif strcmp(Parameter, 'SharpnessLaplacian')
    xlabel('SharpnessLaplacian');
end

ylabel('Diffusion');
title(append('Diffusion vs. phase - ', Parameter, ' - (All Timepoints Combined)'));
set(gca, 'XScale', 'log')
xlim([0.010 0.035])
grid on;
saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhase.png'))
saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhase.svg'))


% Fig = figure;
% hold on;
% 
% colors = flipud(lines(length(timepoints)));
% legendEntries = strings(1, length(timepoints));
% 
% for i = 1:length(timepoints)
%     d = eval(['Diff' num2str(timepoints(i))]);
%     p = eval(['phase' num2str(timepoints(i))]);
% 
%     scatter(p, d, 25, 'filled', 'MarkerFaceColor', colors(i,:));
%     legendEntries(i) = ['Time ' num2str(timepoints(i)) ' min'];
% end
% 
% if strcmp(Parameter, 'Phase')
%     xlabel('Phase');
% elseif strcmp(Parameter, 'IntPhaseCh')
%     xlabel('IntPhaseCh');
% elseif strcmp(Parameter, 'GradientMagnitude')
%     xlabel('GradientMagnitude');
% elseif strcmp(Parameter, 'LocalVariance')
%     xlabel('LocalVariance');
% elseif strcmp(Parameter, 'SharpnessLaplacian')
%     xlabel('SharpnessLaplacian');
% end
% ylabel('Diffusion');
% title(append('Diffusion vs. phase - ', Parameter, ' - (by timepoint)'));
% legend(legendEntries);
% set(gca, 'XScale', 'log')
% xlim([0.010 0.055])
% grid on;
% saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhasePerTimepoint.png'))
% saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhasePerTimepoint.svg'))
% hold off;

Fig = figure;
hold on;

% Collect all data in arrays
xAll = []; % parameter values (e.g. LocalVariance)
yAll = []; % diffusion values
gAll = []; % group labels (categorical)

for i = 1:length(timepoints)
    d = eval(['Diff' num2str(timepoints(i))]);
    p = eval(['phase' num2str(timepoints(i))]);   % replace if Parameter ≠ phase
    
    xAll = [xAll; p(:)];
    yAll = [yAll; d(:)];
    gAll = [gAll; repmat(timepoints(i), numel(d), 1)];
end

% Convert group to categorical for labeling
gAll = categorical(gAll, timepoints, "Time " + string(timepoints) + " min");

% Define colors
colors = flipud(lines(length(timepoints)));

% Plot scatter with marginal histograms
h = scatterhist(xAll, yAll, 'Group', gAll, ...
    'Kernel', 'on', ...          % smooth density instead of hist bars
    'Color', colors, ...
    'Marker', '.', ...
    'MarkerSize', 10, ...
    'LineWidth', 1.2);

% Axis labels
if strcmp(Parameter, 'Phase')
    xlabel('Phase');
elseif strcmp(Parameter, 'IntPhaseCh')
    xlabel('IntPhaseCh');
elseif strcmp(Parameter, 'GradientMagnitude')
    xlabel('GradientMagnitude');
elseif strcmp(Parameter, 'LocalVariance')
    xlabel('LocalVariance');
elseif strcmp(Parameter, 'SharpnessLaplacian')
    xlabel('SharpnessLaplacian');
end
ylabel('Diffusion');

% Title
title(['Diffusion vs. phase - ', Parameter, ' - (by timepoint)']);

% Legend
legend('show');

% % Log scale for x-axis
% set(h(1), 'XScale', 'log');
% xlim(h(1), [0.075 0.35]);
% ylim(h(1), [0 2.1]);
% 
% % Also fix the bottom histogram to log scale
% set(h(2), 'XScale', 'log');
% xlim(h(2), [0.075 0.35]);
% 
% ylim(h(3), [0 2.1]);

% Labels and title
xlabel(h(1), Parameter);
ylabel(h(1), 'Diffusion');
title(h(1), ['Diffusion vs. phase - ', Parameter, ' - (by timepoint)']);
grid(h(1), 'on');

% Save
saveas(Fig, fullfile(OutputFolder, 'DiffvsPhasePerTimepoint.png'));
saveas(Fig, fullfile(OutputFolder, 'DiffvsPhasePerTimepoint.svg'));

hold off;