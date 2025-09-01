FolderName = 'D:\Multimodal tracking\20250724\alldata';
OutputFolder = 'D:\Multimodal tracking\20250724\Analysis\Gradient';
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
    ylabel(yLabels{v});
    % if v == 4
    %     ylim([-0.01 0.01])
    % end
    grid on;
    hold off;
    
    saveas(Fig, append(OutputFolder, filesep, varNames{v}, '.png'))

end



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
title('Diffusion vs. Phase/Intensity (All Timepoints Combined)');
grid on;
saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhase.png'))



Fig = figure;
hold on;

colors = lines(length(timepoints));
legendEntries = strings(1, length(timepoints));

for i = 1:length(timepoints)
    d = eval(['Diff' num2str(timepoints(i))]);
    p = eval(['phase' num2str(timepoints(i))]);
    
    scatter(p, d, 25, 'filled', 'MarkerFaceColor', colors(i,:));
    legendEntries(i) = ['Time ' num2str(timepoints(i)) ' min'];
end

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
title('Diffusion vs. Phase/Intensity (by Timepoint)');
legend(legendEntries);
grid on;
saveas(Fig, append(OutputFolder, filesep, 'DiffvsPhasePerTimepoint.png'))
hold off;