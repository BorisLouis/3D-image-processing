FolderName = 'D:\Multimodal tracking\20250708\AllSamples';
Folder = dir(FolderName);

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
            if strcmp(File(end-3:end), 'e0_1')
                Diff0 = [Diff0; [allRes.DR]'];
                a0 = [a0; [allRes.aR]'];
                n0 = [n0; [allRes.nR]'];
                phase0 = [phase0; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), 'e2_1')
                Diff2 = [Diff2; [allRes.DR]'];
                a2 = [a2; [allRes.aR]'];
                n2 = [n2; [allRes.nR]'];
                phase2 = [phase2; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), 'e4_1')
                Diff4 = [Diff4; [allRes.DR]'];
                a4 = [a4; [allRes.aR]'];
                n4 = [n4; [allRes.nR]'];
                phase4 = [phase4; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), 'e6_1')
                Diff6 = [Diff6; [allRes.DR]'];
                a6 = [a6; [allRes.aR]'];
                n6 = [n6; [allRes.nR]'];
                phase6 = [phase6; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), 'e8_1')
                Diff8 = [Diff8; [allRes.DR]'];
                a8 = [a8; [allRes.aR]'];
                n8 = [n8; [allRes.nR]'];
                phase8 = [phase8; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), '10_1')
                Diff10 = [Diff10; [allRes.DR]'];
                a10 = [a10; [allRes.aR]'];
                n10 = [n10; [allRes.nR]'];
                phase10 = [phase10; [allRes.Phase]'];
            elseif strcmp(File(end-3:end), '13_1')
                Diff13 = [Diff13; [allRes.DR]'];
                a13 = [a13; [allRes.aR]'];
                n13 = [n13; [allRes.nR]'];
                phase13 = [phase13; [allRes.Phase]'];
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
varNames = {'Diff', 'a', 'n', 'phase'};
yLabels = {'Diffusion', 'Anomalous exponent', 'Viscosity', 'Phase'};

for v = 1:length(varNames)
    figure;
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
    if v == 4
        ylim([-0.01 0.01])
    end
    grid on;
    hold off;
end



allDiff = [];
allPhase = [];

for i = 1:length(timepoints)
    d = eval(['Diff' num2str(timepoints(i))]);
    p = eval(['phase' num2str(timepoints(i))]);
    
    allDiff = [allDiff; d];
    allPhase = [allPhase; p];
end

figure;
scatter(allPhase, allDiff, 20, 'filled');
xlabel('Phase');
ylabel('Diffusion');
title('Diffusion vs. Phase (All Timepoints Combined)');
grid on;



figure;
hold on;

colors = lines(length(timepoints));
legendEntries = strings(1, length(timepoints));

for i = 1:length(timepoints)
    d = eval(['Diff' num2str(timepoints(i))]);
    p = eval(['phase' num2str(timepoints(i))]);
    
    scatter(p, d, 25, 'filled', 'MarkerFaceColor', colors(i,:));
    legendEntries(i) = ['Time ' num2str(timepoints(i)) ' min'];
end

xlabel('Phase');
ylabel('Diffusion');
title('Diffusion vs. Phase (by Timepoint)');
legend(legendEntries);
grid on;
hold off;