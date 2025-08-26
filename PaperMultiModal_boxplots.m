clc
close all

MainFolder = 'S:\Dual Color\20250121_dualcolor\PS_286g_136r\Gelation\17_min\Cutted';
Folder = dir(MainFolder);
Folder([Folder.isdir] ~= 1) = [];

OutputFolder = fullfile(MainFolder, 'Figures');
if ~exist(OutputFolder, 'dir')
    mkdir(OutputFolder);
end

colors = lines(20);
jitterWidth = 0.3;
timePoints = [17.000 17.016, 17.032, 17.050, 17.062, 17.078, 17.083, 17.100]; % x-axis labels

% Preallocate storage for boxplots
AllDiff = cell(2, numel(timePoints));
AllAnExp = cell(2, numel(timePoints));
AllVisc = cell(2, numel(timePoints));

AllxStd = cell(2, numel(timePoints));
AllyStd = cell(2, numel(timePoints));
AllzStd = cell(2, numel(timePoints));
AllrStd = cell(2, numel(timePoints));

% Create figures up front for scatter + legend
legendHandles = [];
legendLabels  = {};

for k = 1:2

    for figIdx = 1:7
        figure(100*k + figIdx); hold on;
    end
    
    for i = 3:size(Folder, 1)-2
        Order = [3, 4, 5, 6, 7, 8, 9, 10];
        tpIndex = i-2; % 1...9
        tpMin = timePoints(tpIndex);
        
        TimeFolder = dir(fullfile(Folder(Order(i-2)).folder, Folder(Order(i-2)).name));
        TimeFolder([TimeFolder.isdir] ~= 1) = [];

        TraceLength = [];
        for j = 3:size(TimeFolder, 1)
            try
                SampleFolder = dir(fullfile(TimeFolder(j).folder, TimeFolder(j).name));

                load(fullfile(SampleFolder(1).folder, ['msdRes_', Folder(Order(i-2)).name, '_', num2str(k), '.mat']));
                % load(fullfile(SampleFolder(1).folder, ['trackResults', Folder(Order(i-2)).name, '_', num2str(k), '.mat']));

                %% Process traces
                % allHeight = cellfun(@height, trackRes.traces(:,1));
                % idx = allHeight>50;
                % Traces = trackRes.traces(idx, 1);
                % 
                % xStd = zeros(size(Traces));
                % yStd = zeros(size(Traces));
                % zStd = zeros(size(Traces));
                % rStd = zeros(size(Traces));
                % 
                % TraceLength = [TraceLength; cellfun(@height, Traces(:,1))];
                % for l = 1:size(Traces, 1)
                %     xStd(l,1) = std(Traces{l}.row);
                %     yStd(l,1) = std(Traces{l}.col);
                %     zStd(l,1) = std(Traces{l}.z);
                %     rStd(l,1) = std(sqrt(Traces{l}.row.^2 + Traces{l}.col.^2 + Traces{l}.z.^2));
                % end

                %% Process diffusion/viscosity/anexp
                Diff = [allRes.DR]';
                Visc = [allRes.nR]';
                AnExp = [allRes.aR]';

                if i == 6
                    Diff = Diff.*3;
                    Visc = Visc./3;
                end

                % Append to global storage for boxplots
                % AllxStd{k,tpIndex} = [AllxStd{k,tpIndex}; xStd];
                % AllyStd{k,tpIndex} = [AllyStd{k,tpIndex}; yStd];
                % AllzStd{k,tpIndex} = [AllzStd{k,tpIndex}; zStd];
                % AllrStd{k,tpIndex} = [AllrStd{k,tpIndex}; rStd];

                AllDiff{k,tpIndex} = [AllDiff{k,tpIndex}; Diff];
                AllAnExp{k,tpIndex} = [AllAnExp{k,tpIndex}; AnExp];
                AllVisc{k,tpIndex} = [AllVisc{k,tpIndex}; Visc];

                % Pick a color per sample j
                col = colors(mod(j-3, size(colors,1))+1, :);

                % Random jittered x-positions
                xjit = @(N) tpMin + (rand(N,1)-0.5)*jitterWidth;

                % Scatter plots
                % figure(100*k+1); scatter(xjit(numel(xStd)), xStd, 30, col);
                % figure(100*k+2); scatter(xjit(numel(yStd)), yStd, 30, col);
                % figure(100*k+3); scatter(xjit(numel(zStd)), zStd, 30, col);
                % figure(100*k+4); scatter(xjit(numel(rStd)), rStd, 30, col);
                figure(100*k+5); scatter(xjit(numel(Diff)), Diff, 30, col);
                figure(100*k+6); scatter(xjit(numel(AnExp)), AnExp, 30, col);
                figure(100*k+7); scatter(xjit(numel(Visc)), Visc, 30, col);

                if k==1 && i==3   % first time we encounter each sample
                    legendHandles(end+1) = scatter(nan, nan, 30, col, 'filled');
                    legendLabels{end+1} = sprintf('Sample %d', j-2);
                end

            catch
            end
        end
        % AvTraceLength(i,k) = mean(TraceLength);
        % FullTraceLength{i,k} = TraceLength;

    end
end

% ✅ Add legends to scatter figures
for k = 1:2
    for figIdx = 1:7
        figure(100*k+figIdx);
        legend(legendHandles, legendLabels, 'Location', 'bestoutside');
        xlabel('Polymerization time (min)');
    end
end

%% ✅ Boxplots with trendlines
boxColors = [0.2 0.6 0.2; 0.8 0.2 0.2]; % k=1 green, k=2 red
fillColors = [0.7 0.9 0.7; 0.95 0.7 0.7]; % softer fill colors
labels = string(timePoints);

function plotBox(allData, ttl, ylbl, logY, timePoints, labels, boxColors, fillColors)
    figure; hold on; title(ttl, 'FontSize', 14, 'FontWeight', 'bold');
    offset = [-0.002, +0.002]; % shift boxes for k=1 and k=2 so they don't overlap

    for k=1:2
        data = allData(k,:);
        grp = [];
        vals = [];
        pos  = [];
        for t=1:numel(data)
            vals = [vals; data{t}];
            grp  = [grp; t*ones(size(data{t}))];
            pos  = [pos; (timePoints(t) + offset(k)) * ones(size(data{t}))];
        end
        if isempty(vals), continue; end

        % Make grouped boxplot
        b = boxplot(vals, pos, 'Colors', boxColors(k,:), 'Symbol','', ...
                    'Widths', 0.004, 'Positions', pos, 'Labels', []);
        set(b,{'LineWidth'},{1.2});

        % Fill boxes with softer colors
        h = findobj(gca,'Tag','Box');
        if k == 1
            for j=1:length(h)
                patch(get(h(j),'XData'), get(h(j),'YData'), fillColors(k,:), ...
                      'FaceAlpha',0.6, 'EdgeColor',boxColors(k,:), 'LineWidth',1.2);
                maxbox = size(h, 1);
            end
        elseif k == 2
            for j=1:maxbox
                patch(get(h(j),'XData'), get(h(j),'YData'), fillColors(k,:), ...
                      'FaceAlpha',0.6, 'EdgeColor',boxColors(k,:), 'LineWidth',1.2);
            end
        end

        % Means + SEM error bars + dashed trendline
        means = cellfun(@mean,data);
        sems  = cellfun(@(d) std(d)./sqrt(numel(d)), data);
        errorbar(timePoints+offset(k), means, sems, 'o--', ...
                 'Color', boxColors(k,:), 'LineWidth',1.5, ...
                 'MarkerFaceColor', boxColors(k,:), 'MarkerSize',6);
    end

    % Axis & style
    set(gca,'XTick',timePoints,'XTickLabel',labels,'FontSize',12);
    xlabel('Polymerization time (min)'); ylabel(ylbl);
    if logY, set(gca,'YScale','log'); end
    grid on; box on;
    legend({'PS 300 nm','', '', '', '', '', '', '', '', '', '', '', '', 'PS 100 nm'}, 'Location','best');
    xlim([16.999 17.101])
    ylim([0 150])
end

% Call plotting function
plotBox(AllDiff, 'Diffusion vs Time', 'Diffusion coefficient (µm^2/s)', false, timePoints, labels, boxColors, fillColors);
plotBox(AllAnExp, 'Anomalous exponent vs Time', 'Anomalous exponent', false, timePoints, labels, boxColors, fillColors);
plotBox(AllVisc, 'Viscosity vs Time', 'Viscosity (cP)', false, timePoints, labels, boxColors, fillColors);
plotBox(AllxStd, 'xStd spread vs Time', 'xStd', false, timePoints, labels, boxColors, fillColors);
plotBox(AllyStd, 'yStd spread vs Time', 'yStd', false, timePoints, labels, boxColors, fillColors);
plotBox(AllzStd, 'zStd spread vs Time', 'zStd', false, timePoints, labels, boxColors, fillColors);
plotBox(AllrStd, 'rStd spread vs Time', 'rStd', false, timePoints, labels, boxColors, fillColors);

%% ✅ Save all figures
figHandles = findall(0,'Type','figure');
for f = 1:numel(figHandles)
    fig = figHandles(f);
    name = get(get(fig,'CurrentAxes'),'Title');
    if (name == ""), name = sprintf('Figure%d', f); else, name = name.String; end
    safeName = regexprep(name,'\W','_');
    saveas(fig, fullfile(OutputFolder, [safeName '.png']));
    saveas(fig, fullfile(OutputFolder, [safeName '.svg']));
end