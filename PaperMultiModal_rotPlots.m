close all
FolderName = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol';
OutputFolder = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\fig';

FileNames = {'glycerol_80', 'glycerol_85', 'glycerol_90', 'glycerol_95'};
Labels = {'80', '85', '90', '95'};

allVals = [];     % all viscosity values
groups  = [];     % group labels
MeanViscs = [88.056 177.842 198.656 275.412];
for i = 1:numel(FileNames)
    % Load .mat file
    Folder = dir(append(FolderName, filesep, FileNames{i}));

    Viscosity = [];
    for j = 3:7
        try
            load(append(Folder(j).folder, filesep, Folder(j).name, filesep, 'msadRes.mat'));
            Viscosity = [Viscosity; [allRes.nr]'];
            disp(append('Median Viscosity in 3D is ', num2str(nanmedian([allRes.nr]))));
            
        catch
            disp(append('Failed to load data from ',Folder(j).folder, filesep, Folder(j).name, filesep, 'msadRes.mat' ));
        end
    end
    % CorrFactor(i) = MeanViscs(i)./nanmean(Viscosity);
    % Viscosity = Viscosity.*(MeanViscs(i)./nanmean(Viscosity));
    % Append values
    allVals  = [allVals; Viscosity];                % add values
    groups   = [groups; repmat(i, numel(Viscosity), 1)]; % add group IDs
end

% Make boxplot
% figure;
% boxplot(allVals, groups, 'Labels', Labels, 'Symbol', "");
% xlabel('Glycerol content (v/v %)','FontSize',12,'FontWeight','bold');
% ylabel('Viscosity (cP)','FontSize',12,'FontWeight','bold');
% title('Comparison of Viscosities')
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',11,'LineWidth',1.2); 
% saveas()


Fig = figure('Color','w');  % white background
% Make boxplot
boxplot(allVals, groups, 'Labels', Labels, ...
    'Whisker', 1.5, ...             % standard whisker length
    'Symbol', '', ...              % outlier marker
    'Widths', 0.6, ...              % box width
    'Colors', lines(numel(FileNames)));  % use a colormap for boxes
% Labeling
xlabel('Glycerol content (v/v %)','FontSize',12,'FontWeight','bold');
ylabel('Viscosity (cP)','FontSize',12,'FontWeight','bold');
title('Comparison of Viscosities','FontSize',14);

% Log scale
% set(gca,'YScale','log');
ylim([0 1200])
% Aesthetics
set(gca,'FontSize',11,'LineWidth',1.2);   % thicker axis lines
% grid on;                                  % light grid for readability
box on;

% Make boxes filled (MATLAB doesn’t fill by default)
h = findobj(gca,'Tag','Box');
colors = lines(numel(FileNames));
for j=1:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), colors(j,:), ...
          'FaceAlpha',0.4, 'EdgeColor',colors(j,:));
end
saveas(Fig, append('S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\Figures', filesep, 'ViscGLycerol.png'));
saveas(Fig, append('S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol\Figures', filesep, 'ViscGLycerol.svg'));
% Preallocate arrays
allVals  = [];
groups   = [];
sampleID = [];
figure()
colors = lines(10);
for i = 1:numel(FileNames)
    % Load .mat file
    Folder = dir(append(FolderName, filesep, FileNames{i}));

    Viscosity = [];
    for j = 3:7%idx{i}
        try
            load(append(Folder(j).folder, filesep, Folder(j).name, filesep, 'msadRes.mat'));
            vals = [allRes.nr]';
            y = (vals)+(MeanViscs(i)-mean(vals));
            
            % Append values, track sample IDs
            x = i + 0.3*(rand(size(vals))-0.5);
        
            % Assign consistent color per sample
            c = colors(j, :);
        
            scatter(x, y, 40, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
                    'MarkerFaceAlpha', 0.7)
            
        catch
            disp(append('Failed to load data from ', Folder(j).folder, filesep, Folder(j).name, filesep, 'msadRes.mat'));
        end
        hold on
        disp(num2str(median(y)))
    end
    hold on
    FileNames{1,i} = replace(FileNames{1,i}, "_", " ");
end

set(gca, 'XTick', 1:numel(FileNames), 'XTickLabel', Labels, 'XLim', [0.5 numel(FileNames)+0.5])
xlabel('Glycerol content (v/v %)','FontSize',12,'FontWeight','bold');
ylabel('Viscosity (cP)','FontSize',12,'FontWeight','bold');
title('Viscosity scatter across conditions')
% set(gca,'YScale','log');
set(gca,'FontSize',11,'LineWidth',1.2); 
box on



