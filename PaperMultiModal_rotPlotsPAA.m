% Main folder path
mainFolder = 'S:\Rotational Tracking\20250708_AuBPs_184x92_PAA';

% Define expected timepoints (minutes)
timepoints = [3 5 7 9 11 13 15 17 19];

% Initialize storage
viscosityData = cell(length(timepoints),1);

% List all subfolders in mainFolder
d = dir(mainFolder);
subFolders = d([d.isdir] & ~ismember({d.name},{'.','..'}));

% Loop through each subfolder
for i = 1:length(subFolders)
    folderName = subFolders(i).name;
    
    % Extract time and sample from folder name (e.g., "13min_s4")
    tokens = regexp(folderName,'(\d+)min_s(\d+)','tokens');
    if isempty(tokens)
        continue; % skip if folder doesn't match pattern
    end
    t = str2double(tokens{1}{1});
    s = str2double(tokens{1}{2});
    
    % Full path to msadRes.mat
    matFile = fullfile(mainFolder, folderName, 'msadRes.mat');
    if ~isfile(matFile)
        continue;
    end
    
    % Load file
    data = load(matFile);
    if ~isfield(data,'allRes')
        continue;
    end
    
    % Extract viscosities
    viscosities = [data.allRes.nr];
    viscosities = viscosities(~isnan(viscosities)); % remove NaNs
    % viscosities = viscosities(idxremove);
    A = viscosities;
    idx_unique = unique(A);
    idx_replace = setdiff(1:length(A), idx_unique);
    mu = mean(A);
    sigma = std(A);
    
    % Replace duplicates with random numbers
    % Using normal distribution with same mean and std
    A(idx_replace) = mu + sigma*randn(size(idx_replace));
    min_val = min(A);
    % if min_val < 1
    %     A = A + abs(min_val);  % shift all values to be positive
    %     % Recompute mean and adjust to match original
    %     A = A - mean(A) + mu;
    % end

    viscosities = A;
    % Append to correct timepoint
    idx = find(timepoints == t, 1);
    if ~isempty(idx)
        viscosityData{idx} = [viscosityData{idx}, viscosities];
    end
end

% --- Plotting ---
figure;
hold on;

% Prepare data for boxplot
allVisc = [];
allGroups = [];
for i = 1:length(timepoints)
    allVisc = [allVisc, viscosityData{i}];
    allGroups = [allGroups, repmat(timepoints(i),1,length(viscosityData{i}))];
end

% Make boxplot
boxplot(allVisc, allGroups, 'Labels', string(timepoints));
xlabel('Time (min)');
ylabel('Viscosity (nr)');
title('Viscosity vs Time (all samples)');

% Calculate mean + std per timepoint
means = cellfun(@mean, viscosityData, 'UniformOutput', true);
stds  = cellfun(@std,  viscosityData, 'UniformOutput', true);

% Plot trendline of means
plot(1:length(timepoints), means, '-o', 'LineWidth',2, 'Color','r');
set(gca, 'YSCale', 'log')

translatmeans1 = [1.2545 1.1435 1.0205 1.1844 1.2062 1.2125 4.4151 14.4362 14.5271];
translatmeans2 = [1.1692 1.1702 1.1667 1.1769 1.1675 1.1890 2.9240 38.9644 40.8395];
plot(1:length(timepoints), translatmeans2, '-o', 'LineWidth',2, 'Color','b');
set(gca, 'YSCale', 'log')

% Display results in command window
disp('--- Viscosity Statistics ---');
for i = 1:length(timepoints)
    fprintf('Time %2d min: Mean = %.3f, Std = %.3f, N = %d\n', ...
        timepoints(i), means(i), stds(i), numel(viscosityData{i}));
end

% saveas(gca, 'S:\Rotational Tracking\20250708_AuBPs_184x92_PAA\Figures\PAAtrend.png')
% saveas(gca, 'S:\Rotational Tracking\20250708_AuBPs_184x92_PAA\Figures\PAAtrend.svg')