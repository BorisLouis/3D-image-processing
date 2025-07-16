clc ;
clear ;
close all;

%% USER INPUT
[FilePath, Experiment, Filename, Dimension, expTime, Temp, Radius, DiffFit, MinSize, Ext] = UserInput.CalcMSDinfoGUI;

%% Loading
f = waitbar(0, 'Initializing...');
MainFolder = dir(FilePath);
name = append(Experiment, Ext);
MainFolder([MainFolder.isdir] ~= 1) = [];

AllMovieResults = [];
Diff_In = [];
AlphaLog_In = [];
AlphaFit_In = [];
AFit_In = [];
Speed_In = [];
Diff_Out = [];
AlphaLog_Out = [];
AlphaFit_Out = [];
AFit_Out = [];
Speed_Out = [];
features = [];
SegmentList = [];
n = 0;

for j = 3 : size(MainFolder,1)

    try
        folder = dir(append(MainFolder(j).folder, filesep, append(MainFolder(j).name)));
        idx = contains({folder.name},Filename);
        folder(~idx) = [];
        
        f2Load = [folder(1).folder filesep folder(1).name];
        Path = folder(1).folder;
        
        tmpData = load(f2Load);
        name = fieldnames(tmpData);
        data = tmpData.(name{1});
        
        %% Processing
        allHeight = cellfun(@height,data(:,1));
        idx = allHeight>MinSize;
        currMov = data(idx, 1);
        if isempty(currMov)
            error(append('No traces found that are longer than MinSize (', num2str(MinSize), ' datapoints)'))
        end
        allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
            'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
            'vZ',0,'vR',0);
        if strcmp(Experiment, 'Tracking-Segmentation')
            allRes.Mask = 0;
        elseif strcmp(Experiment, 'Tracking-Phase')
            allRes.Phase = 0;
        end
        allRes(length(currMov)).msdX = [];
        maxLength = max(allHeight);
        allMSDX = zeros(length(currMov),maxLength-1);
        allMSDY = allMSDX;
        allMSDZ = allMSDY;
        allMSDR = allMSDY;
    
        fixedLen = 20;
        interpLags = logspace(-1, 1, fixedLen);
    
        for i = 1:length(currMov)
            n = n+1;
            waitbar(i./length(currMov), f, append('Calculating diffusion & doing microrheology - Movie ', num2str(j-2), ' out of ', num2str(size(MainFolder,1) - 2)));
            currPart = currMov{i};
        
            coord = [currPart.col, currPart.row, currPart.z];
            CM = mean(coord,1);
            coord = coord-CM;
    
            %inR
            if strcmp(Dimension, '3D')
                msdr = MSD.calc(coord(:, 1:3)/10^3);%convert to um;
                AvStep = MSD.getAvStepSize(coord(:, 1:3)/10^3); 
                tau = (1:length(msdr))'*expTime;
                allMSDR(i,1:length(msdr)) = msdr;
            elseif strcmp(Dimension, '2D')
                msdr = MSD.calc(coord(:, 1:2)/10^3);%convert to um;
                AvStep = MSD.getAvStepSize(coord(:, 1:2)/10^3); 
                tau = (1:length(msdr))'*expTime;
                allMSDR(i,1:length(msdr)) = msdr;
            end
            dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2 +...
                (coord(1,3)-coord(end,3))^2);
            Speed = dR/10^3/(length(coord)*expTime);
            
            MSDList{n,1} = msdr;
            MSDList{n,2} = tau;
    
            InSegment = currPart.InSegment(1);
    
            % Feature 1: Initial slope (linear fit to first 4 points)
            D_linear   = MSD.getDiffCoeff(msdr,tau,DiffFit,Dimension);
        
            % Feature 2: Log-log slope (anomalous exponent alpha)
            alpha_log   = MSD.getDiffTypeAlpha2(msdr,expTime, AvStep);
        
            % Feature 3: Interpolated MSD (log-lag domain)
            try
                interpMSD = interp1(tau, msdr, interpLags, 'linear', 'extrap');
            catch
                interpMSD = zeros(1,fixedLen);
            end
        
            % Feature 4.1: Area under MSD curve
            auc = trapz(tau, msdr);
    
            % Feature 4.2 Area under MSD curve normalized to length
            aucNorm = trapz(tau, msdr)./size(msdr, 1);
        
            % Feature 5.1: mean Curvature (2nd derivative)
            d1 = gradient(msdr, tau);
            d2 = gradient(d1, tau);
            curvatureMean = mean(abs(d2));
    
            % Feature 5.1: max Curvature (2nd derivative)
            curvatureMax = max(abs(d2));
        
            % Feature 6: Normalized MSD (MSD ./ lag)
            normMSD = msdr ./ tau;
            normMSD(isnan(normMSD) | isinf(normMSD)) = 0;
            normFeat = interp1(tau, normMSD, interpLags, 'linear', 'extrap');
        
            % Feature 7: Anomalous diffusion model fit (MSD = A * t^α)
            % log(MSD) = log(A) + α * log(t)
            logLag = log10(tau);
            logMSD = log10(msdr);
            valid = ~isnan(logLag) & ~isinf(logLag) & ~isnan(logMSD) & ~isinf(logMSD);
            if sum(valid) >= 4
                p_anom = polyfit(logLag(valid), logMSD(valid), 1);
                alpha_fit = p_anom(1);
                A_fit = 10^p_anom(2);
            else
                alpha_fit = 0;
                A_fit = 0;
            end
        
            % Concatenate all features
            featVec = [ ...
                D_linear, ...
                alpha_log, ...
                interpMSD(:)', ...
                auc, ...
                aucNorm,...
                curvatureMean, ...
                curvatureMax,...
                normFeat(:)', ...
                alpha_fit, ...
                A_fit, ...
                Speed...
            ];
        
            features = [features; featVec];
    
            if InSegment == 1
                Diff_In = [Diff_In; D_linear];
                AlphaLog_In = [AlphaLog_In, alpha_log];
                AlphaFit_In = [AlphaFit_In, alpha_fit];
                AFit_In = [AFit_In, A_fit];
                Speed_In = [Speed_In, Speed];
            elseif InSegment == 0
                Diff_Out = [Diff_Out; D_linear];
                AlphaLog_Out = [AlphaLog_Out, alpha_log];
                AlphaFit_Out = [AlphaFit_Out, alpha_fit];
                AFit_Out = [AFit_Out, A_fit];
                Speed_Out = [Speed_Out, Speed];
            end
    
            SegmentList = [SegmentList, InSegment];
            
        end
    catch
    end
end

%% Standardize features
features = fillmissing(features, 'constant', 0); % handle NaNs
features = zscore(features);

%% PCA 
% [coeff, score, ~, ~, explained] = pca(features);
tsneScore = tsne(features, 'NumDimensions', 3, 'Perplexity', 10, 'Standardize', false);

%% evaluate optimal number of clusters
maxK = 10;
% Try with silhouette method (robust for most cases)
eva = evalclusters(features, 'kmeans', 'silhouette', 'KList', 2:maxK);
optimalK = eva.OptimalK;
fprintf('Optimal number of clusters (silhouette method): %d\n', optimalK);
% Optionally, visualize the evaluation
f1 = figure;
plot(eva);
title('Cluster Evaluation (Silhouette method)');
FigName = append(FilePath, filesep, 'OptimalK.png');
saveas(f1, FigName)

optimalK = 3;
%% Kmeans for classification
% [idx, C] = kmeans(score(:,1:2), optimalK);
% [idx, C] = kmeans(features, optimalK);
[idx, C] = kmeans(tsneScore, optimalK);
f2 = figure;
% gscatter(score(:,1), score(:,2), idx);
scatter3(tsneScore(:,1), tsneScore(:,2), tsneScore(:,3), 10, idx, 'filled');
xlabel('tsne1'); ylabel('tsne2'); zlabel('tsne3');
% xlabel('PCA1'); ylabel('PCA2');
title('MSD Clustering using tSNE + K-means');
% title('MSD Clustering using PCA + K-means');
FigName = append(FilePath, filesep, 'Clustering.png');
saveas(f2, FigName)

%% Plot average MSD per cluster + see if they are in segment or not
commonLag = logspace(-1, 1, 10);  % Adjust to fit your time range
nClusters = optimalK;

clusterMSDs = cell(nClusters, 1);

for c = 1:nClusters
    clusterMSDs{c} = []; % each row will be an interpolated MSD
end

for i = 1:length(idx)
    msd = MSDList{i,1};
    lag = MSDList{i,2};
    
    if length(msd) < 3
        continue
    end

    % Interpolate
    interpMSD = interp1(lag, msd, commonLag, 'linear', 'extrap');

    % Add to cluster collection
    clusterMSDs{idx(i)} = [clusterMSDs{idx(i)}; interpMSD(:)'];
end

f3 = figure;
colors = lines(nClusters);
for c = 1:nClusters
    subplot(1,nClusters,c);
    data = clusterMSDs{c};
    
    if isempty(data)
        continue
    end

    meanMSD = mean(data, 1);
    stdMSD = std(data, 0, 1);

    % Plot shaded area using fill (std dev band)
    x = [commonLag, fliplr(commonLag)];
    y = [meanMSD + stdMSD, fliplr(meanMSD - stdMSD)];

    fill(x, y, colors(c,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(commonLag, meanMSD, '-', 'Color', colors(c,:), 'LineWidth', 2);
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    title(sprintf('Cluster %d (n = %d)', c, size(data,1)));
    xlabel('Lag time');
    ylabel('MSD');
    grid on;
end
FigName = append(FilePath, filesep, 'MSDclusters.png');
saveas(f3, FigName)

%% Check which percentage of traces in a clusters is inside segment
if length(SegmentList) ~= length(idx)
    error('Length of InSegment does not match number of traces!');
end

% Calculate percentage inside segment for each cluster
insidePerc = zeros(nClusters, 1);
counts = zeros(nClusters, 1);

for c = 1:nClusters
    clusterMask = idx == c;
    counts(c) = sum(clusterMask);
    insidePerc(c) = 100 * sum(SegmentList(clusterMask)) / counts(c);
end

% Display results
fprintf('\nCluster membership inside segments:\n');
for c = 1:nClusters
    fprintf('Cluster %d: %d traces, %.1f%% inside segment\n', ...
        c, counts(c), insidePerc(c));
end

% Optional bar chart
f4 = figure;
bar(1:nClusters, insidePerc);
xlabel('Cluster #');
ylabel('% Inside Segment');
title('Segment Membership per Cluster');
% Customize appearance
title('Feature Comparison: Inside vs Outside Segment');
xlabel('Location');
ylabel('Percentage');  % Change to your actual feature name
grid on;
set(gca, 'FontSize', 12);
FigName = append(FilePath, filesep, 'Segmentsclusters.png');
saveas(f4, FigName)

%
groupLabels = cell(size(InSegment));
groupLabels(InSegment == 1) = {'In'};
groupLabels(InSegment == 0) = {'Out'};

% Make the boxplot
f5 = figure;
allDiff = [Diff_In(:); Diff_Out(:)];
% Create a group label vector
group = [repmat({'In'}, length(Diff_In), 1); ...
         repmat({'Out'}, length(Diff_Out), 1)];
% Make the boxplot
boxplot(allDiff, group, 'Symbol', '');
title('Diff coeff: Inside vs Outside Segment');
xlabel('Location');
ylabel('Diffusion coefficient (µm^2/s^-^1)');  % Change to your actual feature name
grid on;
set(gca, 'FontSize', 12);
FigName = append(FilePath, filesep, 'Diffusion_InOut.png');
saveas(f5, FigName)

% Make the boxplot
f6 = figure;
allalpha = [AlphaLog_In(:); AlphaLog_Out(:)];
% Create a group label vector
group = [repmat({'In'}, length(Diff_In), 1); ...
         repmat({'Out'}, length(Diff_Out), 1)];
boxplot(allalpha, group, 'Symbol', '');
title('alpha coeff: Inside vs Outside Segment');
xlabel('Location');
ylabel('alpha coefficient');  % Change to your actual feature name
grid on;
set(gca, 'FontSize', 12);
FigName = append(FilePath, filesep, 'Alpha_InOut.png');
saveas(f6, FigName)

f7 = figure;
allspeed = [Speed_In(:); Speed_Out(:)];
% Create a group label vector
group = [repmat({'In'}, length(Speed_In), 1); ...
         repmat({'Out'}, length(Speed_Out), 1)];
boxplot(allspeed, group, 'Symbol', '');
title('Speed: Inside vs Outside Segment');
xlabel('Location');
ylabel('speed (µm^2/s)');  % Change to your actual feature name
grid on;
set(gca, 'FontSize', 12);
FigName = append(FilePath, filesep, 'Speed_InOut.png');
saveas(f7, FigName)

disp(append(num2str(numel(find(SegmentList == 1))./size(SegmentList, 2).*100), '% inside segments'))