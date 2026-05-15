%% =========================================================================
%  generateSyntheticTraces.m  –  v3
%
%  PART 1  - Learn statistical model from real manually labelled windows.
%  PART 2  - Generate synthetic windows with the same distributions.
%  PART 3  - Plot per-channel feature distributions (real vs synthetic).
%  PART 4  - Plot 100 random example XY traces per class (synthetic).
%  PART 4b - Plot 100 random example XY traces per class (real).
%
%  KEY FIX vs v2:
%  Ch2 (DAC) noise fraction reduced from 8% to 2% for the active class.
%  Real active DAC is a very sharp spike near +1 with almost no spread —
%  8% noise was enough to broaden the synthetic peak around 0 visibly.
%  All other channels remain at 8% noise.
%
%  Output files:
%    syntheticModel.mat   - learned model (reusable)
%    simulated_labels.mat - synthetic windows in manual_labels_v3 format
%
%  Author : Steven Huysecom - 2026
%% =========================================================================

clear; close all; clc;
rng(2026);

%% =========================================================================
%  USER SETTINGS
%% =========================================================================

LABELS_FILE = ...
    'E:\MultiColor - lysosome tracking\20260311-new analysis\GoodTrainingData\manual_labels_v3.mat';

MODEL_FILE  = 'syntheticModel.mat';

OUTPUT_FILE = ...
    'E:\MultiColor - lysosome tracking\20260311-new analysis\test2\simulated_labels.mat';

N_ACTIVE    = 1000;
N_NONACTIVE = 3500;

MAX_COORD = 512 * 81;   % 41472 nm  (for XY trace display only)

% Per-channel noise fraction for active resampling (fraction of per-channel std).
% Ch2 (DAC) gets a much tighter noise because real active DAC is very sharply
% concentrated near +1 — even 8% noise broadens the peak too much.
% All other channels use 8% which gives enough diversity without distorting shape.
%
%   Index:  1      2      3      4      5      6      7      8      9     10
RESAMPLE_NOISE_FRAC = ...
          [0,  0,  0,  0,  0,  0,  0,  0,  0,  0.0];
%                 ^
%                 Ch2 DAC: 2% noise only — keeps the sharp +1 peak intact


%% =========================================================================
%  PART 1 - LOAD REAL DATA AND LEARN MODEL
%% =========================================================================

fprintf('=== PART 1: Learning distributions from real data ===\n');

L = load(LABELS_FILE, 'labelSeqs', 'labelTargets', 'labelCount', 'activeCount');

seqs   = L.labelSeqs;
labels = L.labelTargets;

[nCh, WIN_SIZE] = size(seqs{1});
nSeq = numel(seqs);

fprintf('  Loaded %d windows  (active: %d,  non-active: %d)\n', ...
    nSeq, L.activeCount, nSeq - L.activeCount);
fprintf('  Channels: %d   WIN_SIZE: %d\n', nCh, WIN_SIZE);

classes    = {'nonactive', 'active'};
classLabel = [0, 1];
model      = struct();

for c = 1:2

    idx       = labels == classLabel(c);
    classSeqs = seqs(idx);
    nClass    = sum(idx);

    fprintf('\n  Class: %s  (%d windows)\n', classes{c}, nClass);

    % Pool all frames [nCh x totalFrames]
    allFrames = cell2mat(cellfun(@double, classSeqs, 'UniformOutput', false));

    model.(classes{c}).mu  = mean(allFrames, 2);
    model.(classes{c}).cov = cov(allFrames');

    % Store all raw windows as a (nClass x nCh x WIN_SIZE) array
    % for direct resampling in Part 2
    rawArr = zeros(nClass, nCh, WIN_SIZE);
    for k = 1:nClass
        rawArr(k,:,:) = double(classSeqs{k});
    end
    model.(classes{c}).rawWindows = rawArr;   % used for active resampling
    model.(classes{c}).nWindows   = nClass;

    for ch = 1:nCh

        vals = allFrames(ch, :);

        model.(classes{c}).channel(ch).values = vals(:);

        % Per-channel std (for noise scaling)
        model.(classes{c}).channel(ch).std = std(vals);

        % Empirical quantile grid
        q = linspace(0, 1, 1000);
        model.(classes{c}).channel(ch).quantiles = quantile(vals, q);
        model.(classes{c}).channel(ch).q         = q;

        % AR(1) temporal autocorrelation
        x1 = vals(1:end-1);
        x2 = vals(2:end);
        a  = corr(x1(:), x2(:), 'Rows', 'complete');
        if isnan(a), a = 0; end
        residual = x2 - a * x1;

        model.(classes{c}).channel(ch).AR1    = a;
        model.(classes{c}).channel(ch).resStd = std(residual);

        % Window-level statistics
        winMeans = cellfun(@(w) mean(double(w(ch,:))), classSeqs);
        winStds  = cellfun(@(w) std( double(w(ch,:))), classSeqs);

        model.(classes{c}).channel(ch).windowMeans = winMeans;
        model.(classes{c}).channel(ch).windowStds  = winStds;
    end
end

model.nCh      = nCh;
model.WIN_SIZE = WIN_SIZE;

save(MODEL_FILE, 'model', '-v7.3');
fprintf('\n  Saved model -> %s\n', MODEL_FILE);


%% =========================================================================
%  PART 2 - GENERATE SYNTHETIC WINDOWS
%% =========================================================================

fprintf('\n=== PART 2: Generating synthetic windows ===\n');

N_TOTAL         = N_ACTIVE + N_NONACTIVE;
syntheticSeqs   = cell(1, N_TOTAL);
syntheticLabels = zeros(N_TOTAL, 1);

classCounts = [N_NONACTIVE, N_ACTIVE];

ptr = 0;
for c = 1:2

    className  = classes{c};
    N          = classCounts(c);
    labelValue = classLabel(c);
    isActive   = (labelValue == 1);

    fprintf('  Generating %d %s windows...\n', N, className);

    mu     = model.(className).mu;
    covMat = model.(className).cov + 1e-10 * eye(nCh);

    nReal     = model.(className).nWindows;
    rawArr    = model.(className).rawWindows;   % (nReal x nCh x WIN_SIZE)

    for s = 1:N

        % ------------------------------------------------------------------
        % ACTIVE CLASS: generate by resampling real windows
        % ------------------------------------------------------------------
        if isActive

            % Pick one real active window as the base for ALL channels
            % This preserves cross-channel correlations exactly
            srcIdx = randi(nReal);
            X = squeeze(rawArr(srcIdx, :, :));   % (nCh x WIN_SIZE)
            X = double(X);

            % Add small per-frame Gaussian noise to each channel,
            % using a per-channel noise fraction so Ch2 (DAC) stays tight
            for ch = 1:nCh
                chStd  = model.(className).channel(ch).std;
                noise  = randn(1, WIN_SIZE) * chStd * RESAMPLE_NOISE_FRAC(ch);
                X(ch,:) = X(ch,:) + noise;
            end

            % Clamp bounded channels after noise
            X(2,:) = min(max(X(2,:), -1), 1);   % DAC in [-1,1]
            if nCh >= 4
                X(4,:) = min(max(X(4,:), -2), 2);
            end

        % ------------------------------------------------------------------
        % NON-ACTIVE CLASS: use full AR(1) + quantile-remap pipeline
        % (unchanged from v1 — non-active distributions matched well)
        % ------------------------------------------------------------------
        else

            % Step 1: correlated Gaussian seed
            X = mvnrnd(mu, covMat, WIN_SIZE)';   % (nCh x WIN_SIZE)

            % Step 2: AR(1) temporal dynamics per channel
            for ch = 1:nCh
                a      = model.(className).channel(ch).AR1;
                resStd = model.(className).channel(ch).resStd;
                for t = 2:WIN_SIZE
                    X(ch,t) = a * X(ch,t-1) ...
                        + sqrt(max(1e-12, 1 - a^2)) * X(ch,t) ...
                        + resStd * randn();
                end
            end

            % Step 3: quantile remapping to empirical distribution
            for ch = 1:nCh
                vals      = X(ch,:);
                q_vals    = tiedrank(vals) ./ (numel(vals) + 1);
                qGrid     = model.(className).channel(ch).q;
                quantVals = model.(className).channel(ch).quantiles;
                X(ch,:)   = interp1(qGrid, quantVals, q_vals, 'linear', 'extrap');
            end

            % Step 4: clamp bounded channels
            X(2,:) = min(max(X(2,:), -1), 1);
            if nCh >= 4
                X(4,:) = min(max(X(4,:), -2), 2);
            end

            % Step 5: match window-level mean and std to a random real window
            for ch = 1:nCh
                targetMean = datasample(model.(className).channel(ch).windowMeans, 1);
                targetStd  = datasample(model.(className).channel(ch).windowStds,  1);
                currMean   = mean(X(ch,:));
                currStd    = std(X(ch,:));
                if currStd < 1e-12, currStd = 1e-12; end
                X(ch,:) = (X(ch,:) - currMean) / currStd * targetStd + targetMean;
            end

        end

        ptr = ptr + 1;
        syntheticSeqs{ptr}   = single(X);
        syntheticLabels(ptr) = labelValue;
    end
end

% Shuffle
perm            = randperm(N_TOTAL);
syntheticSeqs   = syntheticSeqs(perm);
syntheticLabels = syntheticLabels(perm);

labelSeqs    = syntheticSeqs;
labelTargets = syntheticLabels(:);
labelCount   = N_TOTAL;
activeCount  = sum(labelTargets == 1);

save(OUTPUT_FILE, 'labelSeqs', 'labelTargets', 'labelCount', 'activeCount', '-v7.3');
fprintf('  Saved %d synthetic windows -> %s\n', N_TOTAL, OUTPUT_FILE);
fprintf('  Active: %d   Non-active: %d\n', activeCount, N_TOTAL - activeCount);


%% =========================================================================
%  PART 3 - COMPARE FEATURE DISTRIBUTIONS: REAL vs SYNTHETIC
%% =========================================================================

fprintf('\n=== PART 3: Plotting feature distributions ===\n');

featNames = { ...
    'Ch1  Step size', ...
    'Ch2  DAC', ...
    'Ch3  Delta step size', ...
    'Ch4  Delta DAC', ...
    'Ch5  Local step variance', ...
    'Ch6  Local alpha / 10', ...
    'Ch7  Raw Rg', ...
    'Ch8  Normalised Rg', ...
    'Ch9  Span / nSteps', ...
    'Ch10 Context score'};

colNon = [0.15, 0.40, 0.85];
colAct = [0.85, 0.15, 0.10];

% Collect per-channel values from real data
real_act_vals = cell(nCh,1);
real_non_vals = cell(nCh,1);
for ch = 1:nCh
    real_act_vals{ch} = model.active.channel(ch).values(:)';
    real_non_vals{ch} = model.nonactive.channel(ch).values(:)';
end

% Collect per-channel values from synthetic data
sim_act_vals = cell(nCh,1);
sim_non_vals = cell(nCh,1);
for i = 1:N_TOTAL
    X = double(syntheticSeqs{i});
    for ch = 1:nCh
        if syntheticLabels(i) == 1
            sim_act_vals{ch} = [sim_act_vals{ch}, X(ch,:)];
        else
            sim_non_vals{ch} = [sim_non_vals{ch}, X(ch,:)];
        end
    end
end

datasetLabels = {'Real labelled data', 'Synthetic data  (v2 - resampled active)'};
actValsDS     = {real_act_vals, sim_act_vals};
nonValsDS     = {real_non_vals, sim_non_vals};
actN          = {L.activeCount, N_ACTIVE};
nonN          = {nSeq - L.activeCount, N_NONACTIVE};

for ds = 1:2
    ttl = sprintf('%s  (active: %d,  non-active: %d)', ...
        datasetLabels{ds}, actN{ds}, nonN{ds});

    figure('Name', ttl, 'Color','w');
    tl = tiledlayout(2, 5, 'Padding','compact', 'TileSpacing','compact');

    for ch = 1:nCh
        nexttile;

        va = actValsDS{ds}{ch};
        vn = nonValsDS{ds}{ch};

        lo = min(prctile(va, 0.5),  prctile(vn, 0.5));
        hi = max(prctile(va, 99.5), prctile(vn, 99.5));
        if hi <= lo, hi = lo + eps; end
        edges = linspace(lo, hi, 61);

        histogram(vn(vn>=lo & vn<=hi), edges, ...
            'Normalization','probability', ...
            'FaceColor',colNon,'FaceAlpha',0.55,'EdgeColor','none');
        hold on;
        histogram(va(va>=lo & va<=hi), edges, ...
            'Normalization','probability', ...
            'FaceColor',colAct,'FaceAlpha',0.55,'EdgeColor','none');

        xline(median(vn),'--','Color',colNon*0.7,'LineWidth',1.5);
        xline(median(va),'--','Color',colAct*0.7,'LineWidth',1.5);

        title(featNames{ch},'FontSize',9,'FontWeight','bold');
        ylabel('Probability','FontSize',8);
        grid on;

        if ch == 1
            legend({'Non-active','Active'},'Location','northeast','FontSize',8);
        end
    end
    title(tl, ttl, 'FontSize',11,'FontWeight','bold');
end


%% =========================================================================
%  PART 4 - PLOT 100 EXAMPLE XY TRACES PER CLASS  (SYNTHETIC)
%% =========================================================================

fprintf('\n=== PART 4: Plotting 100 synthetic XY traces per class ===\n');

N_PLOT       = 100;
nCols        = 10;
nRows        = ceil(N_PLOT / nCols);
classColors  = {colNon, colAct};

for c = 1:2

    className  = classes{c};
    labelValue = classLabel(c);
    col        = classColors{c};

    classIdx = find(syntheticLabels == labelValue);
    if numel(classIdx) >= N_PLOT
        selIdx = classIdx(randperm(numel(classIdx), N_PLOT));
    else
        selIdx = classIdx;
    end

    figure('Name', sprintf('100 synthetic %s traces', className), 'Color','w');
    tl = tiledlayout(nRows, nCols, 'Padding','none', 'TileSpacing','none');

    for k = 1:numel(selIdx)

        X   = double(syntheticSeqs{selIdx(k)});
        ss  = X(1,:) * MAX_COORD;
        dac = min(max(X(2,:), -1), 1);

        angle = 2*pi * rand();   % random initial direction
        xp = zeros(1, WIN_SIZE+1);
        yp = zeros(1, WIN_SIZE+1);

        for t = 1:WIN_SIZE
            turn  = acos(dac(t)) * sign(randn());
            angle = angle + turn;
            xp(t+1) = xp(t) + ss(t) * cos(angle);
            yp(t+1) = yp(t) + ss(t) * sin(angle);
        end

        nexttile;
        plot(xp, yp, '-', 'Color', col, 'LineWidth', 0.8);
        axis equal off;
    end

    title(tl, sprintf('100 synthetic %s traces', className), ...
        'FontSize', 11, 'FontWeight', 'bold');
end


%% =========================================================================
%  PART 4b - PLOT 100 EXAMPLE XY TRACES FROM REAL LABELLED DATA
%% =========================================================================

fprintf('\n=== PART 4b: Plotting 100 real XY traces per class ===\n');

real_act_seqs = seqs(labels == 1);
real_non_seqs = seqs(labels == 0);

realSets      = {real_non_seqs, real_act_seqs};
realSetNames  = {'non-active', 'active'};
realSetColors = {colNon, colAct};

for c = 1:2

    classSeqsPlot = realSets{c};
    col           = realSetColors{c};
    className     = realSetNames{c};

    nAvail = numel(classSeqsPlot);
    selIdx = randperm(nAvail, min(N_PLOT, nAvail));

    figure('Name', sprintf('100 REAL %s traces', className), 'Color','w');
    tl = tiledlayout(nRows, nCols, 'Padding','none', 'TileSpacing','none');

    for k = 1:numel(selIdx)

        X   = double(classSeqsPlot{selIdx(k)});
        ss  = X(1,:) * MAX_COORD;
        dac = min(max(X(2,:), -1), 1);

        angle = 2*pi * rand();
        xp = zeros(1, WIN_SIZE+1);
        yp = zeros(1, WIN_SIZE+1);

        for t = 1:WIN_SIZE
            turn  = acos(dac(t)) * sign(randn());
            angle = angle + turn;
            xp(t+1) = xp(t) + ss(t) * cos(angle);
            yp(t+1) = yp(t) + ss(t) * sin(angle);
        end

        nexttile;
        plot(xp, yp, '-', 'Color', col, 'LineWidth', 0.8);
        axis equal off;
    end

    title(tl, sprintf('100 REAL %s traces', className), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

fprintf('\nAll done.\n');