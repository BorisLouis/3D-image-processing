%% ========================================================================
% learnFeatureDistributions.m
%
% Learn statistical models from manually labelled microscopy traces.
%
% INPUT:
%   seqs        : cell array, each cell = [nCh x WIN_SIZE]
%   labels      : logical or numeric vector
%                 1 = active
%                 0 = non-active
%
% OUTPUT:
%   syntheticModel.mat
%
% The script learns:
%   - empirical per-channel distributions
%   - channel covariance
%   - temporal AR(1) dynamics
%   - window-level statistics
%% ------------------------------------------------------------------------
% LOAD LABELLED DATA
% -------------------------------------------------------------------------

LABELS_FILE = ...
'E:\MultiColor - lysosome tracking\20260311-new analysis\GoodTrainingData\manual_labels_v3.mat';

L = load(LABELS_FILE, ...
    'labelSeqs', ...
    'labelTargets', ...
    'labelCount', ...
    'activeCount');

labelSeqs    = L.labelSeqs;
labelTargets = L.labelTargets;

SAVE_NAME = 'syntheticModel.mat';

%% ------------------------------------------------------------------------
% BASIC INFO
% -------------------------------------------------------------------------

seqs   = labelSeqs;
labels = labelTargets;

nSeq = numel(seqs);

[nCh, WIN_SIZE] = size(seqs{1});

classes = {'nonactive', 'active'};

%% ------------------------------------------------------------------------
% BUILD MODEL
% -------------------------------------------------------------------------

model = struct();

for c = 1:2

    if c == 1
        idx = labels == 0;
    else
        idx = labels == 1;
    end

    classSeqs = seqs(idx);

    fprintf('\nProcessing class: %s\n', classes{c});

    % ---------------------------------------------------------------------
    % Pool all frames
    % ---------------------------------------------------------------------

    allFrames = [];

    for i = 1:numel(classSeqs)
        allFrames = [allFrames, double(classSeqs{i})];
    end

    % allFrames = [nCh x totalFrames]

    model.(classes{c}).mu  = mean(allFrames, 2);
    model.(classes{c}).cov = cov(allFrames');

    % ---------------------------------------------------------------------
    % Empirical distributions
    % ---------------------------------------------------------------------

    for ch = 1:nCh

        vals = allFrames(ch, :);

        model.(classes{c}).channel(ch).values = vals(:);

        % empirical quantiles
        q = linspace(0, 1, 1000);
        model.(classes{c}).channel(ch).quantiles = ...
            quantile(vals, q);

        model.(classes{c}).channel(ch).q = q;

        % -----------------------------------------------------------------
        % AR(1) temporal dynamics
        % -----------------------------------------------------------------

        x1 = vals(1:end-1);
        x2 = vals(2:end);

        a = corr(x1(:), x2(:), 'Rows', 'complete');

        if isnan(a)
            a = 0;
        end

        residual = x2 - a * x1;

        model.(classes{c}).channel(ch).AR1 = a;
        model.(classes{c}).channel(ch).resStd = std(residual);

        % -----------------------------------------------------------------
        % Window-level stats
        % -----------------------------------------------------------------

        winMeans = zeros(numel(classSeqs),1);
        winStds  = zeros(numel(classSeqs),1);

        for k = 1:numel(classSeqs)

            tmp = double(classSeqs{k}(ch,:));

            winMeans(k) = mean(tmp);
            winStds(k)  = std(tmp);
        end

        model.(classes{c}).channel(ch).windowMeans = winMeans;
        model.(classes{c}).channel(ch).windowStds  = winStds;

    end
end

model.nCh = nCh;
model.WIN_SIZE = WIN_SIZE;

save(SAVE_NAME, 'model', '-v7.3');

%% ------------------------------------------------------------------------
% PLOT REAL FEATURE DISTRIBUTIONS
% -------------------------------------------------------------------------

featureNames = { ...
    'Ch1 Step size', ...
    'Ch2 DAC (turn-angle cosine)', ...
    'Ch3 \Delta step size', ...
    'Ch4 \Delta DAC', ...
    'Ch5 Local step variance', ...
    'Ch6 Local \alpha / 10', ...
    'Ch7 Raw Rg (norm)', ...
    'Ch8 Normalised Rg', ...
    'Ch9 Span / nSteps', ...
    'Ch10 Context score'};

figure('Name','Real labelled data','Color','w');

tiledlayout(2,5,'Padding','compact','TileSpacing','compact');

for ch = 1:nCh

    nexttile;

    valsNon = model.nonactive.channel(ch).values;
    valsAct = model.active.channel(ch).values;

    histogram(valsNon, ...
        60, ...
        'Normalization','probability', ...
        'FaceAlpha',0.5);

    hold on;

    histogram(valsAct, ...
        60, ...
        'Normalization','probability', ...
        'FaceAlpha',0.5);

    xline(mean(valsNon),'--','LineWidth',1.5);
    xline(mean(valsAct),'--','LineWidth',1.5);

    title(featureNames{ch}, 'Interpreter','tex');

    xlabel('Feature value');
    ylabel('Probability');

    legend({'Non-active','Active'}, ...
        'Location','best');

    grid on;

end

sgtitle('Feature distributions - Real labelled data');

fprintf('\nSaved model to: %s\n', SAVE_NAME);
