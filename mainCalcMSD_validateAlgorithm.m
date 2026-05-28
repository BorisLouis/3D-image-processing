%% ======================================================================
%  validateNetwork.m
%
%  Validates the trained BiLSTM classifier against manually labelled
%  windows from manual_labels_v3.mat.
%
%  Each window in labelSeqs is a (10 x WIN_SIZE) feature matrix that was
%  already extracted by extractWindowFeatures during manual labelling —
%  so we pass it directly to predict(), with no re-extraction needed.
%
%  Outputs
%  -------
%   • Confusion matrix (counts + normalised)
%   • Per-class precision / recall / F1 in the console
%   • Per-window P(active) histogram, split by true class
%   • Scatter: P(active) vs window index, coloured by true label
%   • Hard cases: windows where human and network disagree most strongly
%
%  Author : Steven Huysecom – 2026
%% ======================================================================

clear; close all; clc;

%% ======================================================================
%  USER SETTINGS  ← adjust paths only
%% ======================================================================

LABELS_FILE  = 'D:\MultiColor - lysosome tracking\20260311-new analysis\GoodTrainingData\manual_labels_v3.mat';
NETWORK_FILE = 'D:\MultiColor - lysosome tracking\Dna_NB\network_v3.mat';

% Classification threshold (must match mainCalcMSD_underConstruct.m)
FRAME_ACTIVE_THR = 0.90;

% How many hard-case windows to plot for inspection
N_HARD_CASES = 12;


%% ======================================================================
%  LOAD
%% ======================================================================

fprintf('Loading labels:  %s\n', LABELS_FILE);
L = load(LABELS_FILE, 'labelSeqs', 'labelTargets', 'labelCount', 'activeCount');

labelSeqs    = L.labelSeqs;       % {1 x N}  each (10 x 50) single
labelTargets = L.labelTargets;    % (N x 1)  0=nonactive, 1=active
N            = numel(labelSeqs);
WIN_SIZE     = size(labelSeqs{1}, 2);
nCh          = size(labelSeqs{1}, 1);

nActive    = sum(labelTargets == 1);
nNonActive = sum(labelTargets == 0);
fprintf('  %d windows total  (%d active,  %d non-active)\n\n', ...
    N, nActive, nNonActive);

fprintf('Loading network: %s\n', NETWORK_FILE);
load(NETWORK_FILE, 'netCNN');
fprintf('  Network loaded.\n\n');

% Find which output column corresponds to 'active'
activeCol = find(string(netCNN.Layers(end).Classes) == 'active', 1);
if isempty(activeCol)
    error('Could not find "active" class in network output. Classes: %s', ...
        strjoin(string(netCNN.Layers(end).Classes), ', '));
end


%% ======================================================================
%  CLASSIFY
%% ======================================================================

fprintf('Classifying %d windows...\n', N);

probs   = predict(netCNN, labelSeqs);   % (N x 2)
pActive = probs(:, activeCol);          % (N x 1)  P(active) per window

% Hard classification at threshold
predLabel = double(pActive >= FRAME_ACTIVE_THR);   % 1=active, 0=nonactive

fprintf('  Done.\n\n');


%% ======================================================================
%  PERFORMANCE METRICS
%% ======================================================================

trueLabel = double(labelTargets(:));

TP = sum(trueLabel == 1 & predLabel == 1);
TN = sum(trueLabel == 0 & predLabel == 0);
FP = sum(trueLabel == 0 & predLabel == 1);
FN = sum(trueLabel == 1 & predLabel == 0);

accuracy  = (TP + TN) / N;

prec_act  = TP / max(TP + FP, 1);
rec_act   = TP / max(TP + FN, 1);
f1_act    = 2 * prec_act * rec_act / max(prec_act + rec_act, 1e-15);

prec_non  = TN / max(TN + FN, 1);
rec_non   = TN / max(TN + FP, 1);
f1_non    = 2 * prec_non * rec_non / max(prec_non + rec_non, 1e-15);

fprintf('==================== VALIDATION RESULTS ====================\n');
fprintf('  Threshold         :  %.2f\n', FRAME_ACTIVE_THR);
fprintf('  Overall accuracy  :  %.1f %%\n\n', accuracy * 100);
fprintf('  %-12s  %8s  %8s  %8s\n', 'Class', 'Precision', 'Recall', 'F1');
fprintf('  %s\n', repmat('-', 1, 44));
fprintf('  %-12s  %8.1f%%  %8.1f%%  %8.3f\n', 'active',    prec_act*100, rec_act*100, f1_act);
fprintf('  %-12s  %8.1f%%  %8.1f%%  %8.3f\n', 'non-active',prec_non*100, rec_non*100, f1_non);
fprintf('\n  Confusion matrix (rows=true, cols=predicted):\n');
fprintf('                  Pred non-act  Pred active\n');
fprintf('  True non-act    %6d        %6d\n', TN, FP);
fprintf('  True active     %6d        %6d\n', FN, TP);
fprintf('\n');


%% ======================================================================
%  FIGURE 1 – CONFUSION MATRIX  (MATLAB confusionchart)
%% ======================================================================

trueLabCat = categorical(trueLabel, [0 1], {'non-active','active'});
predLabCat = categorical(predLabel, [0 1], {'non-active','active'});

figure('Name','Confusion Matrix – Network vs Manual Labels', ...
       'Color','w','Position',[40 40 550 480]);
cm = confusionchart(trueLabCat, predLabCat, ...
    'Title', sprintf('Network vs Manual Labels  (threshold = %.2f)', FRAME_ACTIVE_THR), ...
    'RowSummary',    'row-normalized', ...
    'ColumnSummary', 'column-normalized');
cm.FontSize = 11;


%% ======================================================================
%  FIGURE 2 – P(active) DISTRIBUTION  split by true class
%% ======================================================================

colNon = [0.15, 0.40, 0.85];
colAct = [0.85, 0.15, 0.10];

figure('Name','P(active) distribution by true class', ...
       'Color','w','Position',[60 60 700 420]);

edges = linspace(0, 1, 51);
histogram(pActive(trueLabel == 0), edges, ...
    'Normalization','probability', ...
    'FaceColor',colNon,'FaceAlpha',0.6,'EdgeColor','none');
hold on;
histogram(pActive(trueLabel == 1), edges, ...
    'Normalization','probability', ...
    'FaceColor',colAct,'FaceAlpha',0.6,'EdgeColor','none');
xline(FRAME_ACTIVE_THR, 'k--', 'LineWidth', 2, ...
      'Label', sprintf('threshold = %.2f', FRAME_ACTIVE_THR), ...
      'LabelVerticalAlignment','bottom','FontSize',10);

xlabel('P(active)  from network','FontSize',12);
ylabel('Probability','FontSize',12);
legend({'True non-active','True active'},'Location','north','FontSize',10);
title('Network P(active) separated by manual label','FontSize',12);
grid on;


%% ======================================================================
%  FIGURE 3 – P(active) PER WINDOW  coloured by true label
%% ======================================================================

figure('Name','P(active) per window – coloured by true label', ...
       'Color','w','Position',[80 80 1100 380]);

actIdx = find(trueLabel == 1);
nonIdx = find(trueLabel == 0);

scatter(nonIdx, pActive(nonIdx), 8, colNon, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(actIdx, pActive(actIdx), 8, colAct, 'filled', 'MarkerFaceAlpha', 0.5);
yline(FRAME_ACTIVE_THR, 'k--', 'LineWidth', 1.8);

xlabel('Window index','FontSize',11);
ylabel('P(active)','FontSize',11);
ylim([0 1]);
legend({'True non-active','True active', ...
        sprintf('Threshold = %.2f', FRAME_ACTIVE_THR)}, ...
       'Location','northeast','FontSize',9);
title('Network P(active) for every labelled window','FontSize',11);
grid on;


%% ======================================================================
%  FIGURE 4 – HARD CASES
%  Windows where network and human disagree most strongly.
%  Sorted by |P(active) - trueLabel| descending.
%  Each subplot shows the 10-channel feature heatmap of that window,
%  plus the P(active) score, true label, and predicted label.
%% ======================================================================

% Score disagreement: distance of P(active) from the true binary label
disagreement = abs(pActive - trueLabel);

% False positives: true=non-active, predicted=active
fpIdx  = find(trueLabel == 0 & predLabel == 1);
[~, o] = sort(disagreement(fpIdx), 'descend');
fpIdx  = fpIdx(o);

% False negatives: true=active, predicted=non-active
fnIdx  = find(trueLabel == 1 & predLabel == 0);
[~, o] = sort(disagreement(fnIdx), 'descend');
fnIdx  = fnIdx(o);

nFP = min(N_HARD_CASES/2, numel(fpIdx));
nFN = min(N_HARD_CASES/2, numel(fnIdx));

featNames = {'StepSize','DAC','dStepSize','dDAC','StepVar', ...
             'Alpha','RgRaw','RgNorm','SpanNorm','Ctx'};

plotHardCases(fpIdx(1:nFP), pActive, trueLabel, labelSeqs, ...
    featNames, WIN_SIZE, colNon, colAct, 'False Positives  (true=non-active, predicted=active)');

plotHardCases(fnIdx(1:nFN), pActive, trueLabel, labelSeqs, ...
    featNames, WIN_SIZE, colNon, colAct, 'False Negatives  (true=active, predicted=non-active)');


%% ======================================================================
%  THRESHOLD SWEEP – F1 vs threshold
%% ======================================================================

threshVec = linspace(0, 1, 101);
f1_act_v  = zeros(size(threshVec));
f1_non_v  = zeros(size(threshVec));
acc_v     = zeros(size(threshVec));

for ti = 1:numel(threshVec)
    thr = threshVec(ti);
    pL  = double(pActive >= thr);
    tp_ = sum(trueLabel==1 & pL==1);
    tn_ = sum(trueLabel==0 & pL==0);
    fp_ = sum(trueLabel==0 & pL==1);
    fn_ = sum(trueLabel==1 & pL==0);
    pr  = tp_ / max(tp_+fp_, 1);
    re  = tp_ / max(tp_+fn_, 1);
    f1_act_v(ti) = 2*pr*re / max(pr+re, 1e-15);
    pn  = tn_ / max(tn_+fn_, 1);
    rn  = tn_ / max(tn_+fp_, 1);
    f1_non_v(ti) = 2*pn*rn / max(pn+rn, 1e-15);
    acc_v(ti) = (tp_+tn_) / N;
end

[~, bestIdx] = max(f1_act_v);
bestThr = threshVec(bestIdx);

figure('Name','F1 vs threshold','Color','w','Position',[100 100 700 400]);
plot(threshVec, f1_act_v,  '-',  'Color',colAct, 'LineWidth',2); hold on;
plot(threshVec, f1_non_v,  '-',  'Color',colNon, 'LineWidth',2);
plot(threshVec, acc_v,     'k-', 'LineWidth',1.5);
xline(FRAME_ACTIVE_THR, 'k--', 'LineWidth',1.5, ...
      'Label',sprintf('current = %.2f', FRAME_ACTIVE_THR), ...
      'LabelVerticalAlignment','bottom','FontSize',9);
xline(bestThr, '--', 'Color',colAct, 'LineWidth',1.5, ...
      'Label',sprintf('best F1_{act} = %.2f', bestThr), ...
      'LabelVerticalAlignment','top','FontSize',9);
xlabel('Threshold','FontSize',11);
ylabel('Score','FontSize',11);
legend({'F1 active','F1 non-active','Accuracy'},'Location','southwest','FontSize',10);
title('Classifier performance vs classification threshold','FontSize',11);
ylim([0 1]); grid on;

fprintf('Threshold sweep: best active F1 = %.3f at threshold = %.2f\n', ...
    f1_act_v(bestIdx), bestThr);
fprintf('=============================================================\n');


%% ======================================================================
%  HELPER
%% ======================================================================

function plotHardCases(idxList, pActive, trueLabel, labelSeqs, ...
    featNames, WIN_SIZE, colNon, colAct, titleStr)

    if isempty(idxList), return; end
    n = numel(idxList);

    nCols = min(n, 4);
    nRows = ceil(n / nCols);

    figure('Name', titleStr, 'Color','w', ...
           'Position', [120 120 nCols*280 nRows*260]);
    tl = tiledlayout(nRows, nCols, 'Padding','compact', 'TileSpacing','compact');
    title(tl, titleStr, 'FontSize', 11, 'FontWeight','bold');

    for k = 1:n
        wi  = idxList(k);
        X   = double(labelSeqs{wi});   % (10 x 50)
        p   = pActive(wi);
        lbl = trueLabel(wi);

        if lbl == 1
            borderCol = [0.85, 0.15, 0.10];   % red = true active
        else
            borderCol = [0.15, 0.40, 0.85];   % blue = true non-active
        end

        nexttile;
        imagesc(X);
        colormap(gca, 'parula');
        colorbar('FontSize',7);
        set(gca,'YTick',1:10,'YTickLabel',featNames,'FontSize',7);
        xlabel('Frame','FontSize',8);
        title(sprintf('Win %d  |  P=%.2f  |  True=%d  Pred=%d', ...
            wi, p, lbl, double(p >= 0.8)), 'FontSize',8);

        % Coloured border indicating true class
        ax = gca;
        ax.XColor = borderCol;
        ax.YColor = borderCol;
        ax.LineWidth = 2.5;
    end
end