%% ======================================================================
%  PLOT_AND_SIMULATE_LABELS.m
%
%  PART 1 – Plot per-frame and window-mean feature distributions from
%            the real manually labelled data (manual_labels_v3.mat).
%
%  PART 2 – Generate synthetic windows with EXACTLY the same
%            distributional shape as the real data, using bootstrap
%            resampling with small controlled perturbation:
%              • Sample a real window at random (with replacement)
%              • Add independent Gaussian noise per channel per frame,
%                scaled to NOISE_FRAC × that channel's empirical std
%              • This preserves: marginal distributions, cross-channel
%                correlations (up to r=0.99 in real data), and within-
%                window temporal autocorrelation structure
%
%  PART 3 – Replot the exact same figures for the simulated data so you
%            can visually verify the match to the real distributions.
%
%  Author : Steven Huysecom – 2026
%% ======================================================================

clear; close all; clc;

%% ======================================================================
%  USER PARAMETERS
%% ======================================================================

LABELS_FILE = 'E:\MultiColor - lysosome tracking\20260311-new analysis\test2\manual_labels_v3.mat';
OUT_FILE    = 'E:\MultiColor - lysosome tracking\20260311-new analysis\test2\simulated_labels.mat';

% Number of simulated windows to generate per class.
% These are multiples of the real counts to give a large augmented set.
N_SIM_NONACTIVE = 6000;   % real non-active count: 1222
N_SIM_ACTIVE    = 2400;   % real active count:     400

% Perturbation noise level: fraction of each channel's per-frame std.
% 0.12 keeps distributions very close to real; increase for more diversity.
% Values above 0.25 will visibly widen the distributions.
NOISE_FRAC = 0.12;

rng(2026);   % reproducibility


%% ======================================================================
%  LOAD REAL DATA
%% ======================================================================

fprintf('Loading %s ...\n', LABELS_FILE);
L = load(LABELS_FILE, 'labelSeqs', 'labelTargets', 'labelCount', 'activeCount');

labelSeqs    = L.labelSeqs;      % {1 × N} cell of (10 × 50) single arrays
labelTargets = L.labelTargets;   % (N × 1) double: 0=nonactive, 1=active
labelCount   = L.labelCount;
activeCount  = L.activeCount;
nNonActive   = labelCount - activeCount;

nCh      = size(labelSeqs{1}, 1);   % 10
WIN_SIZE = size(labelSeqs{1}, 2);   % 50

fprintf('Loaded %d windows  (%d active,  %d non-active)\n', ...
    labelCount, activeCount, nNonActive);
fprintf('Channels: %d   WIN_SIZE: %d\n\n', nCh, WIN_SIZE);

act_mask = logical(labelTargets == 1);   % (N×1) logical
non_mask = logical(labelTargets == 0);

% Separate into class cell arrays
act_wins_real = labelSeqs(act_mask);   % {1 × nActive}
non_wins_real = labelSeqs(non_mask);   % {1 × nNonActive}
nAct_real     = numel(act_wins_real);
nNon_real     = numel(non_wins_real);


%% ======================================================================
%  PART 1 – PLOT REAL DATA DISTRIBUTIONS
%% ======================================================================

fprintf('=== PART 1: Plotting real data distributions ===\n');
plotDistributions(labelSeqs, act_mask, non_mask, nCh, ...
    'Real labelled data', activeCount, nNonActive);


%% ======================================================================
%  PART 2 – SIMULATE BY BOOTSTRAP RESAMPLING + PERTURBATION
%% ======================================================================

fprintf('=== PART 2: Simulating %d non-active + %d active windows ===\n', ...
    N_SIM_NONACTIVE, N_SIM_ACTIVE);

% Compute per-channel noise scale from empirical std of real per-frame values
ch_std_act = zeros(nCh, 1);
ch_std_non = zeros(nCh, 1);
for ch = 1:nCh
    vals_act = cellfun(@(w) w(ch,:), act_wins_real, 'UniformOutput', false);
    vals_non = cellfun(@(w) w(ch,:), non_wins_real, 'UniformOutput', false);
    ch_std_act(ch) = std(double([vals_act{:}]));
    ch_std_non(ch) = std(double([vals_non{:}]));
end
noise_scale_act = NOISE_FRAC .* ch_std_act;   % (nCh × 1)
noise_scale_non = NOISE_FRAC .* ch_std_non;

fprintf('  Noise scale per channel (non-active): ');
fprintf('%.2e  ', noise_scale_non); fprintf('\n');
fprintf('  Noise scale per channel (active):     ');
fprintf('%.2e  ', noise_scale_act); fprintf('\n\n');

% Pre-allocate output
N_SIM_TOTAL = N_SIM_NONACTIVE + N_SIM_ACTIVE;
sim_seqs    = cell(1, N_SIM_TOTAL);
sim_targets = zeros(N_SIM_TOTAL, 1);

% ---- Non-active windows ----------------------------------------------
for i = 1:N_SIM_NONACTIVE
    % Bootstrap: sample one real non-active window with replacement
    src = double(non_wins_real{randi(nNon_real)});   % (nCh × WIN_SIZE)

    % Add independent Gaussian noise, broadcast noise_scale across frames
    noise = randn(nCh, WIN_SIZE) .* noise_scale_non;   % (nCh × WIN_SIZE)
    sim_seqs{i}    = single(src + noise);
    sim_targets(i) = 0;
end

% ---- Active windows ---------------------------------------------------
for i = 1:N_SIM_ACTIVE
    idx = N_SIM_NONACTIVE + i;
    src = double(act_wins_real{randi(nAct_real)});

    noise = randn(nCh, WIN_SIZE) .* noise_scale_act;
    sim_seqs{idx}    = single(src + noise);
    sim_targets(idx) = 1;
end

% Shuffle so active and non-active are interleaved
perm        = randperm(N_SIM_TOTAL);
sim_seqs    = sim_seqs(perm);
sim_targets = sim_targets(perm);

labelCount_sim  = N_SIM_TOTAL;
activeCount_sim = N_SIM_ACTIVE;
act_mask_sim    = logical(sim_targets == 1);
non_mask_sim    = logical(sim_targets == 0);

fprintf('  Done: %d windows (%d active, %d non-active)\n\n', ...
    N_SIM_TOTAL, N_SIM_ACTIVE, N_SIM_NONACTIVE);


%% ======================================================================
%  PART 3 – PLOT SIMULATED DISTRIBUTIONS (same format as Part 1)
%% ======================================================================

fprintf('=== PART 3: Plotting simulated data distributions ===\n');
plotDistributions(sim_seqs, act_mask_sim, non_mask_sim, nCh, ...
    sprintf('Simulated  (bootstrap + %.0f%% noise)', NOISE_FRAC*100), ...
    activeCount_sim, N_SIM_NONACTIVE);


%% ======================================================================
%  SAVE SIMULATED LABELS
%% ======================================================================

labelSeqs    = sim_seqs;
labelTargets = sim_targets;
labelCount   = labelCount_sim;
activeCount  = activeCount_sim;

save(OUT_FILE, 'labelSeqs', 'labelTargets', 'labelCount', 'activeCount', '-v7.3');
fprintf('\nSaved simulated labels → %s\n', OUT_FILE);
fprintf('  %d windows  (%d active,  %d non-active)\n\n', ...
    labelCount, activeCount, N_SIM_NONACTIVE);

fprintf('To merge with real labels before training in Section 4:\n');
fprintf('  L = load(LABELS_FILE);\n');
fprintf('  S = load(OUT_FILE);\n');
fprintf('  labelSeqs    = [L.labelSeqs,    S.labelSeqs];\n');
fprintf('  labelTargets = [L.labelTargets; S.labelTargets];\n');
fprintf('  labelCount   = numel(labelSeqs);\n');
fprintf('  activeCount  = sum(labelTargets == 1);\n');


%% ======================================================================
%  HELPER FUNCTION
%  plotDistributions(seqs, act_mask, non_mask, nCh, titleStr, nAct, nNon)
%
%  Produces two figures:
%    Fig A – per-frame distributions (every frame value pooled)
%    Fig B – window-mean distributions (one mean per window)
%  Both with non-active in blue, active in red, dashed median lines.
%% ======================================================================

function plotDistributions(seqs, act_mask, non_mask, nCh, titleStr, nAct, nNon)

    colAct = [0.85, 0.15, 0.10];   % red
    colNon = [0.15, 0.40, 0.85];   % blue
    nBins  = 60;

    featNames = { ...
        'Ch1  Step size', ...
        'Ch2  DAC  (turn-angle cosine)', ...
        'Ch3  \Delta step size', ...
        'Ch4  \Delta DAC', ...
        'Ch5  Local step variance', ...
        'Ch6  Local \alpha / 10', ...
        'Ch7  Raw Rg  (norm)', ...
        'Ch8  Normalised Rg', ...
        'Ch9  Span / nSteps', ...
        'Ch10 Context score'};

    % Pre-separate windows by class (outside the channel loop for speed)
    act_wins = seqs(act_mask);
    non_wins = seqs(non_mask);

    % ================================================================
    % Figure A: per-frame distributions
    % ================================================================
    figure('Name', ['Per-frame – ' titleStr], ...
           'Color', 'w');

    for ch = 1:nCh
        ax = subplot(2, 5, ch);

        % Pool all per-frame values for this channel across all windows
        tmpA = cellfun(@(w) w(ch,:), act_wins, 'UniformOutput', false);
        tmpN = cellfun(@(w) w(ch,:), non_wins, 'UniformOutput', false);
        
        va = double([tmpA{:}]);
        vn = double([tmpN{:}]);

        % Clip to [0.5th, 99.5th] percentile for display
        lo = min(prctile(va, 0.5),  prctile(vn, 0.5));
        hi = max(prctile(va, 99.5), prctile(vn, 99.5));
        if hi <= lo, hi = lo + eps; end
        edges = linspace(lo, hi, nBins+1);

        % Non-active (blue) then active (red) on top
        histogram(ax, vn(vn>=lo & vn<=hi), edges, ...
            'Normalization', 'probability', ...
            'FaceColor', colNon, 'FaceAlpha', 0.55, 'EdgeColor', 'none');
        hold(ax, 'on');
        histogram(ax, va(va>=lo & va<=hi), edges, ...
            'Normalization', 'probability', ...
            'FaceColor', colAct, 'FaceAlpha', 0.55, 'EdgeColor', 'none');

        % Dashed median lines
        xline(ax, median(vn), '--', 'Color', colNon * 0.7, 'LineWidth', 1.8, ...
              'Label', sprintf('%.3g', median(vn)), ...
              'LabelHorizontalAlignment', 'left', 'FontSize', 7);
        xline(ax, median(va), '--', 'Color', colAct * 0.7, 'LineWidth', 1.8, ...
              'Label', sprintf('%.3g', median(va)), ...
              'LabelHorizontalAlignment', 'right', 'FontSize', 7);

        title(ax, featNames{ch}, 'FontSize', 9, 'FontWeight', 'bold');
        ylabel(ax, 'Probability', 'FontSize', 8);
        grid(ax, 'on'); ax.GridAlpha = 0.25;

        if ch == 1
            legend(ax, ...
                sprintf('Non-active  (n=%d)', nNon), ...
                sprintf('Active  (n=%d)', nAct), ...
                'Location', 'northeast', 'FontSize', 8);
        end
    end

    sgtitle(['Per-frame feature distributions  –  ' titleStr], ...
        'FontSize', 11, 'FontWeight', 'bold');

    % ================================================================
    % Figure B: window-mean distributions
    % ================================================================
    figure('Name', ['Window-mean – ' titleStr], ...
            'Color', 'w');

    for ch = 1:nCh
        ax = subplot(2, 5, ch);

        % One mean value per window
        wm_act = cellfun(@(w) mean(double(w(ch,:))), act_wins);
        wm_non = cellfun(@(w) mean(double(w(ch,:))), non_wins);

        lo = min(prctile(wm_act, 1),  prctile(wm_non, 1));
        hi = max(prctile(wm_act, 99), prctile(wm_non, 99));
        if hi <= lo, hi = lo + eps; end
        edges = linspace(lo, hi, 45);

        histogram(ax, wm_non(wm_non>=lo & wm_non<=hi), edges, ...
            'Normalization', 'probability', ...
            'FaceColor', colNon, 'FaceAlpha', 0.55, 'EdgeColor', 'none');
        hold(ax, 'on');
        histogram(ax, wm_act(wm_act>=lo & wm_act<=hi), edges, ...
            'Normalization', 'probability', ...
            'FaceColor', colAct, 'FaceAlpha', 0.55, 'EdgeColor', 'none');

        xline(ax, median(wm_non), '--', 'Color', colNon * 0.7, 'LineWidth', 1.8, ...
              'Label', sprintf('%.3g', median(wm_non)), ...
              'LabelHorizontalAlignment', 'left', 'FontSize', 7);
        xline(ax, median(wm_act), '--', 'Color', colAct * 0.7, 'LineWidth', 1.8, ...
              'Label', sprintf('%.3g', median(wm_act)), ...
              'LabelHorizontalAlignment', 'right', 'FontSize', 7);

        title(ax, featNames{ch}, 'FontSize', 9, 'FontWeight', 'bold');
        ylabel(ax, 'Probability', 'FontSize', 8);
        xlabel(ax, 'Window mean', 'FontSize', 7);
        grid(ax, 'on'); ax.GridAlpha = 0.25;

        if ch == 1
            legend(ax, 'Non-active', 'Active', ...
                'Location', 'northeast', 'FontSize', 8);
        end
    end

    sgtitle(['Window-mean feature distributions  –  ' titleStr], ...
        'FontSize', 11, 'FontWeight', 'bold');

end