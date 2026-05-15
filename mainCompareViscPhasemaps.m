%% comparePhaseViscosity.m
%
% Compares PhaseMap (3D+time, multiplane) to ViscosityMap (2D, widefield)
% for every cell, to assess whether phase structure correlates with viscosity
% structure inside the cytoplasm ROI.
%
% Folder convention assumed:
%   mainFolder_Multiplane/
%       <cellID>/
%           PhaseMovie/
%               PhaseMap.mat   -> variable: QPmapTimeAv  (x, y, z, t) or (x, y, z)
%
%   mainFolder_Widefield/
%       <cellID>/
%           ViscosityMap.mat   -> variable: ViscMap  (x, y) with NaN outside ROI
%
% Outputs per cell:
%   - Figure with side-by-side maps + overlay
%   - Saved PNG to resultsDir
%   - Metrics printed to console and saved in results struct
% -------------------------------------------------------------------------

clear; clc; close all;

%% ---- USER CONFIG -------------------------------------------------------
mainFolder_Multiplane = 'D:\Data Hannah\testdata Hannah\Multiplane';
mainFolder_Widefield  = 'D:\Data Hannah\testdata Hannah\Widefield';
resultsDir            = 'D:\Data Hannah\testdata Hannah\Results';
% -------------------------------------------------------------------------

if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

%% Find matching cell folders
mpCells = getCellFolders(mainFolder_Multiplane);
wfCells = getCellFolders(mainFolder_Widefield);

% Canonicalise names by stripping known channel-specific tags that appear
% in one modality but not the other (e.g. "_GEMview_", "_widefield_", ...).
% Add more patterns to stripPatterns if your naming evolves.
stripPatterns = {'_GEMview_', '_Widefield_', '_widefield_', '_Multiplane_', '_multiplane_'};

mpCanon = canonicalNames(mpCells.names, stripPatterns);
wfCanon = canonicalNames(wfCells.names, stripPatterns);

% Match canonical names
[~, idx_mp, idx_wf] = intersect(mpCanon, wfCanon);
fprintf('Found %d matching cells.\n\n', numel(idx_mp));

% Log any unmatched folders so the user knows
unmatched_mp = setdiff(1:numel(mpCells.names), idx_mp);
unmatched_wf = setdiff(1:numel(wfCells.names), idx_wf);
for k = unmatched_mp
    fprintf('[WARNING] No widefield match for multiplane folder: %s\n', mpCells.names{k});
end
for k = unmatched_wf
    fprintf('[WARNING] No multiplane match for widefield folder: %s\n', wfCells.names{k});
end

allResults = [];

for ci = 1:numel(idx_mp)
    cellID  = mpCells.names{idx_mp(ci)};   % use multiplane name as canonical ID
    mpPath  = mpCells.paths{idx_mp(ci)};
    wfPath  = wfCells.paths{idx_wf(ci)};

    fprintf('=== Cell: %s ===\n', cellID);

    %% Load PhaseMap
    phaseFile = findFile(mpPath, 'PhaseMap.mat', 'PhaseMovie');
    if isempty(phaseFile)
        warning('No PhaseMap.mat found for %s, skipping.', cellID); continue
    end
    phData = load(phaseFile);
    phVarName = fieldnames(phData); phVarName = phVarName{1};
    P = double(phData.(phVarName));   % (x, y, z) or (x, y, z, t)

    %% Load ViscosityMap
    viscFile = findFile(wfPath, 'ViscosityMap.mat', '');
    if isempty(viscFile)
        warning('No ViscosityMap.mat found for %s, skipping.', cellID); continue
    end
    vData = load(viscFile);
    vVarName = fieldnames(vData); vVarName = vVarName{1};
    V = vData.(vVarName);
    % Unwrap nested cell/object BEFORE casting to double
    while iscell(V) && numel(V) == 1, V = V{1}; end
    if isstruct(V) && numel(V) == 1
        fn = fieldnames(V); V = V.(fn{1});
    end
    V = squeeze(double(V));


    %% Reduce PhaseMap to 2D
    ph2D = reducePhase(P);

    %% ROI mask from ViscosityMap (non-NaN pixels)
    roi = ~isnan(V);
    fprintf('  ViscMap size: %dx%d  |  ROI pixels: %d (%.1f%%)\n', ...
        size(V,1), size(V,2), sum(roi(:)), 100*mean(roi(:)));

    %% Binarize phase map via local variance (Otsu threshold)
    % Phase values inside the cell are not reliably above background,
    % so we threshold local std dev: high inside cell, flat outside.
    ph_smooth  = imgaussfilt(ph2D, 2);
    ph_var     = stdfilt(ph_smooth, ones(9));
    thresh     = graythresh(mat2gray(ph_var));
    ph2D_bin   = mat2gray(ph_var) > thresh/2;
    ph2D_bin   = imfill(ph2D_bin, 'holes');
    ph2D_bin   = imopen(ph2D_bin,  strel('disk', 3));
    ph2D_bin   = imclose(ph2D_bin, strel('disk', 5));
    ph2D_bin   = imfill(ph2D_bin,  'holes');
    cc_bin = bwconncomp(ph2D_bin);
    if cc_bin.NumObjects > 1
        sizes_bin = cellfun(@numel, cc_bin.PixelIdxList);
        [~, mi]   = max(sizes_bin);
        tmp       = false(size(ph2D_bin));
        tmp(cc_bin.PixelIdxList{mi}) = true;
        ph2D_bin  = tmp;
    end
    fprintf('  Phase binary mask: %d pixels (%.1f%% of image)\n', ...
        sum(ph2D_bin(:)), 100*mean(ph2D_bin(:)));

    %% Register: binary mask -> binary ROI (deterministic, no optimizer)
    [~, tform] = registerPhaseToVisc(ph2D_bin, roi, size(V));
    % Apply the same transform to the raw phase values for metrics/figures
    ph2D_reg   = applyTform(ph2D, tform, size(V), median(ph2D(:)));

    %% Compute structural similarity metrics (inside ROI only)
    metrics = computeMetrics(ph2D_reg, V, roi);
    fprintf('  SSIM  (structural similarity): %.4f\n', metrics.ssim);
    fprintf('  Spearman r (rank correlation): %.4f  (p=%.4g)\n', ...
        metrics.spearman_r, metrics.spearman_p);
    fprintf('  Pearson r  (linear)          : %.4f  (p=%.4g)\n', ...
        metrics.pearson_r,  metrics.pearson_p);
    fprintf('  Gradient similarity          : %.4f\n', metrics.grad_sim);
    fprintf('\n');

    entry.cellID  = cellID;
    entry.metrics = metrics;
    entry.tform   = tform;
    allResults    = [allResults, entry];

    %% Figure - overview (saved to PNG)
    figPath = fullfile(resultsDir, sprintf('%s_comparison.png', cellID));
    makeFigure(ph2D, ph2D_reg, V, roi, metrics, cellID, figPath);

    %% Scatter plot - pixel phase vs viscosity (displayed + saved)
    scatPath = fullfile(resultsDir, sprintf('%s_scatter.png', cellID));
    makeScatterFigure(ph2D_reg, V, roi, metrics, cellID, scatPath);
end

%% Summary table
if ~isempty(allResults)
    fprintf('\n========== SUMMARY ==========\n');
    fprintf('%-20s  %7s  %10s  %10s  %10s\n', ...
        'Cell','SSIM','Spearman_r','Pearson_r','Grad_sim');
    for ci = 1:numel(allResults)
        m = allResults(ci).metrics;
        fprintf('%-20s  %7.4f  %10.4f  %10.4f  %10.4f\n', ...
            allResults(ci).cellID, m.ssim, m.spearman_r, m.pearson_r, m.grad_sim);
    end
    save(fullfile(resultsDir, 'allResults.mat'), 'allResults');
    fprintf('\nResults saved to %s\n', resultsDir);
end


%% =======================================================================
%  LOCAL FUNCTIONS
%  =======================================================================

function ph2D = reducePhase(P)
% Reduce 4D (x,y,z,t) or 3D (x,y,z) phase volume to a single 2D image.
%
% Strategy: use the z-mean projection, then optionally compare with
% middle-slice and max-variance-slice approaches.
%
% z-mean projection is preferred because:
%   - It integrates signal across all focal planes (improves SNR)
%   - Structures present throughout z appear consistently
%   - Avoids arbitrary slice selection
%   - For a cell this is the integrated phase column, analogous to
%     widefield integration along optical axis

    if ndims(P) == 4
        P = mean(P, 4);   % time-average first if needed (already done if QPmapTimeAv)
    end
    % P is now (x, y, z)
    ph2D = mean(P, 3);    % z mean projection
end


function [ph_bin_reg, tform] = registerPhaseToVisc(ph2D_bin, roi, V_size)
% Binary shape registration: align phase cell mask onto viscosity ROI.
%
% Inputs (both logical binary):
%   ph2D_bin : cell mask from phase map  (size of ph2D)
%   roi      : cytoplasm ROI from visc map (size of V)
%   V_size   : [rows cols] of the viscosity image
%
% Alignment: area ratio -> isotropic scale; centroid difference -> translation.
% Fully deterministic, no optimizer.

    Rout = imref2d(V_size);

    [ph_rows, ph_cols] = find(ph2D_bin);
    ph_cy   = mean(ph_rows);
    ph_cx   = mean(ph_cols);
    ph_area = numel(ph_rows);

    [v_rows, v_cols] = find(roi);
    v_cy   = mean(v_rows);
    v_cx   = mean(v_cols);
    v_area = numel(v_rows);

    scale = sqrt(v_area / ph_area);
    tx = v_cx - scale * ph_cx;   % affine2d: x=col
    ty = v_cy - scale * ph_cy;   % affine2d: y=row
    T     = [scale 0 0; 0 scale 0; tx ty 1];
    tform = affine2d(T);

    ph_bin_reg = imwarp(double(ph2D_bin), tform, 'OutputView', Rout) > 0.5;
    overlap    = sum(ph_bin_reg(:) & roi(:)) / sum(roi(:));
    fprintf('    Binary mask overlap after registration: %.1f%%\n', overlap*100);
    if overlap < 0.5
        warning('Registration overlap <50%% - check phase binarization manually.');
    end
end


function ph_reg = applyTform(ph2D, tform, V_size, fillVal)
% Apply precomputed affine transform to raw (non-binary) phase map.
    Rout   = imref2d(V_size);
    ph_reg = imwarp(ph2D, tform, 'OutputView', Rout, 'FillValues', fillVal);
end


function metrics = computeMetrics(ph_reg, V, roi)
% Compute structural similarity metrics between registered phase and viscosity,
% evaluated ONLY inside the cytoplasm ROI.
%
% Metrics:
%  1. SSIM  - Structural Similarity Index (captures luminance, contrast,
%             structure). Insensitive to global offset/scale -> good for
%             comparing structural features with different absolute values.
%  2. Spearman rank correlation - non-parametric, robust to monotonic
%             transforms of either image. Best for comparing 'order' of
%             pixel values (i.e., which pixels are high/low) regardless
%             of units.
%  3. Pearson correlation - linear, sensitive to absolute scale differences
%             but shows whether there is a linear relationship.
%  4. Gradient similarity - compares the spatial gradient magnitude maps.
%             Edges/boundaries present in both images will score high.
%             Particularly useful for structural feature comparison.

    % Extract ROI pixels
    ph_roi = ph_reg(roi);
    v_roi  = V(roi);

    % Normalize each to zero-mean unit-variance for fair comparison
    ph_z = (ph_roi - mean(ph_roi)) / (std(ph_roi) + eps);
    v_z  = (v_roi  - mean(v_roi))  / (std(v_roi)  + eps);

    % 1. SSIM on normalized full images (masked outside ROI = 0)
    ph_norm_img = zeros(size(V));
    v_norm_img  = zeros(size(V));
    ph_norm_img(roi) = ph_z;
    v_norm_img(roi)  = v_z;
    % Scale to [0,1] for ssim function
    ph_ssim = mat2gray(ph_norm_img);
    v_ssim  = mat2gray(v_norm_img);
    [ssim_val, ~] = ssim(ph_ssim, v_ssim);
    metrics.ssim = ssim_val;

    % 2. Spearman rank correlation
    n = numel(ph_roi);
    ph_rank = tiedrank(ph_roi);
    v_rank  = tiedrank(v_roi);
    d2 = (ph_rank - v_rank).^2;
    rs = 1 - 6*sum(d2) / (n*(n^2-1));
    % p-value via t approximation
    t_stat = rs * sqrt((n-2)/(1-rs^2+eps));
    p_s = 2 * (1 - tcdf(abs(t_stat), n-2));
    metrics.spearman_r = rs;
    metrics.spearman_p = p_s;

    % 3. Pearson correlation
    rp = corr(ph_roi, v_roi, 'type', 'Pearson');
    t_stat_p = rp * sqrt((n-2)/(1-rp^2+eps));
    p_p = 2 * (1 - tcdf(abs(t_stat_p), n-2));
    metrics.pearson_r = rp;
    metrics.pearson_p = p_p;

    % 4. Gradient similarity
    % Fill NaN outside ROI with local mean before gradient so NaN->0
    % boundary doesn't create spurious edges
    V_filled = V;
    V_filled(~roi) = mean(V(roi));
    ph_filled = ph_reg;
    ph_filled(~roi) = mean(ph_reg(roi));
    [gx_ph, gy_ph] = imgradientxy(ph_filled, 'sobel');
    [gx_v,  gy_v]  = imgradientxy(V_filled,  'sobel');
    gmag_ph = sqrt(gx_ph.^2 + gy_ph.^2);
    gmag_v  = sqrt(gx_v.^2  + gy_v.^2);
    % Restrict to ROI pixels only
    ph_g_roi = gmag_ph(roi);
    v_g_roi  = gmag_v(roi);
    % Pearson correlation of gradient magnitudes (more interpretable than cosine)
    n_ph = norm(ph_g_roi); n_v = norm(v_g_roi);
    if n_ph < eps || n_v < eps
        grad_sim = NaN;
        warning('grad_sim: zero-norm gradient vector, returning NaN.');
    else
        grad_sim = dot(ph_g_roi - mean(ph_g_roi), v_g_roi - mean(v_g_roi)) / ...
                   (n_ph * n_v + eps);
    end
    metrics.grad_sim = grad_sim;
end


function makeFigure(ph2D_orig, ph2D_reg, V, roi, metrics, cellID, savePath)
% Generate a 5-panel figure:
%   1. Original phase map (z-mean projection, original size)
%   2. Registered phase map (same frame as viscosity)
%   3. Viscosity map with ROI
%   4. Overlay: registered phase (R channel) + viscosity (G channel), ROI only
%   5. Pixel scatter plot: phase vs viscosity with trendline

    fig = figure('Visible','off','Position',[100 100 1500 320]);

    % Panel 1 – original phase
    subplot(1,5,1);
    imagesc(ph2D_orig); axis image off; colormap(gca, 'parula');
    colorbar; title('Phase map (z-mean)', 'Interpreter','none');
    xlabel('[px]');

    % Panel 2 – registered phase
    subplot(1,5,2);
    ph_disp = ph2D_reg;
    ph_disp(~roi) = NaN;
    imagesc(ph_disp); axis image off; colormap(gca, 'parula');
    colorbar; title('Phase (registered, ROI)', 'Interpreter','none');

    % Panel 3 – viscosity
    subplot(1,5,3);
    v_disp = V; v_disp(~roi) = NaN;
    imagesc(v_disp); axis image off; colormap(gca, 'hot');
    colorbar; title('Viscosity map (ROI)', 'Interpreter','none');

    % Panel 4 – overlay (normalised)
    subplot(1,5,4);
    ph_n = normalise(ph2D_reg, roi);
    v_n  = normalise(V, roi);
    overlay = zeros([size(V), 3]);
    overlay(:,:,1) = ph_n .* roi;   % phase -> red
    overlay(:,:,2) = v_n  .* roi;   % visc  -> green
    image(overlay); axis image off;
    title(sprintf('Overlay  SSIM=%.3f  \\rho_s=%.3f', ...
        metrics.ssim, metrics.spearman_r), 'Interpreter','tex');

    % Panel 5 – pixel scatter: phase vs viscosity with trendline
    subplot(1,5,5);
    ph_roi_vals = ph2D_reg(roi);
    v_roi_vals  = V(roi);
    % Subsample for speed if many pixels
    maxPts = 5000;
    if numel(ph_roi_vals) > maxPts
        idx = randperm(numel(ph_roi_vals), maxPts);
        ph_plt = ph_roi_vals(idx);
        v_plt  = v_roi_vals(idx);
    else
        ph_plt = ph_roi_vals;
        v_plt  = v_roi_vals;
    end
    scatter(ph_plt, v_plt, 4, [0.4 0.6 0.9], 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    % Linear trendline
    p_fit = polyfit(ph_roi_vals, v_roi_vals, 1);
    x_fit = linspace(min(ph_roi_vals), max(ph_roi_vals), 100);
    plot(x_fit, polyval(p_fit, x_fit), 'r-', 'LineWidth', 2);
    xlabel('Phase (a.u.)'); ylabel('Viscosity (cP)');
    title(sprintf('Pixel scatter\n\\rho_s=%.3f  r=%.3f', ...
        metrics.spearman_r, metrics.pearson_r), 'Interpreter','tex');
    box on; grid on;

    sgtitle(sprintf('Cell: %s', strrep(cellID,'_','\_')), 'FontSize', 13);

    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
    fprintf('  Figure saved: %s\n', savePath);
end


function out = normalise(img, mask)
    out = zeros(size(img));
    vals = img(mask);
    mn = min(vals); mx = max(vals);
    out(mask) = (vals - mn) / (mx - mn + eps);
end


function canon = canonicalNames(names, stripPatterns)
% Remove channel-specific substrings so both modalities match on the
% common base name. E.g. "_GEMview_" is stripped from widefield names.
% Add extra patterns to stripPatterns in the USER CONFIG section as needed.
    canon = names;
    for i = 1:numel(canon)
        for p = 1:numel(stripPatterns)
            canon{i} = strrep(canon{i}, stripPatterns{p}, '_');
        end
        % Collapse double underscores left behind by removal
        while contains(canon{i}, '__')
            canon{i} = strrep(canon{i}, '__', '_');
        end
        % Strip trailing underscore
        if endsWith(canon{i}, '_')
            canon{i} = canon{i}(1:end-1);
        end
    end
end


function cellList = getCellFolders(rootPath)
% Return struct with .names and .paths of immediate child folders
    d = dir(rootPath);
    d = d([d.isdir] & ~startsWith({d.name},'.'));
    cellList.names = {d.name};
    cellList.paths = fullfile(rootPath, {d.name});
end



function makeScatterFigure(ph2D_reg, V, roi, metrics, cellID, savePath)
% Dedicated scatter plot: phase pixel values (x) vs viscosity pixel values (y)
% Shows all ROI pixels (subsampled for rendering) with linear trendline,
% Spearman and Pearson annotations, and marginal histograms.

    ph_roi_vals = ph2D_reg(roi);
    v_roi_vals  = V(roi);

    % Subsample for scatter rendering only (trendline uses all points)
    maxPts = 8000;
    if numel(ph_roi_vals) > maxPts
        idx    = randperm(numel(ph_roi_vals), maxPts);
        ph_plt = ph_roi_vals(idx);
        v_plt  = v_roi_vals(idx);
    else
        ph_plt = ph_roi_vals;
        v_plt  = v_roi_vals;
    end

    fig = figure('Visible','on','Position',[200 200 560 500]);

    % Main scatter axes
    ax_main = axes('Position',[0.12 0.12 0.62 0.62]);
    scatter(ax_main, ph_plt, v_plt, 8, [0.25 0.55 0.85], ...
        'filled', 'MarkerFaceAlpha', 0.25);
    hold(ax_main, 'on');

    % Linear trendline (fit on ALL ROI pixels)
    p_fit = polyfit(ph_roi_vals, v_roi_vals, 1);
    x_fit = linspace(min(ph_roi_vals), max(ph_roi_vals), 200);
    plot(ax_main, x_fit, polyval(p_fit, x_fit), 'r-', 'LineWidth', 2.5);

    xlabel(ax_main, 'Phase shift (a.u.)', 'FontSize', 11);
    ylabel(ax_main, 'Viscosity (cP)',      'FontSize', 11);
    box(ax_main, 'on'); grid(ax_main, 'on');


    % Top marginal histogram – phase distribution
    ax_top = axes('Position',[0.12 0.76 0.62 0.16]);
    histogram(ax_top, ph_roi_vals, 50, 'FaceColor',[0.25 0.55 0.85], ...
        'EdgeColor','none', 'FaceAlpha', 0.6);
    set(ax_top, 'XTickLabel',[], 'YTickLabel',[]); box(ax_top,'on');
    xlim(ax_top, ax_main.XLim);

    % Right marginal histogram – viscosity distribution
    ax_right = axes('Position',[0.76 0.12 0.16 0.62]);
    histogram(ax_right, v_roi_vals, 50, 'FaceColor',[0.85 0.35 0.25], ...
        'EdgeColor','none', 'FaceAlpha', 0.6, 'Orientation','horizontal');
    set(ax_right, 'XTickLabel',[], 'YTickLabel',[]); box(ax_right,'on');
    ylim(ax_right, ax_main.YLim);

    sgtitle(sprintf('Phase vs Viscosity — %s', strrep(cellID,'_','\_')), ...
        'FontSize', 11);

    exportgraphics(fig, savePath, 'Resolution', 150);
    fprintf('  Scatter saved: %s\n', savePath);
end


function filePath = findFile(cellPath, fileName, subFolder)
% Search for fileName inside cellPath, optionally first trying subFolder.
    filePath = '';
    if ~isempty(subFolder)
        candidate = fullfile(cellPath, subFolder, fileName);
        if exist(candidate, 'file'), filePath = candidate; return; end
    end
    % Recursive search fallback
    hits = dir(fullfile(cellPath, '**', fileName));
    if ~isempty(hits)
        filePath = fullfile(hits(1).folder, hits(1).name);
    end
end