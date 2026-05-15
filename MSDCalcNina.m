%% MSD Analysis Script — Enhanced
% Calculates MSD per trace, fits MSD(tau) = A * tau^alpha,
% classifies motion, and computes 8 parameters per trace.
%
% Active motion:       alpha > 1  (superdiffusion / directed)
% Normal diffusion:    alpha ~ 1  (Brownian)
% Subdiffusion:        alpha < 1  (confined / anomalous)
%
% Only traces with MORE than 5 data points are analysed.
%
% OUTPUT  (ResultsPath .xlsx)
%   Sheet 1 : Active traces   (alpha > 1)
%   Sheet 2 : Confined traces (alpha < 1)
%   Columns : track_id | D (um2/min) | alpha | speed (um/min) |
%             total_path (um) | end_to_end (um) |
%             max_pairwise_dist (um) | Rg (um)

clear; clc; close all;

%% -------------------------------------------------------------------------
%  1. Load data
% --------------------------------------------------------------------------
filename    = 'C:\Users\steve\OneDrive\Documenten\Time exp\Time 0\1k 1mg PIC\1K 1mg_well1_pos1_XY.xlsx';
tif_file    = 'C:\Users\steve\OneDrive\Documenten\Time exp\Time 0\1k 1mg PIC\1K 1mg_well1_pos1.tif';
VideoPath   = 'C:\Users\steve\OneDrive\Documenten\Time exp\Time 0\1k 1mg PIC\Movie';
FigPath     = 'C:\Users\steve\OneDrive\Documenten\Time exp\Time 0\1k 1mg PIC\ResultsFig.png';
ResultsPath = 'C:\Users\steve\OneDrive\Documenten\Time exp\Time 0\1k 1mg PIC\Results.xlsx';

MIN_POINTS  = 5;   % traces with <= this many points are skipped

data = readtable(filename);

track_ids = data.track_id;
times_min = data.time_min;
x_um      = data.x_um;
y_um      = data.y_um;

%% -------------------------------------------------------------------------
%  2. Loop over every trace — compute MSD, alpha, D, and motility params
% --------------------------------------------------------------------------
unique_tracks = unique(track_ids);
n_tracks      = numel(unique_tracks);

% Pre-allocate result arrays
alpha_values   = NaN(n_tracks, 1);
D_values       = NaN(n_tracks, 1);   % diffusion coefficient  [um^2/min]
speed_values   = NaN(n_tracks, 1);   % mean speed             [um/min]
pathlen_values = NaN(n_tracks, 1);   % total path length      [um]
e2e_values     = NaN(n_tracks, 1);   % end-to-end distance    [um]
maxdist_values = NaN(n_tracks, 1);   % max pairwise distance  [um]
Rg_values      = NaN(n_tracks, 1);   % radius of gyration     [um]
is_valid       = false(n_tracks, 1);

for i = 1:n_tracks

    % --- Extract and sort this trace --------------------------------------
    idx = (track_ids == unique_tracks(i));

    t = times_min(idx);
    x = x_um(idx);
    y = y_um(idx);

    [t, order] = sort(t);
    x = x(order);
    y = y(order);

    n_pts = numel(t);

    if n_pts <= MIN_POINTS
        continue;
    end

    % --- MSD for lags 1 .. floor(n_pts/4) --------------------------------
    max_lag  = max(floor(n_pts / 4), 2);
    msd      = zeros(max_lag, 1);
    lag_time = zeros(max_lag, 1);

    for lag = 1:max_lag
        dx = x(lag+1:end) - x(1:end-lag);
        dy = y(lag+1:end) - y(1:end-lag);
        msd(lag)      = mean(dx.^2 + dy.^2);
        lag_time(lag) = mean(t(lag+1:end) - t(1:end-lag));
    end

    % --- Power-law fit in log-log space -----------------------------------
    keep = isfinite(log(lag_time)) & isfinite(log(msd)) & (msd > 0) & (lag_time > 0);

    if sum(keep) < 2
        continue;
    end

    p = polyfit(log(lag_time(keep)), log(msd(keep)), 1);

    alpha_values(i) = p(1);       % slope  = alpha
    A_fit           = exp(p(2));  % intercept in linear space

    % Diffusion coefficient: MSD = 4*D*tau for 2D Brownian motion.
    % Generalised from the power-law intercept A: D = A / 4
    D_values(i) = A_fit / 4;

    % --- Motility parameters ----------------------------------------------
    step_dist = sqrt(diff(x).^2 + diff(y).^2);   % step distances [um]
    step_time = diff(t);                           % step durations [min]

    % Total path length
    pathlen_values(i) = sum(step_dist);

    % Mean speed (total path length / total elapsed time)
    total_time = t(end) - t(1);
    if total_time > 0
        speed_values(i) = pathlen_values(i) / total_time;
    end

    % End-to-end distance (first point to last point)
    e2e_values(i) = sqrt((x(end) - x(1))^2 + (y(end) - y(1))^2);

    % Maximum pairwise distance (furthest apart any two points on the trace)
    % Uses pdist from the Statistics Toolbox.
    % Fallback (no toolbox): dmat = sqrt((x-x').^2+(y-y').^2); max(dmat(:))
    all_dists = pdist([x, y]);
    if ~isempty(all_dists)
        maxdist_values(i) = max(all_dists);
    else
        maxdist_values(i) = 0;
    end

    % Radius of gyration: Rg = sqrt( mean( |r_i - r_cm|^2 ) )
    x_cm         = mean(x);
    y_cm         = mean(y);
    Rg_values(i) = sqrt(mean((x - x_cm).^2 + (y - y_cm).^2));

    is_valid(i) = true;
end

%% -------------------------------------------------------------------------
%  3. Summarise and display results
% --------------------------------------------------------------------------
valid_alphas = alpha_values(is_valid);
n_valid      = numel(valid_alphas);

if n_valid == 0
    fprintf('No valid traces found (all traces have <= %d data points).\n', MIN_POINTS);
    return
end

n_active   = sum(valid_alphas > 1);
n_subdiff  = sum(valid_alphas < 1);
pct_active = 100 * n_active / n_valid;

fprintf('============================================================\n');
fprintf('  MSD Analysis Results\n');
fprintf('============================================================\n');
fprintf('  Total traces in file          : %d\n',   n_tracks);
fprintf('  Traces analysed (> %d pts)    : %d\n',   MIN_POINTS, n_valid);
fprintf('  Traces skipped (too short)    : %d\n',   n_tracks - n_valid);
fprintf('------------------------------------------------------------\n');
fprintf('  Active motion  (alpha > 1)    : %d  (%.1f %%)\n', n_active,  pct_active);
fprintf('  Subdiffusive   (alpha < 1)    : %d  (%.1f %%)\n', n_subdiff, 100*n_subdiff/n_valid);
fprintf('------------------------------------------------------------\n');
fprintf('  Alpha  -  mean   : %.3f\n', mean(valid_alphas));
fprintf('  Alpha  -  median : %.3f\n', median(valid_alphas));
fprintf('  Alpha  -  std    : %.3f\n', std(valid_alphas));
fprintf('============================================================\n');

%% -------------------------------------------------------------------------
%  4. Export results to Excel (2 sheets)
% --------------------------------------------------------------------------
col_headers = {'track_id', 'D_um2_per_min', 'alpha', ...
               'speed_um_per_min', 'total_path_um', ...
               'end_to_end_um', 'max_pairwise_dist_um', 'Rg_um'};

% Helper: build a result table for a given logical index
build_table = @(sel) array2table( ...
    [unique_tracks(sel), ...
     D_values(sel),      alpha_values(sel), ...
     speed_values(sel),  pathlen_values(sel), ...
     e2e_values(sel),    maxdist_values(sel), ...
     Rg_values(sel)], ...
    'VariableNames', col_headers);

active_mask   = is_valid & (alpha_values > 1);
confined_mask = is_valid & (alpha_values < 1);

T_active   = build_table(active_mask);
T_confined = build_table(confined_mask);

% Delete existing file to avoid stale sheets from previous runs
if isfile(ResultsPath)
    delete(ResultsPath);
end

writetable(T_active,   ResultsPath, 'Sheet', 'Active traces');
writetable(T_confined, ResultsPath, 'Sheet', 'Confined traces');

fprintf('Results saved to:\n  %s\n', ResultsPath);
fprintf('  Sheet "Active traces"   : %d rows\n', height(T_active));
fprintf('  Sheet "Confined traces" : %d rows\n', height(T_confined));

%% -------------------------------------------------------------------------
%  5. Plots  (histogram + bar chart)
% --------------------------------------------------------------------------
FigGraph = figure('Name', 'MSD Analysis', 'Color', 'w', 'Position', [100 100 1000 420]);

subplot(1, 2, 1);
histogram(valid_alphas, 25, 'FaceColor', [0.25 0.50 0.80], 'EdgeColor', 'w');
hold on;
xline(1, 'r--', 'LineWidth', 2);
text(1.02, ylim * [0; 0.95], '\alpha = 1', 'Color', 'r', 'FontSize', 10);
xlabel('\alpha  (anomalous diffusion exponent)', 'FontSize', 11);
ylabel('Number of traces',                       'FontSize', 11);
title('Distribution of \alpha values',           'FontSize', 12);
box off;

subplot(1, 2, 2);
counts  = [n_subdiff, n_active];
labels  = {sprintf('Subdiffusive\n(\\alpha < 1)'), sprintf('Active\n(\\alpha > 1)')};
colours = [0.85 0.40 0.40; 0.30 0.70 0.45];

b = bar(counts, 'FaceColor', 'flat');
b.CData = colours;
set(gca, 'XTickLabel', labels, 'FontSize', 11);
ylabel('Number of traces', 'FontSize', 11);
title('Motion type classification', 'FontSize', 12);

for k = 1:2
    text(k, counts(k) + 0.3, sprintf('%.1f %%', 100*counts(k)/n_valid), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end
ylim([0, max(counts) * 1.2]);
box off;

sgtitle(sprintf('MSD Power-Law Fit  -  %.1f %% active motion  (\\alpha > 1)', pct_active), ...
        'FontSize', 13, 'FontWeight', 'bold');
saveas(FigGraph, FigPath);

%% -------------------------------------------------------------------------
%  6. Video: raw frames + growing, alpha-colour-coded traces
% --------------------------------------------------------------------------
if ~isfile(tif_file)
    warning('TIF file "%s" not found - skipping video generation.', tif_file);
else

% Blue-to-orange colormap (256 steps)
N_COL  = 256;
blue   = [0.09, 0.46, 0.71];
orange = [1.00, 0.55, 0.00];
cmap_bo = [linspace(blue(1), orange(1), N_COL)', ...
           linspace(blue(2), orange(2), N_COL)', ...
           linspace(blue(3), orange(3), N_COL)'];

% Map each trace alpha to a colour (fixed display range [0, 2])
ALPHA_RANGE = [0, 2];
alpha_norm  = (alpha_values - ALPHA_RANGE(1)) / diff(ALPHA_RANGE);
alpha_norm  = max(0, min(1, alpha_norm));
cidx        = max(1, round(alpha_norm * (N_COL-1)) + 1);

track_col = repmat([0.45 0.45 0.45], n_tracks, 1);
for i = 1:n_tracks
    if is_valid(i)
        track_col(i, :) = cmap_bo(cidx(i), :);
    end
end

% Pre-group pixel positions by frame for fast lookup
track_px = cell(n_tracks, 1);
for i = 1:n_tracks
    idx_t       = (track_ids == unique_tracks(i));
    fr_t        = data.frame(idx_t);
    xp_t        = data.x_px(idx_t);
    yp_t        = data.y_px(idx_t);
    [fr_t, ord] = sort(fr_t);
    track_px{i} = [fr_t, xp_t(ord), yp_t(ord)];
end

% TIF metadata
tif_info     = imfinfo(tif_file);
n_tif_frames = numel(tif_info);
img_w        = tif_info(1).Width;
img_h        = tif_info(1).Height;

all_frames   = unique(data.frame);
f_start      = min(all_frames);
f_end        = min(max(all_frames), f_start + n_tif_frames - 1);
n_vid_frames = f_end - f_start + 1;

% Figure layout: image panel + colorbar strip
CB_W   = 55;
CB_PAD = 10;
fig_w  = img_w + CB_W;
fig_h  = img_h;

fig = figure('Name', 'Cell Traces Video', 'Color', 'k', ...
             'Units', 'pixels', 'Position', [80 80 fig_w fig_h], ...
             'MenuBar', 'none', 'ToolBar', 'none');

ax_img = axes('Parent', fig, 'Units', 'pixels', ...
              'Position', [0, 0, img_w, img_h], ...
              'Color', 'k', 'XColor', 'none', 'YColor', 'none');
axis(ax_img, 'image', 'off');
hold(ax_img, 'on');

cb_h  = img_h - 2*CB_PAD;
ax_cb = axes('Parent', fig, 'Units', 'pixels', ...
             'Position', [img_w + CB_PAD, CB_PAD, 16, cb_h]);

imagesc(ax_cb, 1, linspace(ALPHA_RANGE(1), ALPHA_RANGE(2), N_COL)', ...
        permute(cmap_bo, [1 3 2]));
set(ax_cb, 'YDir', 'normal', 'XTick', [], ...
    'YAxisLocation', 'right', 'TickDir', 'out', ...
    'FontSize', 8, 'FontWeight', 'bold', ...
    'XColor', 'w', 'YColor', 'w', 'Color', 'k', 'Box', 'off');
ylabel(ax_cb, '\alpha', 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');

% First frame: create image handle + all line/dot handles
img0   = imread(tif_file, 1);
h_img  = imagesc(ax_img, img0);
colormap(ax_img, gray);
axis(ax_img, 'image', 'off');
hold(ax_img, 'on');

h_line = gobjects(n_tracks, 1);
h_dot  = gobjects(n_tracks, 1);
for i = 1:n_tracks
    if ~is_valid(i); continue; end
    c         = track_col(i, :);
    h_line(i) = plot(ax_img, NaN, NaN, '-',  'Color', [c, 0.85], 'LineWidth', 1.4);
    h_dot(i)  = plot(ax_img, NaN, NaN, 'o',  'Color', c, ...
                     'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 0.5);
end

h_txt = text(ax_img, img_w * 0.02, img_h * 0.03, '', ...
             'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', ...
             'VerticalAlignment', 'top', 'Interpreter', 'none');

% VideoWriter
vid_out           = VideoWriter(VideoPath, 'MPEG-4');
vid_out.FrameRate = 15;
vid_out.Quality   = 92;
open(vid_out);

fprintf('Rendering video  (%d frames) ...\n', n_vid_frames);

for f = f_start:f_end
    tif_idx = f - f_start + 1;
    set(h_img, 'CData', imread(tif_file, tif_idx));

    for i = 1:n_tracks
        if ~is_valid(i); continue; end
        td   = track_px{i};
        mask = td(:, 1) <= f;
        if sum(mask) < 1
            set(h_line(i), 'XData', NaN, 'YData', NaN);
            set(h_dot(i),  'XData', NaN, 'YData', NaN);
        else
            xp_tr = td(mask, 2);
            yp_tr = td(mask, 3);
            set(h_line(i), 'XData', xp_tr,      'YData', yp_tr);
            set(h_dot(i),  'XData', xp_tr(end), 'YData', yp_tr(end));
        end
    end

    set(h_txt, 'String', sprintf('frame %d', f));
    drawnow limitrate;
    writeVideo(vid_out, getframe(fig));

    if mod(tif_idx, 50) == 0
        fprintf('  ... %d / %d frames done\n', tif_idx, n_vid_frames);
    end
end

close(vid_out);
close(fig);
fprintf('Video saved  -->  %s.mp4\n', VideoPath);

end  % isfile check