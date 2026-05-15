BASE_DIR = 'E:\MultiColor - lysosome tracking\20260311-new analysis\A549 mSi';
Folder   = dir(BASE_DIR);
FOLDERS  = {Folder([Folder.isdir] & ~ismember({Folder.name}, {'.','..'})).name};

fprintf('=== Rendering trace videos ===\n');
nFolders = numel(FOLDERS);

for fi = 1:nFolders

    fName   = FOLDERS{fi};
    matPath = fullfile(BASE_DIR, fName, 'Traces3D.mat');
    tifPath = fullfile(BASE_DIR, fName, 'calibrated2', 'calibratedPlane2.tif');

    % --- Skip if required files are missing ----------------------------------
    if ~isfile(matPath)
        warning('MAT file not found – skipping: %s', matPath);
        continue
    end
    if ~isfile(tifPath)
        warning('TIF file not found – skipping: %s', tifPath);
        continue
    end

    % --- Load & filter traces (keep only traces > 50 datapoints) -------------
    S           = load(matPath, 'TrackedData');
    TrackedData = S.TrackedData;
    keep        = cellfun(@(t) size(t,1), TrackedData) > 50;
    TrackedData = TrackedData(keep);
    nTraces     = numel(TrackedData);
    fprintf('\n[%d/%d]  %-32s  →  %d traces (>50 pts)\n', fi, nFolders, fName, nTraces);

    if nTraces == 0
        warning('No traces longer than 50 points – skipping: %s', fName);
        continue
    end

    % --- Build track_px: cell array of [frame, col_px, row_px] ---------------
    track_px = cell(nTraces, 1);
    for i = 1:nTraces
        T          = TrackedData{i};          % m×12 table
        frames     = T.t;                     % column 10: frame number
        xp         = T.colM;                  % column 5:  x (col in image)
        yp         = T.rowM;                  % column 4:  y (row in image)
        [frames, ord] = sort(frames);
        track_px{i}   = [frames, xp(ord), yp(ord)];
    end

    % --- Assign a unique colour to every trace --------------------------------
    cmap       = lines(nTraces);             % nTraces × 3 RGB
    track_col  = cmap;

    % --- TIF metadata ---------------------------------------------------------
    tif_info     = imfinfo(tifPath);
    n_tif_frames = numel(tif_info);
    img_w        = tif_info(1).Width;
    img_h        = tif_info(1).Height;

    % Frame range: union of all trace frames, clamped to available TIF frames
    all_frames   = cell2mat(cellfun(@(t) t(:,1), track_px, 'UniformOutput', false));
    f_start      = min(all_frames);
    f_end        = min(max(all_frames), f_start + n_tif_frames - 1);
    n_vid_frames = f_end - f_start + 1;

    % --- Figure ---------------------------------------------------------------
    fig = figure('Name', fName, 'Color', 'k', ...
                 'Units', 'pixels', 'Position', [80 80 img_w img_h], ...
                 'MenuBar', 'none', 'ToolBar', 'none', 'Visible', 'off');

    ax_img = axes('Parent', fig, 'Units', 'pixels', ...
                  'Position', [0 0 img_w img_h], ...
                  'Color', 'k', 'XColor', 'none', 'YColor', 'none');
    axis(ax_img, 'image', 'off');
    hold(ax_img, 'on');

    % First background frame
    img0  = imread(tifPath, 1);
    h_img = imagesc(ax_img, img0);
    colormap(ax_img, gray);
    clim(ax_img, [min(img0(:)), max(img0(:))]);   % auto-contrast on frame 1
    axis(ax_img, 'image', 'off');
    hold(ax_img, 'on');

    % Pre-create line + dot handles for every trace
    h_line = gobjects(nTraces, 1);
    h_dot  = gobjects(nTraces, 1);
    for i = 1:nTraces
        c         = track_col(i, :);
        h_line(i) = plot(ax_img, NaN, NaN, '-',  'Color', [c, 0.85], 'LineWidth', 1.4);
        h_dot(i)  = plot(ax_img, NaN, NaN, 'o',  'Color', c, ...
                         'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 0.5);
    end

    % Frame counter label
    h_txt = text(ax_img, img_w*0.02, img_h*0.04, '', ...
                 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', ...
                 'VerticalAlignment', 'top', 'Interpreter', 'none');

    % --- VideoWriter ----------------------------------------------------------
    vidPath           = fullfile(BASE_DIR, fName, 'tracevideo.mp4');
    vid_out           = VideoWriter(vidPath, 'MPEG-4');
    vid_out.FrameRate = 15;
    vid_out.Quality   = 92;
    open(vid_out);

    % --- Render loop ----------------------------------------------------------
    fprintf('  Rendering %d frames ...\n', n_vid_frames);

    for f = f_start:f_end
        tif_idx = f - f_start + 1;

        % Background image
        set(h_img, 'CData', imread(tifPath, tif_idx));

        % Update every trace up to current frame
        for i = 1:nTraces
            td   = track_px{i};
            mask = td(:,1) <= f;
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
            fprintf('    ... %d / %d frames done\n', tif_idx, n_vid_frames);
        end
    end

    close(vid_out);
    close(fig);
    fprintf('  Video saved  →  %s\n', vidPath);

end

fprintf('\n=== All done ===\n');