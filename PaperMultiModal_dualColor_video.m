% ========================================================================
% Script: make_dualcolor_traces_video.m
% Purpose: Render growing 3D traces from dual-color recordings as videos.
% ========================================================================

clear; clc; close all;

% === USER INPUT ===
mainFolder = 'S:\Dual Color\20250121_dualcolor\Multicolor_particles\In_water';
sampleFolders = {'0_min1', '0_min2', '0_min3', '0_min4', '0_min5'};
fps = 100;  % frames per second
minLength = 25; % Only keep traces longer than this

% === LOOP OVER SAMPLE FOLDERS ===
for f = 1:length(sampleFolders)
    samplePath = fullfile(mainFolder, sampleFolders{f});
    fprintf('Processing folder: %s\n', samplePath);

    % Load data
    data1 = load(fullfile(samplePath, 'traces3D_1.mat'));
    data2 = load(fullfile(samplePath, 'traces3D_2.mat'));
    traces1 = data1.traces;
    traces2 = data2.traces;

    % === Filter by trace length ===
    len1 = cellfun(@(x) height(x), traces1);
    len2 = cellfun(@(x) height(x), traces2);

    validIdx1 = find(len1 > minLength);
    validIdx2 = find(len2 > minLength);

    % Assume correspondence by index (shorter set defines count)
    validIdx = intersect(validIdx1, validIdx2);
    traces1 = traces1(validIdx);
    traces2 = traces2(validIdx);
    nTraces = length(traces1);

    if nTraces == 0
        warning('No traces longer than %d found in %s. Skipping.', minLength, sampleFolders{f});
        continue;
    end

    % === Assign consistent colors ===
    cmap = lines(nTraces);  % distinct colors

    % === Find global time range ===
    allTimes = [];
    for i = 1:nTraces
        allTimes = [allTimes; traces1{i}.t; traces2{i}.t];
    end
    tMin = min(allTimes);
    tMax = max(allTimes);
    totalDuration = tMax - tMin;
    nFrames = ceil(totalDuration * fps / (1 / fps));  % ensures 100 frames/sec

    % === Setup Video Writers ===
    aviName = fullfile(samplePath, 'dualcolor_traces.avi');
    mp4Name = fullfile(samplePath, 'dualcolor_traces.mp4');
    v1 = VideoWriter(aviName, 'Uncompressed AVI');
    v1.FrameRate = fps;
    open(v1);
    v2 = VideoWriter(mp4Name, 'MPEG-4');
    v2.FrameRate = fps;
    open(v2);

    % === Prepare Figure ===
    fig = figure('Position', [100 100 1600 800], 'Color', 'w');
    tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
    ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
    title(ax1, 'Channel 1'); title(ax2, 'Channel 2');
    xlabel(ax1, 'X'); ylabel(ax1, 'Y'); zlabel(ax1, 'Z');
    xlabel(ax2, 'X'); ylabel(ax2, 'Y'); zlabel(ax2, 'Z');
    view(ax1, 3); view(ax2, 3);

    % Determine axis limits
    allX = []; allY = []; allZ = [];
    for i = 1:nTraces
        allX = [allX; traces1{i}.col; traces2{i}.col];
        allY = [allY; traces1{i}.row; traces2{i}.row];
        allZ = [allZ; traces1{i}.z; traces2{i}.z];
    end
    lims = [min([allX allY allZ]), max([allX allY allZ])];
    xlim(ax1, [lims(1) lims(4)]);
    ylim(ax1, [lims(2) lims(5)]);
    zlim(ax1, [lims(3) lims(6)]);
    xlim(ax2, [lims(1) lims(4)]);
    ylim(ax2, [lims(2) lims(5)]);
    zlim(ax2, [lims(3) lims(6)]);

    % === RENDER VIDEO ===
    fprintf('Rendering video (%d frames, %d traces)...\n', nFrames, nTraces);
    for frame = 1:nFrames
        tNow = tMin + (frame-1)/fps;  % current time

        % Clear axes for new frame
        cla(ax1); cla(ax2);

        % Plot accumulated traces up to tNow
        for i = 1:nTraces
            % Channel 1
            tmask1 = traces1{i}.t <= tNow;
            if any(tmask1)
                plot3(ax1, traces1{i}.col(tmask1), traces1{i}.row(tmask1), ...
                    traces1{i}.z(tmask1), 'Color', cmap(i,:), 'LineWidth', 1.5);
            end

            % Channel 2
            tmask2 = traces2{i}.t <= tNow;
            if any(tmask2)
                plot3(ax2, traces2{i}.col(tmask2), traces2{i}.row(tmask2), ...
                    traces2{i}.z(tmask2), 'Color', cmap(i,:), 'LineWidth', 1.5);
            end
        end

        % Capture and write frame
        frameImg = getframe(fig);
        writeVideo(v1, frameImg);
        writeVideo(v2, frameImg);
    end

    close(v1); close(v2);
    fprintf('Saved video to:\n  %s\n  %s\n\n', aviName, mp4Name);
end

disp('✅ All videos successfully generated!');
