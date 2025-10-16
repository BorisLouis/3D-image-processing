%% dualcolor_error_analysis.m
% Compares raw vs SR-corrected traces between two channels.
% Computes mismatch statistics and plots one example trace.

clear; clc; close all;

%% -------------------- USER SETTINGS --------------------
rootFolder = 'S:\Dual Color\20250121_dualcolor\Multicolor_particles\In_water';

file_raw1 = fullfile(rootFolder,'traces3D_noSRCal1.mat');
file_raw2 = fullfile(rootFolder,'traces3D_noSRCal2.mat');
% file_sr1  = fullfile(rootFolder,'traces3D_1.mat');
% file_sr2  = fullfile(rootFolder,'traces3D_2.mat');

minTraceLenStats = 30;   % only traces longer than this used for stats
minTraceLenPlot  = 30;   % pick one trace longer than this for plotting
distThresh = 500;       % nm, max allowed distance when pairing detections
%% --------------------------------------------------------

%% --- Helper: extract [x y z t] from trace entry ---
function pts = extract_xyzt(entry)
    if istable(entry)
        row = entry{:,1}; col = entry{:,2}; z = entry{:,3}; t = entry{:,10};
    else
        row = entry(:,1); col = entry(:,2); z = entry(:,3); t = entry(:,10);
    end
    pts = [col, row, z, t]; % [x y z t]
end

%% --- Helper: match two trace sets (channel1 vs channel2) ---
function [allDiffs, perTrace] = match_and_diff(allTraces1, allTraces2, minLen, distThresh)
    allDiffs = []; perTrace = {};
    h = waitbar(0, 'initalizing');
    for i = 1:size(allTraces1,1)
        waitbar(i./size(allTraces1,1), h, append('Matching ', num2str(i), '/', num2str(size(allTraces1,1))))
        T1 = extract_xyzt(allTraces1{i,1});
        if size(T1,1) < minLen, continue; end
        bestMatch = [];
        for j = 1:size(allTraces2,1)
            T2 = extract_xyzt(allTraces2{j,1});
            commonT = intersect(T1(:,4), T2(:,4));
            if numel(commonT) < minLen, continue; end
            % align by common timestamps
            [~,ia,ib] = intersect(T1(:,4), T2(:,4));
            P1 = T1(ia,1:3);
            P2 = T2(ib,1:3);
            % nearest neighbor within frame
            diffs = P1 - P2;
            dists = sqrt(sum(diffs.^2,2));
            if median(dists) < distThresh
                if isempty(bestMatch) || numel(commonT) > size(bestMatch.P1,1)
                    bestMatch = struct('P1',P1,'P2',P2,'diffs',diffs,'dists',dists);
                end
            end
        end
        
        if ~isempty(bestMatch)
            allDiffs = [allDiffs; bestMatch.diffs]; %#ok<AGROW>
            perTrace{end+1} = bestMatch; %#ok<AGROW>
        end
    end
    close(h)
end

%% --- Load traces ---
s1 = load(file_raw1,'allTraces'); raw1 = s1.allTraces;
s2 = load(file_raw2,'allTraces'); raw2 = s2.allTraces;
% s3 = load(file_sr1,'allTraces');  sr1  = s3.allTraces;
% s4 = load(file_sr2,'allTraces');  sr2  = s4.allTraces;

% fprintf('Loaded %d raw ch1, %d raw ch2, %d SR ch1, %d SR ch2 traces\n', ...
%     size(raw1,1), size(raw2,1), size(sr1,1), size(sr2,1));

%% --- Compute mismatches ---
[diffs_raw, traces_raw] = match_and_diff(raw1, raw2, minTraceLenStats, distThresh);
% [diffs_sr,  traces_sr ] = match_and_diff(sr1,  sr2,  minTraceLenStats, distThresh);

% stats
stats_raw.mean3D = mean(sqrt(sum(diffs_raw.^2,2)));
stats_raw.median3D = median(sqrt(sum(diffs_raw.^2,2)));
stats_raw.std3D = std(sqrt(sum(diffs_raw.^2,2)));
% stats_sr.mean3D = mean(sqrt(sum(diffs_sr.^2,2)));
% stats_sr.median3D = median(sqrt(sum(diffs_sr.^2,2)));
% stats_sr.std3D = std(sqrt(sum(diffs_sr.^2,2)));

fprintf('\n--- 3D mismatch statistics ---\n');
fprintf('RAW: mean %.2f, median %.2f, std %.2f nm\n', ...
    stats_raw.mean3D, stats_raw.median3D, stats_raw.std3D);
fprintf(' SR: mean %.2f, median %.2f, std %.2f nm\n', ...
    stats_sr.mean3D, stats_sr.median3D, stats_sr.std3D);

fprintf('\nPer-axis (mean ± std) [nm]:\n');
fprintf('RAW dx: %.2f ± %.2f, dy: %.2f ± %.2f, dz: %.2f ± %.2f\n', ...
    mean(diffs_raw(:,1)), std(diffs_raw(:,1)), ...
    mean(diffs_raw(:,2)), std(diffs_raw(:,2)), ...
    mean(diffs_raw(:,3)), std(diffs_raw(:,3)));
% fprintf(' SR dx: %.2f ± %.2f, dy: %.2f ± %.2f, dz: %.2f ± %.2f\n', ...
%     mean(diffs_sr(:,1)), std(diffs_sr(:,1)), ...
%     mean(diffs_sr(:,2)), std(diffs_sr(:,2)), ...
%     mean(diffs_sr(:,3)), std(diffs_sr(:,3)));

%% --- Pick one long raw trace for plotting ---
chosen = [];
for i = 1:size(raw1,1)
    T1_raw = extract_xyzt(raw1{i,1});
    if size(T1_raw,1) < minTraceLenPlot, continue; end
    for j = 1:size(raw2,1)
        T2_raw = extract_xyzt(raw2{j,1});
        [commonT, ia, ib] = intersect(T1_raw(:,4), T2_raw(:,4));
        if numel(commonT) >= minTraceLenPlot
            chosen.P1 = T1_raw(ia,1:3);
            chosen.P2 = T2_raw(ib,1:3);
            chosen.T  = commonT;
            break;
        end
    end
    if ~isempty(chosen), break; end
end

if isempty(chosen)
    warning('No raw trace with >= %d points found for plotting.', minTraceLenPlot);
else
    %% --- Find the same particle in SR dataset by timestamp + proximity ---
    srTolerance = 500; % nm, allowed deviation
    
    candidateSR1 = [];
    for i = 1:size(sr1,1)
        T1_sr = extract_xyzt(sr1{i,1});
        [commonT, ia, ib] = intersect(chosen.T, T1_sr(:,4));
        if numel(commonT) >= minTraceLenPlot
            diffs = chosen.P1(ismember(chosen.T,commonT),:) - T1_sr(ib,1:3);
            if median(sqrt(sum(diffs.^2,2))) < srTolerance
                candidateSR1 = T1_sr(ib,1:3);
                break;
            end
        end
    end
    
    candidateSR2 = [];
    for j = 1:size(sr2,1)
        T2_sr = extract_xyzt(sr2{j,1});
        [commonT, ia, ib] = intersect(chosen.T, T2_sr(:,4));
        if numel(commonT) >= minTraceLenPlot
            diffs = chosen.P2(ismember(chosen.T,commonT),:) - T2_sr(ib,1:3);
            if median(sqrt(sum(diffs.^2,2))) < srTolerance
                candidateSR2 = T2_sr(ib,1:3);
                break;
            end
        end
    end
    
    %% --- Plot raw ---
    figure;
    plot3(chosen.P1(:,1),chosen.P1(:,2),chosen.P1(:,3),'g.-','LineWidth',1.5); hold on;
    plot3(chosen.P2(:,1),chosen.P2(:,2),chosen.P2(:,3),'r.-');
    xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
    title('Trace overlay BEFORE SR correction');
    legend('Channel 1','Channel 2 raw');
    grid on; axis equal; view(3);
    axLimits = axis; camView = get(gca,'View');
    
    %% --- Plot SR (same trace, matched by coords) ---
    if ~isempty(candidateSR1) && ~isempty(candidateSR2)
        figure;
        plot3(candidateSR1(:,1),candidateSR1(:,2),candidateSR1(:,3),'g.-','LineWidth',1.5); hold on;
        plot3(candidateSR2(:,1),candidateSR2(:,2),candidateSR2(:,3),'r.-');
        xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
        title('Trace overlay AFTER SR correction');
        legend('Channel 1','Channel 2 SR corrected');
        grid on; axis equal;
        axis(axLimits); set(gca,'View',camView);
    else
        warning('Could not find matching SR trace for the chosen raw trace.');
    end
end