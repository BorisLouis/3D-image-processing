%% dualcolor_analysis_with_calibration.m
% Performs mismatch analysis, trace overlay, and correction using precomputed
% calibration (from SRCalibration2.mat) instead of RANSAC or Procrustes fitting.

clear; clc; close all;

%% -------------------- USER SETTINGS --------------------
rootFolder = 'S:\Dual Color\20250121_dualcolor\Multicolor_particles\In_water';
calibrationFile = 'S:\Dual Color\20250121_dualcolor\2DCal\SRCalibration2.mat';
samplePattern = '0_min*';
minTraceCommon = 25; % minimum common timepoints for a trace to plot
%% -------------------------------------------------------

%% --- Load calibration transform ---
cal = load(calibrationFile, 'SRCal');
SRCal = cal.SRCal;
transTable = SRCal.Transformations;
Tcells = transTable.Coords2toCoords1;

% Extract T, b, c from all 8 rows
T_list = zeros(2,2,8);
b_list = zeros(8,1);
c_list = zeros(8,2);
for i = 1:8
    T_list(:,:,i) = Tcells{i,1}.T;
    b_list(i) = Tcells{i,1}.b;
    c_list(i,:) = Tcells{i,1}.c;
end

% Compute mean rotation, scale, and translation
T_mean = mean(T_list, 3);
b_mean = mean(b_list);
c_mean = mean(c_list, 1);

fprintf('Loaded calibration transform (average of 8 Procrustes fits):\n');
disp(T_mean);
fprintf('Mean scale: %.6f\n', b_mean);
fprintf('Mean translation: [%.6f %.6f]\n', c_mean);

%% helper: load SRList from particle.mat -> Nx4 matrix [x y z t]
function pts = load_particle_SRlist(particleMatPath)
    s = load(particleMatPath, 'particle');
    p = s.particle;
    SR = p.SRList;
    if istable(SR)
        row = SR{:,1}; col = SR{:,2}; z = SR{:,3}; t = SR{:,10};
    else
        row = SR(:,1); col = SR(:,2); z = SR(:,3); t = SR(:,10);
    end
    x = col; y = row;
    pts = [x(:), y(:), z(:), t(:)];
end

%% --- Step 1: Collect all particles and compute z offset ---
sampleDirs = dir(fullfile(rootFolder, samplePattern));
allPairs = [];

fprintf('Scanning particle.mat files...\n');
for sd = 1:numel(sampleDirs)
    sampleName = sampleDirs(sd).name;
    samplePath = fullfile(rootFolder, sampleName);
    c1 = fullfile(samplePath, append('data_', samplePath(end)), 'calibrated1', 'particle.mat');
    c2 = fullfile(samplePath, append('data_', samplePath(end)), 'calibrated2', 'particle.mat');
    if ~(exist(c1,'file') && exist(c2,'file')), continue; end
    
    pts1 = load_particle_SRlist(c1);
    pts2 = load_particle_SRlist(c2);

    % Match detections by frame/time
    t_common = intersect(unique(pts1(:,4)), unique(pts2(:,4)));
    for ti = 1:numel(t_common)
        tval = t_common(ti);
        P1 = pts1(pts1(:,4) == tval, :);
        P2 = pts2(pts2(:,4) == tval, :);
        if isempty(P1) || isempty(P2), continue; end

        % Nearest neighbour match (per frame)
        for k = 1:size(P1,1)
            diffs = P2(:,1:3) - P1(k,1:3);
            dists = sqrt(sum(diffs.^2, 2));
            [mind, idx] = min(dists);
            if mind < 1500 % arbitrary loose threshold
                allPairs = [allPairs; P1(k,1:3), P2(idx,1:3)];
            end
        end
    end
end

if isempty(allPairs)
    error('No matched detections found!');
end

% Compute mean z offset between channel2 and channel1
z_diff = allPairs(:,3) - allPairs(:,6);
z_offset = mean(z_diff);
fprintf('Mean z-offset (channel2 -> channel1): %.6f\n', z_offset);

%% --- Apply correction to all channel2 points and compute mismatch ---
X1 = allPairs(:,1:3);  % channel1 detections
X2 = allPairs(:,4:6);  % channel2 detections

% Correct x and y
XY2_corr = b_mean * (X2(:,1:2) * T_mean) + c_mean;
% Correct z
Z2_corr = X2(:,3) - z_offset;

X2_corr = [XY2_corr, Z2_corr];

% Calculate mismatch statistics
diffs_before = X1 - X2; % before correction
diffs_after = X1 - X2_corr; % after correction

dists_before = sqrt(sum(diffs_before.^2,2));
dists_after = sqrt(sum(diffs_after.^2,2));

fprintf('\n--- Mismatch statistics ---\n');
fprintf('Before correction: mean 3D distance = %.4f ± %.4f\n', mean(dists_before), std(dists_before));
fprintf('After correction:  mean 3D distance = %.4f ± %.4f\n', mean(dists_after), std(dists_after));

fprintf('Per-axis before (mean ± std): X: %.4f ± %.4f, Y: %.4f ± %.4f, Z: %.4f ± %.4f\n', ...
    mean(diffs_before(:,1)), std(diffs_before(:,1)), ...
    mean(diffs_before(:,2)), std(diffs_before(:,2)), ...
    mean(diffs_before(:,3)), std(diffs_before(:,3)));

fprintf('Per-axis after (mean ± std):  X: %.4f ± %.4f, Y: %.4f ± %.4f, Z: %.4f ± %.4f\n', ...
    mean(diffs_after(:,1)), std(diffs_after(:,1)), ...
    mean(diffs_after(:,2)), std(diffs_after(:,2)), ...
    mean(diffs_after(:,3)), std(diffs_after(:,3)));

%% --- Step 2: Trace selection and plotting ---
bestTracePair = struct('coords1',[],'coords2_uncorr',[],'coords2_corr',[],'sample','');

for sd = 1:numel(sampleDirs)
    sampleName = sampleDirs(sd).name;
    samplePath = fullfile(rootFolder, sampleName);
    tr1file = fullfile(samplePath, 'trackResults1.mat');
    tr2file = fullfile(samplePath, 'trackResults2.mat');
    if ~(exist(tr1file,'file') && exist(tr2file,'file')), continue; end

    d1 = load(tr1file, 'trackRes');
    d2 = load(tr2file, 'trackRes');
    tr1 = d1.trackRes.traces;
    tr2 = d2.trackRes.traces;

    for i = 1:size(tr1,1)
        T1 = tr1{i,1};
        if istable(T1)
            row = T1{:,1}; col = T1{:,2}; z = T1{:,3}; t = T1{:,10};
        else
            row = T1(:,1); col = T1(:,2); z = T1(:,3); t = T1(:,10);
        end
        coords1 = [col, row, z, t];
        for j = 1:size(tr2,1)
            T2 = tr2{j,1};
            if istable(T2)
                row2 = T2{:,1}; col2 = T2{:,2}; z2 = T2{:,3}; t2 = T2{:,10};
            else
                row2 = T2(:,1); col2 = T2(:,2); z2 = T2(:,3); t2 = T2(:,10);
            end
            coords2 = [col2, row2, z2, t2];
            [commonT, ia, ib] = intersect(coords1(:,4), coords2(:,4));
            if numel(commonT) > minTraceCommon
                P1 = coords1(ia,1:3);
                P2_uncorr = coords2(ib,1:3);
                P2_corr = [b_mean * (P2_uncorr(:,1:2) * T_mean) + c_mean, P2_uncorr(:,3) - z_offset];
                bestTracePair.coords1 = P1;
                bestTracePair.coords2_uncorr = P2_uncorr;
                bestTracePair.coords2_corr = P2_corr;
                bestTracePair.sample = sampleName;
                break;
            end
        end
        if ~isempty(bestTracePair.sample), break; end
    end
    if ~isempty(bestTracePair.sample), break; end
end

if isempty(bestTracePair.sample)
    warning('No trace found with > %d overlapping detections.', minTraceCommon);
else
    fprintf('Selected trace from sample %s\n', bestTracePair.sample);
    
    % Plot before correction
    figure;
    plot3(bestTracePair.coords1(:,1), bestTracePair.coords1(:,2), bestTracePair.coords1(:,3), 'go-', 'LineWidth', 1.5); hold on;
    plot3(bestTracePair.coords2_uncorr(:,1), bestTracePair.coords2_uncorr(:,2), bestTracePair.coords2_uncorr(:,3), 'ro-');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Trace overlay BEFORE correction');
    legend('Channel 1','Channel 2 uncorrected');
    grid on; axis equal; view(3);

    % Plot after correction
    figure;
    plot3(bestTracePair.coords1(:,1), bestTracePair.coords1(:,2), bestTracePair.coords1(:,3), 'go-', 'LineWidth', 1.5); hold on;
    plot3(bestTracePair.coords2_corr(:,1), bestTracePair.coords2_corr(:,2), bestTracePair.coords2_corr(:,3), 'ro-');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Trace overlay AFTER correction');
    legend('Channel 1','Channel 2 corrected');
    grid on; axis equal; view(3);
end

fprintf('\nDone.\n');