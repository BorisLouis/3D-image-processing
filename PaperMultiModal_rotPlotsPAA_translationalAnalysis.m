% compute_viscosity_from_traces.m (Broersma-corrected)
% Given folder structure described by the user, this script will:
% - load Traces3DCommon.mat from each sample folder
% - for every trace longer than minPoints, compute MSD (per channel), diffusion
%   coefficient and viscosity (per channel). The trace viscosity is the mean of
%   the two channel-based viscosities.
% - group results per glycerol content and save statistics, full datasets and a
%   boxplot figure (png + svg) into the "translational analysis" folder.
%
% Uses Broersma translational hydrodynamic correction for prolate ellipsoids.
%
% USAGE: edit the top section for path / physical constants if needed and run
% in MATLAB.

clearvars; close all; clc

%% USER PARAMETERS - edit if needed
rootFolder = 'S:\Rotational Tracking\20250708_AuBPs_184x92_glycerol';
outFolder  = fullfile(rootFolder,'translational analysis');
if ~exist(outFolder,'dir'), mkdir(outFolder); end

dt = 0.01;               % exposure time in seconds (10 ms)
minPoints = 21;          % only take traces strictly longer than 20 points
T = 298.15;              % temperature in Kelvin (set to 298.15 K = 25 C)

% bipyramid approximated as a prolate spheroid with long axis = 182 nm,
% short axis = 91 nm (values taken from the user). We'll convert to meters
long_axis_nm = 182; short_axis_nm = 91;
a = (long_axis_nm/2)*1e-9;   % semi-major axis (m)
b = (short_axis_nm/2)*1e-9;  % semi-minor axis (m)

kB = 1.380649e-23;      % Boltzmann constant (J/K)

%% Broersma translational friction coefficients
p = a/b;  % aspect ratio
% Parallel and perpendicular translational friction coefficients (Broersma)
f_par  = 4*pi*a ./ (log(2*p) - 0.5);
f_perp = 8*pi*a ./ (log(2*p) + 0.5);
% Orientation-averaged friction coefficient
f_avg = (f_par + 2*f_perp)/3;
% Important: these are in units of (Pa·s·m) for multiplication by viscosity

%% Walk through glycerol content folders
dirInfo = dir(rootFolder);
isub = [dirInfo(:).isdir];
subDirs = {dirInfo(isub).name}';
subDirs = subDirs(~ismember(subDirs,{'.','..'}));
glycerolFolders = subDirs(startsWith(subDirs,'glycerol_'));

results = struct();

for gi = 1:numel(glycerolFolders)
    gname = glycerolFolders{gi};
    gpath = fullfile(rootFolder,gname);
    fprintf('Processing glycerol content folder: %s\n', gname);

    sdir = dir(gpath);
    sdir = sdir([sdir.isdir]);
    sdirs = {sdir.name};
    sdirs = sdirs(~ismember(sdirs,{'.','..'}));
    sampleFolders = sdirs(startsWith(sdirs,[gname '_sample_']));

    D_all = []; eta_all = [];

    for si = 1:numel(sampleFolders)
        spath = fullfile(gpath,sampleFolders{si});
        datafile = fullfile(spath,'Traces3DCommon.mat');
        if ~exist(datafile,'file')
            warning('Missing Traces3DCommon.mat in %s -- skipping', spath);
            continue
        end

        S = load(datafile);
        if ~isfield(S,'traces3Dcommon')
            warning('File %s does not contain ''traces3Dcommon'' variable -- skipping', datafile);
            continue
        end
        traces3Dcommon = S.traces3Dcommon;

        if istable(traces3Dcommon)
            if any(strcmp(traces3Dcommon.Properties.VariableNames,'Coord1'))
                coord1_col = traces3Dcommon.Coord1;
            else
                coord1_col = traces3Dcommon{:,1};
            end
            if any(strcmp(traces3Dcommon.Properties.VariableNames,'Coord2'))
                coord2_col = traces3Dcommon.Coord2;
            else
                coord2_col = traces3Dcommon{:,2};
            end
        else
            error('traces3Dcommon is expected to be a table. Please check %s', datafile);
        end

        nTraces = height(traces3Dcommon);
        for ti = 1:nTraces
            c1 = coord1_col{ti};
            c2 = coord2_col{ti};
            if isempty(c1) || isempty(c2), continue; end
            m1 = size(c1,1);
            m2 = size(c2,1);
            if m1 < minPoints || m2 < minPoints
                continue
            end

            xy1 = double(c1)*1e-9; % m
            xy2 = double(c2)*1e-9; % m

            try
                [tau1, msd1] = msd_2d(xy1,dt);
                [tau2, msd2] = msd_2d(xy2,dt);
            catch ME
                warning('MSD computation failed for trace %d in %s: %s', ti, spath, ME.message);
                continue
            end

            nFit = min(4,length(tau1)); nFit = max(nFit,2);
            p1 = polyfit(tau1(1:nFit), msd1(1:nFit), 1);
            slope1 = p1(1);
            nFit2 = min(4,length(tau2)); nFit2 = max(nFit2,2);
            p2 = polyfit(tau2(1:nFit2), msd2(1:nFit2), 1);
            slope2 = p2(1);

            D1 = slope1/4;
            D2 = slope2/4;

            % viscosity via Broersma correction
            eta1 = kB*T / (D1 * f_avg);
            eta2 = kB*T / (D2 * f_avg);

            D_all(end+1,1:2) = [D1, D2]; %#ok<SAGROW>
            eta_all(end+1,1:2) = [eta1, eta2]; %#ok<SAGROW>
        end
    end

    if isempty(D_all)
        results.(gname).D = []; results.(gname).eta = [];
        results.(gname).stats = [];
        continue
    end

    D_mean_per_trace = mean(D_all,2);
    eta_mean_per_trace = mean(eta_all,2);

    stats.D.mean = mean(D_mean_per_trace);
    stats.D.median = median(D_mean_per_trace);
    stats.D.std = std(D_mean_per_trace);

    stats.eta.mean = mean(eta_mean_per_trace);
    stats.eta.median = median(eta_mean_per_trace);
    stats.eta.std = std(eta_mean_per_trace);

    results.(gname).D = D_mean_per_trace;
    results.(gname).eta = eta_mean_per_trace;
    results.(gname).raw.D_channels = D_all;
    results.(gname).raw.eta_channels = eta_all;
    results.(gname).stats = stats;

    fprintf('  collected %d valid traces for %s\n', numel(eta_mean_per_trace), gname);
end

%% Save results
save(fullfile(outFolder,'translational_results.mat'),'results','dt','minPoints','T','a','b','p','f_par','f_perp','f_avg');

%% Prepare summary and boxplot
gNames = fieldnames(results);
nums = zeros(numel(gNames),1);
for i=1:numel(gNames)
    tokens = regexp(gNames{i},'glycerol_(\d+)','tokens');
    if ~isempty(tokens), nums(i)=str2double(tokens{1}{1}); else nums(i)=inf; end
end
[~,ord]=sort(nums);
gNames = gNames(ord);

allEta = {}; gLabels = {};
for i=1:numel(gNames)
    g = gNames{i};
    allEta{i} = results.(g).eta;
    gLabels{i} = g;
end

figure('visible','off');
groupedEta = []; groupIdx = [];
for i=1:numel(allEta)
    if isempty(allEta{i}), continue; end
    groupedEta = [groupedEta; allEta{i}]; %#ok<AGROW>
    groupIdx = [groupIdx; i*ones(numel(allEta{i}),1)]; %#ok<AGROW>
end
boxplot(groupedEta,groupIdx,'Labels',gLabels);
xlabel('Glycerol content folder'); ylabel('Viscosity (Pa·s)');
title('Viscosity per glycerol content (Broersma-corrected)');
hold on
means = zeros(numel(gNames),1);
for i=1:numel(gNames)
    if isempty(allEta{i}), means(i)=NaN; else means(i)=mean(allEta{i}); end
end
xticks = 1:numel(gNames);
plot(xticks,means,'--ok','LineWidth',1.2,'MarkerFaceColor','w');
hold off

pngFile = fullfile(outFolder,'viscosity_boxplot.png');
svgFile = fullfile(outFolder,'viscosity_boxplot.svg');
print(gcf,'-dpng','-r300',pngFile);
try
    exportgraphics(gcf,svgFile,'ContentType','vector');
catch
    saveas(gcf,svgFile);
end
close(gcf);

%% Save summary CSV
rows = {};
for i=1:numel(gNames)
    g = gNames{i};
    st = results.(g).stats;
    if isempty(st)
        rows{end+1,1} = g; rows{end,2:7} = {NaN,NaN,NaN,NaN,NaN,NaN}; %#ok<SAGROW>
    else
        rows{end+1,1} = g;
        rows{end,2} = st.D.mean; rows{end,3}=st.D.median; rows{end,4}=st.D.std;
        rows{end,5} = st.eta.mean; rows{end,6}=st.eta.median; rows{end,7}=st.eta.std; %#ok<SAGROW>
    end
end
Tsummary = cell2table(rows,'VariableNames',{'GlycerolFolder','D_mean','D_median','D_std','eta_mean','eta_median','eta_std'});
writetable(Tsummary,fullfile(outFolder,'translational_summary_stats.csv'));

fprintf('Done. Results saved in %s\n', outFolder);

%% Local helper
function [tau, msd] = msd_2d(xy,dt)
    N = size(xy,1);
    maxlag = N-1;
    msd = zeros(maxlag,1);
    tau = (1:maxlag)'.*dt;
    for lag = 1:maxlag
        diffs = xy(1+lag:end,:) - xy(1:end-lag,:);
        sq = sum(diffs.^2,2);
        msd(lag) = mean(sq);
    end
end
