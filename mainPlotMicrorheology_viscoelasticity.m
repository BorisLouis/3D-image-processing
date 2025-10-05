%% Script to compute storage and loss moduli from MSD (τ-domain)
% Author: ChatGPT
% Date: 2025-09-23

clear; clc;
close all

%% PARAMETERS
mainFolder = 'S:\Dual Color\20250121_dualcolor\PS_286g_136r';

aList  = [0.143, 0.068];           % bead radius [µm]
T  = 296.15;              % temperature [K]
kB = 1.380649e-23;     % Boltzmann constant [J/K]

fracKeep = 0.25;       % fraction of MSD points to keep (most reliable)
minPairs = 10;          % minimum number of displacements per MSD point
nInterp  = 50;         % #points in common τ-grid

%% GET TIMEPOINT FOLDERS
d = dir(mainFolder);
isDir = [d.isdir];
timeFolders = {d(isDir).name};
timeFolders = timeFolders(~ismember(timeFolders,{'.','..'}));
% keep only folders with "_min" in the name
timeFolders = timeFolders(contains(timeFolders,'_min'));
nTimepoints = numel(timeFolders);

results = struct();

%% LOOP OVER TIMEPOINTS AND CHANNELS
for t = 1:nTimepoints
    tFolder = fullfile(mainFolder, timeFolders{t});
    sDirs = dir(tFolder);
    sDirs = sDirs([sDirs.isdir]);
    sNames = {sDirs.name};
    sNames = sNames(~ismember(sNames,{'.','..'}));
    % keep only sample folders (contain time string + '_')
    sNames = sNames(contains(sNames,[timeFolders{t} '_']));
    
    for ch = 1:2
        a = aList(ch);
        allTau = {};
        allGp  = {};
        allGpp = {};
        
        for s = 1:numel(sNames)
            sFolder = fullfile(tFolder, sNames{s});
            fName = fullfile(sFolder, sprintf('msdRes%d.mat', ch));
            if ~isfile(fName), continue; end
            
            S = load(fName); % loads allRes
            if ~isfield(S,'allRes'), continue; end
            
            for tr = 1:numel(S.allRes)
                msdFull = S.allRes(tr).msdR(:);
                tauFull = S.allRes(tr).tau(:);
                if isempty(msdFull) || isempty(tauFull), continue; end
                
                % keep only reliable fraction of MSD
                nPts = numel(msdFull);
                nKeep = max(5, round(fracKeep*nPts));
                msd = msdFull(1:nKeep);
                tau = tauFull(1:nKeep);
                
                % exclude points with very poor statistics
                validMask = (nPts - (1:nKeep) + 1) >= minPairs;
                msd = msd(validMask);
                tau = tau(validMask);

                if numel(msd)<5, continue; end

                logTau = log10(tau);
                logMSD = log10(msd);
                windowLength = min(11, numel(tau)); % must be odd, <= number of points
                if mod(windowLength,2)==0
                    windowLength = windowLength-1; % make odd
                end
                polyOrder = 2;
                smoothLogMSD = sgolayfilt(logMSD, polyOrder, windowLength);
                alpha = gradient(smoothLogMSD) ./ gradient(logTau);
                % alpha(alpha < 0) = 0;
                % alpha(alpha > 1) = 1;


                % alpha = gradient(log(msd*(10^-15))) ./ gradient(log(tau));
                
                % complex modulus magnitude (Mason–Weitz GSER)
                Gmag = (kB*T) ./ (pi*a*10^(-6)*msd*10^-(15) .* gamma(1+alpha));
                
                % storage and loss moduli
                Gp  = Gmag .* cos(pi*alpha/2);
                Gpp = Gmag .* sin(pi*alpha/2);
                
                allTau{end+1} = tau;
                allGp{end+1}  = Gp;
                allGpp{end+1} = Gpp;
            end
        end
        
        % average across traces on a common τ-grid
        if ~isempty(allGp)
            minTau = max(cellfun(@(v) min(v), allTau));
            maxTau = min(cellfun(@(v) max(v), allTau));
            if minTau < maxTau
                tauGrid = logspace(log10(minTau), log10(maxTau), nInterp);
                
                GpMat  = nan(numel(allGp), nInterp);
                GppMat = nan(numel(allGpp), nInterp);
                for ii = 1:numel(allGp)
                    GpMat(ii,:)  = interp1(allTau{ii}, allGp{ii},  tauGrid, 'linear', NaN);
                    GppMat(ii,:) = interp1(allTau{ii}, allGpp{ii}, tauGrid, 'linear', NaN);
                end
                
                results(t,ch).tau = tauGrid;
                results(t,ch).Gp  = nanmean(GpMat,1);
                results(t,ch).Gpp = nanmean(GppMat,1);
            else
                results(t,ch).tau = [];
                results(t,ch).Gp  = [];
                results(t,ch).Gpp = [];
            end
        else
            results(t,ch).tau = [];
            results(t,ch).Gp  = [];
            results(t,ch).Gpp = [];
        end
    end
end

%% DETERMINE GLOBAL AXIS LIMITS
allTauVals = [];
allModVals = [];
for t = 1:nTimepoints
    for ch = 1:2
        allTauVals = [allTauVals, results(t,ch).tau];
        allModVals = [allModVals, results(t,ch).Gp, results(t,ch).Gpp];
    end
end
allTauVals = allTauVals(~isnan(allTauVals) & allTauVals>0);
allModVals = allModVals(~isnan(allModVals) & allModVals>0);

xRange = [min(allTauVals), max(allTauVals)];
yRange = [min(allModVals), max(allModVals)];

%% PLOT RESULTS
figure('Units','normalized','Position',[0.05 0.08 0.9 0.82]);
Order = [3 4 5 6 1 2];
nTimepoints = 6;
tselect = [2 3 4 5 6 7];
n = 0;
for t = tselect
    n = n+1;
    for ch = 1:2
        axIdx = (ch-1)*nTimepoints + Order(n);
        subplot(2, nTimepoints, axIdx);
        
        tau = results(t,ch).tau;
        Gp  = results(t,ch).Gp;
        Gpp = results(t,ch).Gpp;
        
        if isempty(tau)
            title(sprintf('%s - Ch%d (no data)', timeFolders{t}, ch));
            axis off;
            continue;
        end
        
        semilogy(tau, Gp, '-o', 'DisplayName', 'G'' (storage)'); hold on;
        semilogy(tau, Gpp, '-s', 'DisplayName', 'G'''' (loss)');
        
        xlabel('Lag time \tau [s]');
        ylabel('Modulus [Pa]');
        title(sprintf('%s - Ch%d', strrep(timeFolders{t}, '_', ' '), ch));
        legend('Location','best');
        grid on;
        ylim([0.1, 5000]);
        xlim([0 0.15]);
        % if ch == 1
        %     ylim(yRange);
        % elseif ch == 2
        %     ylim(yRange);
        % end
    end
end
sgtitle('Averaged Storage (G'') and Loss (G'''') moduli per timepoint and channel (τ-domain)');