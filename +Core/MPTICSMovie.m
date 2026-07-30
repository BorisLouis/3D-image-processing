classdef MPTICSMovie < Core.MPMovie
    %MPTICSMOVIE Summary of this class goes here
    %   Detailed explanation goes here
    %
    %   UNIT CONVENTIONS (obj.info fields) -- enforced throughout this file:
    %       obj.info.Wavelength   : emission wavelength, in NANOMETERS (nm)
    %       obj.info.NA           : numerical aperture, dimensionless
    %       obj.info.PxSize       : pixel size, in MICROMETERS (um) per pixel
    %       obj.info.Radius       : tracer/probe radius, in NANOMETERS (nm)
    %       obj.info.Temperature  : in KELVIN (K)
    %       obj.info.ExpTime      : frame exposure / lag time, in SECONDS (s)
    %       obj.info.BeadDiameter : (OPTIONAL) calibration-bead apparent size,
    %                               in MICROMETERS (um), used only as an
    %                               approximate finite-bead-size correction
    %                               in calculateOmega (see comments there).
    %
    %   If your PxSize is stored in nm/pixel instead, convert it
    %   (PxSize_um = PxSize_nm / 1000) before assigning obj.info.PxSize, or
    %   every length derived from pixel coordinates (i.e. the omega
    %   calibration) will be off by a factor of 1000 relative to the
    %   theoretical omega_um used as a fallback in getDiffusionmap.

    properties
        AllFrames
        Omegas
        AutocorrMap
        DiffusionMap
        ViscosityMap
        Results
    end
    
    methods
        function obj = MPTICSMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
        end
        
        function LoadAllFrames(obj)
            h = waitbar(0, 'initializing');
            for c = 1:obj.calibrated{1, 1}.nPlanes
                waitbar(c./(obj.calibrated{1, 1}.nPlanes), h, append('Loading plane ', num2str(c), '/', num2str(obj.calibrated{1, 1}.nPlanes)));
                Path = append(obj.calibrated{1, 1}.mainPath, filesep, 'calibratedPlane', num2str(c), '.tif');
                if strcmp(obj.info.frame2Load, 'all')
                    obj.AllFrames{c, 1} = Load.Movie.tif.getFrame(Path, 1:obj.raw.maxFrame);
                else
                    Load.movie.tif.getFrame(obj.calibrated{1, 1}.mainPath, obj.info.frame2Load);
                end
            end
            close(h)
        end

        function calculateOmega(obj)
            % Computes the spatial PSF/SACF-derived omega (beam waist) per
            % plane from immobilized-bead calibration movies.
            %
            % FIXES APPLIED vs. previous version:
            %   - fitSACF now restricts the fit to a local window around the
            %     peak and excludes the noise-contaminated zero-lag bin, so
            %     the fitted width genuinely reflects the central peak and
            %     not far-field background/noise (see fitSACF for details).
            %   - Results.wAvg is now explicitly documented to be in
            %     MICROMETERS, matching the theoretical omega_um used
            %     downstream in getDiffusionmap, as long as obj.info.PxSize
            %     is supplied in um/pixel (see class-level comment block).
            %   - Optional finite-bead-size correction: if you set
            %     obj.info.BeadDiameter (in um), this attempts a simple
            %     quadrature-subtraction correction
            %         omega_PSF = sqrt(max(omega_measured^2 - BeadDiameter^2, 0))
            %     This treats the bead's apparent profile as approximately
            %     Gaussian with a width equal to BeadDiameter, which is a
            %     simplification -- verify it suits your bead geometry
            %     before trusting it blindly. If BeadDiameter is not set,
            %     no correction is applied and wAvgCorrected == wAvg.
            %
            %   TEMPORAL-FLUCTUATION METHOD (for crowded/complex samples,
            %   e.g. cells with many overlapping particles on top of
            %   static cell architecture):
            %     Computing the SACF directly on a single raw frame mixes
            %     the PSF-scale correlation you want with much larger
            %     correlation lengths from static structure (cell body,
            %     nucleus boundary, slow gradients). A purely spatial
            %     background subtraction (e.g. imopen) cannot cleanly
            %     separate these in real cell images, because there's
            %     often no single length scale dividing "particle" from
            %     "cell structure" -- it's a continuum, not two discrete
            %     scales.
            %
            %     What DOES cleanly separate them is TIME: the cell
            %     architecture is essentially static frame-to-frame,
            %     while the particles you care about are exactly the
            %     thing that's fluctuating. So for every plane:
            %       1. Compute the time-averaged image over the
            %          calibration frames (the static structure).
            %       2. For every frame, subtract that time-average to
            %          get the pure fluctuation signal (particles only,
            %          ~zero-mean) -- this removes static structure
            %          regardless of its spatial scale.
            %       3. Compute the SACF of each fluctuation frame (with
            %          SACF's internal imopen step skipped -- it isn't
            %          meaningful on an already ~zero-mean, signed map),
            %          and AVERAGE the SACF maps themselves across
            %          frames. Averaging the (low-noise) SACF maps before
            %          fitting is far more stable than fitting many noisy
            %          per-frame SACFs and averaging the resulting omegas
            %          afterwards.
            %       4. Fit that single averaged SACF once with fitSACF
            %          -- this is Results.wAvg, the recommended number
            %          to use downstream in getDiffusionmap.
            %     Per-frame fits (List/wList/R2) are still computed too,
            %     but purely as a diagnostic/QC measure of frame-to-frame
            %     spread -- they are no longer the primary estimate.

            if strcmp(obj.info.runMethod, 'run')
                run = 1;
            else
                if exist(append(obj.raw.movInfo.Path, filesep, 'OmegaResults.mat'))
                    run = 0;
                else
                    run = 1;
                end
            end

            if run == 1
                h = waitbar(0, 'initializing');
                for c = 1:obj.calibrated{1, 1}.nPlanes
                    Movie = obj.AllFrames{c,1};
                    MaxFrame = obj.info.SACFframes; %size(obj.AllFrames, 3)

                    nanMask = mean(Movie, 3);
                    nanMask(nanMask == 0) = NaN;
                    nanIdx = isnan(nanMask);
                    
                    Movie(repmat(nanIdx, [1 1 size(Movie,3)])) = NaN;

                    % --- Temporal baseline: the static, non-fluctuating ---
                    % structure (cell body, nucleus, slow gradients),
                    % computed over the same frames analyzed below.
                    TimeBaseline = mean(Movie(:,:,1:MaxFrame), 3, 'omitnan');

                    List = [];   % reset per plane -- avoids stale entries
                    R2   = [];   % leaking across planes when nPlanes > 1
                    AccSACF = [];   % running sum for the averaged-SACF fit
                    nAcc = 0;

                    for frame = 1:MaxFrame
                        waitbar(frame./MaxFrame, h, append('calc sacf frame ', num2str(frame), '/',...
                                num2str(MaxFrame), ' - Plane ', num2str(c)));
                        Frame = double(Movie(:,:,frame));

                        % Pure fluctuation signal: removes static cell/
                        % nucleus architecture regardless of its spatial
                        % scale, leaving (approximately) only the moving/
                        % fluctuating particles.
                        FluctFrame = Frame - TimeBaseline;

                        % skipBackgroundSubtraction = true: the temporal
                        % step above already removed the static baseline,
                        % so SACF's internal imopen step is skipped (not
                        % meaningful on an already ~zero-mean, signed map).
                        sacf = obj.SACF(FluctFrame, true);

                        % Per-frame fit kept only for diagnostics/QC
                        % (frame-to-frame spread) -- not the primary omega.
                        [List(frame), R2(frame)] = obj.fitSACF(sacf);

                        % Accumulate for the averaged-SACF fit (primary method)
                        if isempty(AccSACF)
                            AccSACF = zeros(size(sacf));
                        end
                        sacfNoNan = sacf;
                        sacfNoNan(isnan(sacfNoNan)) = 0;
                        AccSACF = AccSACF + sacfNoNan;
                        nAcc = nAcc + 1;
                    end

                    AvgSACF = AccSACF ./ max(nAcc, 1);

                    % --- Primary estimate: one fit on the time-averaged SACF ---
                    [omegaAvg_px, R2Avg] = obj.fitSACF(AvgSACF);

                    cleanedList = List;
                    cleanedList(R2 < 0.90) = [];

                    % Diagnostic per-frame values, in PIXELS (see fitSACF).
                    % Convert to physical units (um) using the pixel size.
                    % Correct ONLY if obj.info.PxSize is in um/pixel -- see
                    % the class-level unit-convention comment block.
                    Results.wListRaw = List.*obj.info.PxSize;
                    Results.wList = cleanedList.*obj.info.PxSize;
                    Results.R2 = R2;          % per-frame R2, diagnostic only

                    % Primary result: from the averaged-SACF fit (more
                    % robust than median(cleanedList) for crowded/complex
                    % samples -- see header comment).
                    Results.wAvg  = omegaAvg_px .* obj.info.PxSize;
                    Results.R2Avg = R2Avg;

                    % --- Optional finite-probe-size correction (see header) ---
                    if isfield(obj.info, 'BeadDiameter') && ~isempty(obj.info.BeadDiameter) && obj.info.BeadDiameter > 0
                        Results.wAvgCorrected = sqrt(max(Results.wAvg.^2 - obj.info.BeadDiameter.^2, 0));
                        Results.wListCorrected = sqrt(max(Results.wList.^2 - obj.info.BeadDiameter.^2, 0));
                    else
                        Results.wAvgCorrected = Results.wAvg;
                        Results.wListCorrected = Results.wList;
                    end

                    obj.Omegas{c,1} = Results;
                end
                close(h)
    
                OmegaResults = obj.Omegas;
                FilePath = append(obj.raw.movInfo.Path, filesep, 'OmegaResults.mat');
                save(FilePath, "OmegaResults")
            else 
                load(append(obj.raw.movInfo.Path, filesep, 'OmegaResults.mat'));
                obj.Omegas = OmegaResults;
                disp('Found omegas - loaded that one')
            end
        end

        function getAutocorrmap(obj)
            if strcmp(obj.info.runMethod, 'run')
                run = 1;
            else
                if exist(append(obj.raw.movInfo.Path, filesep, 'AutocorrMap.mat'))
                    run = 0;
                else
                    run = 1;
                end
            end

            if run == 1
                h = waitbar(0, 'initializing');
                for c = 1:obj.calibrated{1, 1}.nPlanes
                    [rows, cols, nLags] = size(obj.AllFrames{c, 1});
                    autoCorr = zeros(rows, cols, nLags, 'single');
                    TACS_Matrix = nan(rows, cols, nLags, 'single');
                    frames = obj.AllFrames{c,1};
    
                    nanIm = mean(frames, 3);
                    nanIm(nanIm == 0) = NaN;
                    idx = ~isnan(nanIm);
                    
                    % Get linear indices of valid pixels
                    validPix = find(idx);       % vector of linear indices
                    nValid = numel(validPix);
    
                    for k = 1:nValid
                        lin = validPix(k);
                        [i, j] = ind2sub([rows, cols], lin);
                        waitbar(k / nValid, h, append('computing TACF plane ', num2str(c))) ;
                        ts = double(squeeze(frames(i,j,:)));
                        ts = ts - mean(ts, 'omitnan');
                        tsTrend = medfilt1(ts, 150);
                        ts = ts - tsTrend;
                        [ac] = obj.TACF(ts);
                        TACS_Matrix(i, j,:) = ac;
                    end
                    obj.AutocorrMap{c,1} = TACS_Matrix;
                end
                close(h)
    
                AutocorrMap = obj.AutocorrMap;
                FilePath = append(obj.raw.movInfo.Path, filesep, 'AutocorrMap.mat');
                save(FilePath, "AutocorrMap")
            else
                load(append(obj.raw.movInfo.Path, filesep, 'AutocorrMap.mat'));
                obj.AutocorrMap = AutocorrMap;
                disp('Found Autocorrelation map - loaded that one')
            end           
        end

        function getDiffusionmap(obj)
            % FIXES APPLIED vs. previous version:
            %   1. D = omega^2 / (4*tau_D), not D = omega / (4*tau_D).
            %      The previous code was missing the square on omega, which
            %      is required by the standard ICS/FCS relation
            %      tau_D = omega_0^2 / (4D) for 2D diffusion through a
            %      Gaussian observation area of waist omega_0.
            %   2. The per-plane CALIBRATED omega (from calculateOmega /
            %      obj.Omegas) is now actually used if available, instead
            %      of always falling back to the theoretical
            %      diffraction-limit formula. The theoretical formula is
            %      kept only as a fallback if no calibration was run for
            %      that plane.

            if strcmp(obj.info.runMethod, 'run')
                run = 1;
            else
                if exist(append(obj.raw.movInfo.Path, filesep, 'ViscosityMap.mat'))
                    run = 0;
                else
                    run = 1;
                end
            end

            if run == 1
                j = jet(256);
                ramp = linspace(0,1,256)';   % fades from black → jet
                BlackJet = j .* ramp;
    
                hh = waitbar(0, 'initializing');
                Tau = [0, (1:size(obj.AutocorrMap{1,1},3)).*obj.info.ExpTime]';
                Tau(end) = [];

                for c = 1:obj.calibrated{1, 1}.nPlanes              
                    data = obj.AutocorrMap{c,1};
                    blockSize = obj.info.TICSWindow;
    
                    newX = floor(size(data,1)/blockSize);
                    newY = floor(size(data,2)/blockSize); 
    
                    B = mean(reshape(data(1:newX*blockSize, 1:newY*blockSize, :), ...
                                     blockSize, newX, blockSize, newY, size(data,3)), [1 3]);
                    data = squeeze(B);
    
                    n = 0;
                    rows = size(data, 1);
                    cols = size(data, 2);
                    options = optimoptions('lsqcurvefit','Display','off');
                    eta = nan(rows, cols);
                    diff = nan(rows, cols);
    
                    nanIm = mean(data, 3);
                    nanIm(nanIm == 0) = NaN;
                    idx = ~isnan(nanIm);
                    validPix = find(idx);       % vector of linear indices
                    nValid = numel(validPix);
    
                    % --- Before loop ---
                    Tau = double(Tau(:));
                    model_fun = @(p, x) p(1) + (1 - p(1)) * (1 ./ (1 + (x ./ p(2))));
                    opts = optimoptions('lsqcurvefit', 'Display', 'off');
                    threshold = 0.15;
                    min_pts = obj.info.FitTACF;    % ACF decays in ~3 points, so min is low

                    % --- Resolve omega for this plane: calibrated first, theoretical fallback ---
                    if ~isempty(obj.Omegas) && size(obj.Omegas, 1) >= c && ~isempty(obj.Omegas{c,1}) ...
                            && isfield(obj.Omegas{c,1}, 'wAvgCorrected') && ~isnan(obj.Omegas{c,1}.wAvgCorrected) ...
                            && obj.Omegas{c,1}.wAvgCorrected > 0
                        omega_um = obj.Omegas{c,1}.wAvgCorrected/1000;   % already in um, see calculateOmega
                        omegaSource = 'calibrated (SACF-derived)';
                    else
                        omega_um = (0.61*obj.info.Wavelength)./obj.info.NA.*10^(-3); % nm -> um
                        omegaSource = 'theoretical (0.61*lambda/NA)';
                    end
                    disp(append('Plane ', num2str(c), ': using ', omegaSource, ...
                        ' omega = ', num2str(omega_um), ' um'));

                    for k = 1:nValid
                        lin = validPix(k);
                        [i, j] = ind2sub([rows, cols], lin);
                        waitbar(k./nValid, hh, append('Fitting on TACF ', num2str(k), '/', ...
                            num2str(nValid), ' - plane ', num2str(c)));
                    
                        % --- AutoCorr ---
                        AutoCorr = double(squeeze(data(i,j,:)));
                        ac_zero = AutoCorr(1);                         % lag-0 value
                        if ac_zero <= 0 || ~isfinite(ac_zero)
                            R2_map(i,j) = NaN;
                            diff(i,j)   = NaN;
                            eta(i,j)    = NaN;
                            continue;
                        end
                        AutoCorr = AutoCorr ./ ac_zero;                % normalize by lag-0
                        AutoCorr = AutoCorr(:);
                    
                        % --- Adaptive cutoff: FIRST crossing, not last ---
                        cutoff_idx = find(AutoCorr < threshold, 1, 'first');
                        if isempty(cutoff_idx)
                            cutoff_idx = min_pts;
                        end
                        cutoff_idx = max(cutoff_idx, min_pts);  % enforce minimum
                    
                        Tau_fit      = Tau(1:cutoff_idx);
                        AutoCorr_fit = AutoCorr(1:cutoff_idx);
                    
                        % --- Per-pixel p0 from 1/e crossing ---
                        tau_guess_idx = find(AutoCorr < exp(-1), 1, 'first');
                        if isempty(tau_guess_idx) || tau_guess_idx == 1
                            p0 = Tau(min(2, end));
                        else
                            p0 = Tau(tau_guess_idx);
                        end
                        p0 = max(p0, 1e-10);
                        if ~isfinite(p0)
                            p0 = mean(Tau_fit);
                        end
                        p0_all = [0, p0];
                        lb = [0, 0];
                        ub = [max(AutoCorr_fit), Tau_fit(find(AutoCorr_fit < 0, 1, 'first'))];
                        
                        % --- lsqcurvefit ---
                        try
                            [p_fit, ~] = lsqcurvefit(model_fun, p0_all, Tau_fit, AutoCorr_fit, lb, ub, opts);
                        catch
                            R2_map(i,j) = NaN;
                            diff(i,j) = NaN;
                            eta(i,j) = NaN;
                            continue;
                        end
                        LifeTime = p_fit(2);
                        
                        % Fig = figure; hold on;
                        % plot(Tau_fit, AutoCorr_fit,'o', 'Color', [0.7 0.7 0.7], 'DisplayName','Radial Average (full)');
                        % plot(Tau_fit, model_fun(p_fit, Tau_fit), 'r-', 'LineWidth',2, ...
                        %      'DisplayName','Fit');
                        % xlabel('timelag (ms)');
                        % ylabel('r(\rho)');
                        % legend();
                        % grid on;
                    
                        % --- R² ---
                        r_pred = model_fun(p_fit, Tau_fit);
                        SS_res = sum((AutoCorr_fit - r_pred).^2);
                        SS_tot = sum((AutoCorr_fit - mean(AutoCorr_fit)).^2);
                        R2_map(i,j) = 1 - SS_res / SS_tot;
                        LifeTimeMap(i,j) = LifeTime;
                    
                        % --- Derived quantities ---
                        % FIX: D = omega_0^2 / (4*tau_D)  (omega must be squared)
                        D = (omega_um.^2) ./ (4*LifeTime);
                        eta(i,j) = (1.380649e-23 * obj.info.Temperature) / ...
                                   (6*pi * obj.info.Radius*1e-9 * D*1e-12) * 1e3;
                        diff(i,j) = D;
                    end

                    etaRes = eta;
                    diffRes = diff;
                    try
                        etaRes = imresize(eta, [obj.raw.movInfo.Width, obj.raw.movInfo.Length]);
                        diffRes = imresize(diff, [obj.raw.movInfo.Width, obj.raw.movInfo.Length]);
                    catch
                        etaRes = eta;
                        diffRes = diff;
                    end

                    Fig1 = figure(); 
                    imagesc(etaRes)
                    colormap(BlackJet);          % <-- apply colormap here
                    cb = colorbar;               % <-- no arguments here
                    cb.Label.String = 'Viscosity (cP)';
                    caxis([0 5]); %obj.info.LimitViscMap
                    title(append('Viscosity map - av visc ', num2str(median(etaRes, 'all', 'omitnan')), ' +/- ', num2str(std(etaRes(:), 'omitnan')), ' cP'));
                    Fig1Path = append(obj.raw.movInfo.Path, filesep, 'ViscosityMap_Plane', num2str(c), '.png');
                    saveas(Fig1, Fig1Path);
    
                    Fig2 = figure();
                    imagesc(diffRes)
                    set(gca, 'ColorScale', 'log');
                    colormap(BlackJet);          % <-- apply colormap here
                    cb = colorbar; 
                    cb.Label.String = 'Diffusion coefficient (µm^2/s)';
                    title(append('Diffusion map - av diffusion ', num2str(median(diffRes, 'all', 'omitnan')), ' +/- ', num2str(std(diffRes(:), 'omitnan')), ' µm^2/s'));
                    Fig2Path = append(obj.raw.movInfo.Path, filesep, 'DiffusionMap_Plane', num2str(c), '.png');
                    saveas(Fig2, Fig2Path);

                    Fig3 = figure();
                    imagesc(R2_map)
                    clim([0 1])
                    cb = colorbar;
                    cb.Label.String = 'R^2 fit';
                    cb.Limits = [0 1];
                    title(append('Fit error map'));
                    Fig3Path = append(obj.raw.movInfo.Path, filesep, 'FitErrorMap_Plane', num2str(c), '.png');
                    saveas(Fig3, Fig3Path);

                    Results.ViscMean = nanmean(etaRes, 'all');
                    Results.ViscStd = nanstd(etaRes(:));
                    Results.DiffMean = nanmean(diffRes, 'all');
                    Results.DiffStd = nanstd(diffRes(:));
                    Results.OmegaUsed_um = omega_um;
                    Results.OmegaSource = omegaSource;
    
                    obj.ViscosityMap{c,1} = etaRes;
                    obj.DiffusionMap{c,1} = diffRes;
                    obj.Results{c,1} = Results;
                    close all
                end
                close(hh)
    
                ViscMap = obj.ViscosityMap;
                DiffMap = obj.DiffusionMap;
                TICSResults = obj.Results;
    
                FilePathVisc = append(obj.raw.movInfo.Path, filesep, 'ViscosityMap.mat');
                FilePathDiff = append(obj.raw.movInfo.Path, filesep, 'DiffusionMap.mat');
                FilePathRes = append(obj.raw.movInfo.Path, filesep, 'TICSResults.mat');
                save(FilePathVisc, "ViscMap");
                save(FilePathDiff, "DiffMap");
                save(FilePathRes, "TICSResults");
            else
                load(append(obj.raw.movInfo.Path, filesep, 'ViscosityMap.mat'));
                load(append(obj.raw.movInfo.Path, filesep, 'DiffusionMap.mat'));
                load(append(obj.raw.movInfo.Path, filesep, 'TICSResults.mat'));
                obj.ViscosityMap = ViscMap;
                obj.DiffusionMap = DiffMap;
                obj.Results = TICSResults;
                disp('Found Viscosity map - loaded that one')
                disp('Found Diffusion map - loaded that one')
                disp('Found TICS results - loaded that one')
            end
        end

        function sacf = SACF(obj, Frame, skipBackgroundSubtraction)
            % FIX vs. previous version: the autocorrelation is now computed
            % via zero-padded FFTs with explicit overlap-count normalization,
            % instead of a plain circular (wrap-around) FFT autocorrelation.
            % Without this, opposite edges of the image are implicitly
            % treated as adjacent, which biases the autocorrelation
            % (especially the offset/far-field level that fitSACF's Ginf
            % term tries to capture) and can in turn bias the fitted omega.
            %
            % Output size/format is unchanged: sacf is still a 2D map the
            % same size as Frame, centered the same way, so fitSACF needs
            % no changes to consume it.
            %
            % skipBackgroundSubtraction (optional, default false): set to
            % true when Frame has ALREADY had a slowly-varying baseline
            % removed upstream (e.g. the temporal-fluctuation subtraction
            % in calculateOmega). In that case Frame is an already
            % ~zero-mean, signed map, and running the morphological-
            % opening background subtraction below on it is not
            % meaningful (imopen assumes locally-bright-structure-over-
            % darker-background, which doesn't hold once the mean has
            % already been removed and values can go negative) -- so it
            % is skipped entirely in that case.

            if nargin < 3 || isempty(skipBackgroundSubtraction)
                skipBackgroundSubtraction = false;
            end

            I = double(Frame);
            [ny, nx] = size(I);
            if skipBackgroundSubtraction
                dI = I;
            else
                dI = I - imopen(I, strel('disk', 10));
            end

            % Zero-pad to 2x size in each dimension to get the correct
            % linear (non-circular) autocorrelation, not the circular one.
            dIpad = zeros(2*ny, 2*nx);
            dIpad(1:ny, 1:nx) = dI;

            F = fft2(dIpad);
            R = ifft2(abs(F).^2);
            R = real(fftshift(R));

            % Overlap-count map: how many real (non-padded) pixel pairs
            % actually contributed at each lag. Used to correct for the
            % shrinking number of valid pairs near the edges of the
            % padded/shifted result (equivalent to dividing by N(lag)
            % instead of always dividing by the same N).
            onesMask = zeros(2*ny, 2*nx);
            onesMask(1:ny, 1:nx) = 1;
            Fo = fft2(onesMask);
            overlap = real(fftshift(ifft2(abs(Fo).^2)));
            overlap(overlap < 1) = NaN; % guard against division by ~0 outside valid support

            R = R ./ overlap;

            % Crop back to the original frame size, centered on the
            % zero-lag point, so the returned sacf has the same size/shape
            % as before (fitSACF centers itself automatically based on
            % size(sacf), so this crop is purely for consistency/efficiency).
            cy0 = ny + 1; cx0 = nx + 1; % center index after padding+fftshift
            rowsIdx = (cy0 - floor(ny/2)) : (cy0 + ceil(ny/2) - 1);
            colsIdx = (cx0 - floor(nx/2)) : (cx0 + ceil(nx/2) - 1);
            R = R(rowsIdx, colsIdx);

            mu = mean(I(:));
            sacf = R / (mu * mu);
        
        end

        function [tacf] = TACF(obj, IntTrace)
            x = IntTrace;
            T = length(x);
            
            % --- Step 1: fluctuations (mean already subtracted upstream) ---
            mu = mean(x);
            dx = x - mu;
            
            % --- Normalize by variance at lag=0 (guarantees G(0)=1) ---
            denominator = sum(dx .* dx) / T;   % variance = G(0) numerator
            
            if denominator == 0 || ~isfinite(denominator)
                tacf = zeros(1, T);
                return;
            end
            
            maxTau = T - 1;
            tacf   = zeros(1, maxTau + 1);
            
            for tau = 0:maxTau
                Nt = T - tau;
                tacf(tau+1) = (sum(dx(1:Nt) .* dx(1+tau:Nt+tau)) / Nt) / denominator;
            end
        end

        function [omega, R2] = fitSACF(obj, sacf)
            % FIXES APPLIED vs. previous version:
            %   1. The fit is now restricted to a local radial window around
            %      the peak (a small multiple of the initial omega guess)
            %      instead of extending all the way to r_max. Previously,
            %      the weight w = sqrt(N_bin) grows with radius (more
            %      pixels per annulus at larger rho), so the far field --
            %      which carries no information about the peak width, only
            %      about the background offset and any residual structure
            %      / wrap-around artifacts -- could dominate the cost
            %      function and pull the fitted omega upward.
            %   2. The zero-lag bin (rho ~ 0) is excluded from the fit
            %      window, since it can carry an extra spike from detector
            %      shot/read noise that doesn't belong to the smooth
            %      diffraction-limited shape being fitted.
            %   3. omega is still returned in PIXELS (unchanged output
            %      convention) -- conversion to physical units happens in
            %      calculateOmega via obj.info.PxSize (must be um/pixel,
            %      see class-level comment block).

            sacf(sacf < 0) = 0;
            r_map = sacf;
            [nx, ny] = size(sacf);
            cx = floor(nx/2) + 1;
            cy = floor(ny/2) + 1;
            
            [xi, eta] = meshgrid(1:ny, 1:nx);
            xi = xi - cy;     % shift to center at zero
            eta = eta - cx;
            
            rho = sqrt(xi.^2 + eta.^2);   % radial distance, in PIXELS
            
            % -------------------------------------------------------------
            % 2. Flatten arrays for easier processing
            % -------------------------------------------------------------
            rho_vec = rho(:);
            r_vec   = r_map(:);
            
            % remove NaNs if present
            valid = ~isnan(r_vec);
            rho_vec = rho_vec(valid);
            r_vec   = r_vec(valid);
            
            % -------------------------------------------------------------
            % 3. Define radial bins
            % -------------------------------------------------------------
            dr = 1;   % bin width in pixels
            r_max = max(rho_vec);
            edges = 0:dr:r_max;
            bin_centers = edges(1:end-1) + dr/2;
            
            r_rad = zeros(size(edges));
            N_bin = zeros(size(edges));
            
            for k = 1:length(bin_centers)
                mask = (rho_vec >= edges(k)) & (rho_vec < edges(k+1));
                r_rad(k) = mean(r_vec(mask));
                N_bin(k) = sum(mask);
            end
            N_bin = [N_bin, 0];
            % Remove empty bins
            nonempty = N_bin > 0;
            rho_fit = edges(nonempty);
            r_fit   = r_rad(nonempty);
            w       = sqrt(N_bin(nonempty));   % weights = sqrt(counts)
            % -------------------------------------------------------------
            % 4. Initial parameter estimates (computed on the FULL radial
            %    range, so the far-field offset/baseline guess Ginf0 is
            %    still robust even though the actual fit below is windowed)
            % -------------------------------------------------------------
            Ginf0 = mean(r_fit(end-3:end));             % asymptote guess
            G0_0  = r_fit(1) - Ginf0;                   % amplitude guess
            omega0_0 = rho_fit(find(r_fit <= Ginf0 + G0_0*exp(-1),1));
            if isempty(omega0_0)
                omega0_0 = 2;  % fallback guess
            end
            
            p0 = [G0_0, omega0_0, Ginf0];   % initial parameters
            
            % -------------------------------------------------------------
            % 4b. Restrict the actual fit to a local window around the
            %     peak: a few times the initial omega guess, excluding the
            %     zero-lag bin. This keeps the far field (noise/background/
            %     residual wrap-around) from dominating the cost function.
            % -------------------------------------------------------------
            fitWindow = max(6*omega0_0, 5);   % at least 5 px, else 6x initial guess
            fitWindow = min(fitWindow, 10);   % empirical cap: in complex/
                                               % crowded samples, residual
                                               % slow structure can still
                                               % contaminate the fit beyond
                                               % ~10 px even after temporal-
                                               % fluctuation preprocessing
                                               % (see calculateOmega), so
                                               % don't let the window grow
                                               % past this regardless of
                                               % the initial omega guess.
            fitWindow = min(fitWindow, r_max); % never exceed available data

            keepMask = (rho_fit > 0.5) & (rho_fit <= fitWindow);
            if nnz(keepMask) < 5
                % fall back to using everything if the window left too few points
                keepMask = true(size(rho_fit));
            end

            rho_fit_win = rho_fit(keepMask);
            r_fit_win   = r_fit(keepMask);
            w_win       = w(keepMask);

            % -------------------------------------------------------------
            % 5. Define model function
            % -------------------------------------------------------------
            model_fun = @(p, rho) p(1)*exp(-(rho.^2)/(p(2)^2)) + p(3);

            % -------------------------------------------------------------
            % 6. Perform nonlinear least squares fit (weighted, windowed)
            % -------------------------------------------------------------
            opts = optimoptions('lsqcurvefit','Display','off');

            lb = [0, 0, 0];   % omega0 must be positive
            ub = [Inf, rho_fit_win(end), r_fit(2*size(r_fit_win, 2))];
            
            [p_fit, resnorm] = lsqcurvefit(@(p, rho) w_win .* model_fun(p, rho), ...
                    p0, rho_fit_win, w_win .* r_fit_win, lb, ub, opts);
            
            omega  = p_fit(2);
            Ginf   = p_fit(3);
            G0     = p_fit(1);

            r_pred  = model_fun(p_fit, rho_fit_win);      % model predictions
            SS_res  = sum((r_fit_win - r_pred).^2);       % residual sum of squares
            SS_tot  = sum((r_fit_win - mean(r_fit_win)).^2); % total sum of squares
            R2      = 1 - SS_res / SS_tot;
            
            if strcmp(obj.info.PlotSACFfit, 'on')
                % -------------------------------------------------------------
                % 7. Plot results (full radial range shown for context, fit
                %    window highlighted)
                % -------------------------------------------------------------
                Fig = figure; hold on;
                plot(rho_fit, r_fit, 'o', 'Color', [0.7 0.7 0.7], 'DisplayName','Radial Average (full)');
                plot(rho_fit_win, r_fit_win, 'ko', 'MarkerFaceColor','k','DisplayName','Used in fit');
                rho_plot = linspace(0, max(rho_fit), 200);
                plot(rho_plot, model_fun(p_fit, rho_plot), 'r-', 'LineWidth',2, ...
                     'DisplayName','Fit');
                xline(fitWindow, '--', 'fit window', 'DisplayName', 'fit window cutoff');
                
                xlabel('\rho (pixels)');
                ylabel('r(\rho)');
                legend();
                title(sprintf('Fit Result:  \\omega_0 = %.3f pixels', omega));
                grid on;

                saveas(Fig, append(obj.raw.movInfo.Path, filesep, 'SacfPlot.png'))
            end
        end

        
    end
end