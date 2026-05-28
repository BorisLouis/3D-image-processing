classdef MPTICSMovie < Core.MPMovie
    %MPTICSMOVIE Summary of this class goes here
    %   Detailed explanation goes here
    
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
                    for frame = 1:MaxFrame
                        waitbar(frame./MaxFrame, h, append('calc sacf frame ', num2str(frame), '/',...
                                num2str(MaxFrame), ' - Plane ', num2str(c)));
                        Frame = Movie(:,:,frame);
                        sacf = obj.SACF(Frame);
                        [List(frame), R2(frame)] = obj.fitSACF(sacf);
                    end
                    cleanedList = List;
                    cleanedList(R2 < 0.90) = [];
                    Results.wListRaw = List.*obj.info.PxSize;
                    Results.wList = cleanedList.*obj.info.PxSize;
                    Results.wAvg = median(cleanedList).*obj.info.PxSize;
                    Results.R2 = R2;
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
                    model_fun = @(p, x) 1 ./ (1 + x ./ p(1));
                    lb        = [0];
                    ub        = [Inf];
                    opts      = optimoptions('lsqcurvefit', 'Display', 'off');
                    threshold = 0.15;
                    min_pts   = obj.info.FitTACF;    % ACF decays in ~3 points, so min is low

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
                        % This is the key fix — the ACF decays in ~3 points, then
                        % the tail is pure noise that was being included before
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
                    
                        % --- lsqcurvefit ---
                        try
                            [p_fit, ~] = lsqcurvefit(model_fun, p0, Tau_fit, AutoCorr_fit, lb, ub, opts);
                        catch
                            R2_map(i,j) = NaN;
                            diff(i,j)   = NaN;
                            eta(i,j)    = NaN;
                            continue;
                        end
                        LifeTime = p_fit(1);
                    
                        % --- R² ---
                        r_pred = model_fun(p_fit, Tau_fit);
                        SS_res = sum((AutoCorr_fit - r_pred).^2);
                        SS_tot = sum((AutoCorr_fit - mean(AutoCorr_fit)).^2);
                        R2_map(i,j) = 1 - SS_res / SS_tot;
                    
                        % --- Derived quantities ---
                        omega_um = (0.61*obj.info.Wavelength)./obj.info.NA.*10^(-3);
                        D = omega_um ./ (4*LifeTime);
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
                    caxis([0 obj.info.LimitViscMap])
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

        function sacf = SACF(obj, Frame)

            I = double(Frame);
            [ny, nx] = size(I);
            dI = I - mean(I(:));
            F = fft2(dI);
            R = ifft2(abs(F).^2);
            R = fftshift(R);
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
            sacf(sacf < 0) = 0;
            r_map = sacf;
            [nx, ny] = size(sacf);
            cx = floor(nx/2) + 1;
            cy = floor(ny/2) + 1;
            
            [xi, eta] = meshgrid(1:ny, 1:nx);
            xi = xi - cy;     % shift to center at zero
            eta = eta - cx;
            
            rho = sqrt(xi.^2 + eta.^2);   % radial distance
            
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
            
            r_rad = zeros(size(bin_centers));
            N_bin = zeros(size(bin_centers));
            
            for k = 1:length(bin_centers)
                mask = (rho_vec >= edges(k)) & (rho_vec < edges(k+1));
                r_rad(k) = mean(r_vec(mask));
                N_bin(k) = sum(mask);
            end
            
            % Remove empty bins
            nonempty = N_bin > 0;
            rho_fit = bin_centers(nonempty);
            r_fit   = r_rad(nonempty);
            w       = sqrt(N_bin(nonempty));   % weights = sqrt(counts)
            
            % -------------------------------------------------------------
            % 4. Initial parameter estimates
            % -------------------------------------------------------------
            Ginf0 = mean(r_fit(end-3:end));             % asymptote guess
            G0_0  = r_fit(1) - Ginf0;                   % amplitude guess
            omega0_0 = rho_fit(find(r_fit <= Ginf0 + G0_0*exp(-1),1));
            if isempty(omega0_0)
                omega0_0 = 2;  % fallback guess
            end
            
            p0 = [G0_0, omega0_0, Ginf0];   % initial parameters
            
            % -------------------------------------------------------------
            % 5. Define model function
            % -------------------------------------------------------------
            model_fun = @(p, rho) p(1)*exp(-(rho.^2)/(p(2)^2)) + p(3);
            
            % -------------------------------------------------------------
            % 6. Perform nonlinear least squares fit (weighted)
            % -------------------------------------------------------------
            opts = optimoptions('lsqcurvefit','Display','off');
            
            lb = [-Inf, 0, -Inf];   % omega0 must be positive
            ub = [ Inf, Inf, Inf];
            
            [p_fit, resnorm] = lsqcurvefit(@(p, rho) w .* model_fun(p,rho), ...
                                p0, rho_fit, w .* r_fit, lb, ub, opts);
            
            G0     = p_fit(1);
            omega  = p_fit(2);
            Ginf   = p_fit(3);

            r_pred  = model_fun(p_fit, rho_fit);          % model predictions
            SS_res  = sum((r_fit - r_pred).^2);           % residual sum of squares
            SS_tot  = sum((r_fit - mean(r_fit)).^2);      % total sum of squares
            R2      = 1 - SS_res / SS_tot;
            
            if strcmp(obj.info.PlotSACFfit, 'on')
                % -------------------------------------------------------------
                % 7. Plot results
                % -------------------------------------------------------------
                Fig = figure; hold on;
                plot(rho_fit, r_fit, 'ko', 'MarkerFaceColor','k','DisplayName','Radial Average');
                rho_plot = linspace(0, max(rho_fit), 200);
                plot(rho_plot, model_fun(p_fit, rho_plot), 'r-', 'LineWidth',2, ...
                     'DisplayName','Fit');
                
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

