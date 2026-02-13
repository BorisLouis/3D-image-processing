classdef MPDDMMovie < Core.MPMovie
    %MPDDMMOVIE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        AllFrames
        DDMOutput
        FitResults
        MSDResults
        IntResults
        SizeResults
    end
    
    methods
        function obj = MPDDMMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
        end
        
        function getFullFrames(obj)
            for c = 1:obj.calibrated{1, 1}.nPlanes
                h = waitbar(0, 'initializing');
                if strcmp(obj.info.frame2Load, 'all')
                    FramesToLoad = 1:obj.raw.movInfo.maxFrame;
                    MaxFrame = obj.raw.movInfo.maxFrame;
                else
                    FramesToLoad = obj.info.frame2Load;
                    MaxFrame = numel(FramesToLoad);
                end
                n = 0;
                for frame = FramesToLoad
                    n = n+1;
                    FilePath = obj.calibrated{1, 1}.filePath.(append('plane', num2str(c)));
                    FilePath(1) = obj.raw.fullPath(1);
                    Frame = double(Load.Movie.tif.getFrame(FilePath, frame));
                    if strcmp(obj.info.ddmParam.CorrectBleaching, 'on')
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame + bleaching correction -- frame ',...
                            num2str(n), '/', num2str(MaxFrame), ' plane ', num2str(c)));
                        if frame == 1
                            RefInt = mean(Frame(Frame ~= 0));
                        else
                            Frame = Frame * (RefInt ./ mean(Frame(Frame ~= 0)));
                        end
                    else
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame - no bleaching correction -- frame ',...
                            num2str(n), '/', num2str(MaxFrame), ' plane ', num2str(c)));
                    end
                    obj.AllFrames(:,:,n,c) = Frame;

                    FrameList = Frame(:);
                    FrameList(FrameList == 0) = [];
                    obj.IntResults.MeanInt = mean(FrameList);
                    obj.IntResults.VarInt = var(FrameList);
                    obj.SizeResults = numel(FrameList).*(obj.info.pxSize*10^(-3)).^2;
                    
                end
                close(h)
            end
        end

        function [Idx] = getROI(obj)
            if strcmp(obj.info.frame2Load, 'all')
                FramesToLoad = 1:obj.raw.movInfo.maxFrame;
            else
                FramesToLoad = obj.info.frame2Load;
            end

            frame = FramesToLoad(1);
            Frame = double(Load.Movie.tif.getframes(obj.calibrated{1, 1}.filePath.(append('plane', num2str(1))), frame));

            Idx = find(Frame ~= 0);
        end

        function getFrameParts(obj, Px)

            for c = 1:obj.calibrated{1, 1}.nPlanes
                h = waitbar(0, 'initializing');
                if strcmp(obj.info.frame2Load, 'all')
                    FramesToLoad = 1:obj.raw.movInfo.maxFrame;
                    MaxFrame = obj.raw.movInfo.maxFrame;
                else
                    FramesToLoad = obj.info.frame2Load;
                    MaxFrame = numel(FramesToLoad);
                end
                n = 0;
                for frame = FramesToLoad
                    n = n+1;
                    Frame = double(Load.Movie.tif.getframes(obj.calibrated{1, 1}.filePath.(append('plane', num2str(c))), frame));
                    ROISize = obj.info.ddmParam.ROISize;
                    [PxRow,PxCol] = ind2sub(size(Frame),Px);
                    ROI = [PxRow - round(ROISize./2), PxRow + round(ROISize./2);...
                        PxCol - round(ROISize./2), PxCol + round(ROISize./2)];
                    Frame = Frame(ROI(1,1):ROI(1,2), ROI(2,1):ROI(2,2));

                    if strcmp(obj.info.ddmParam.CorrectBleaching, 'on')
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame + bleaching correction -- frame ',...
                            num2str(n), '/', num2str(MaxFrame), ' plane ', num2str(c)));
                        if frame == 1
                            RefInt = mean(Frame(Frame ~= 0));
                        else
                            Frame = Frame * (RefInt ./ mean(Frame(Frame ~= 0)));
                        end
                    else
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame - no bleaching correction -- frame ',...
                            num2str(n), '/', num2str(MaxFrame), ' plane ', num2str(c)));
                    end
                    obj.AllFrames(:,:,n,c) = Frame;
                end
                close(h)
            end
        end

        function mainDDM(obj,varargin)
            for c = obj.calibrated{1,1}.nPlanes
                DDMOutputfile = append(obj.calibrated{1, 1}.mainPath, filesep, 'DDMOutput' , num2str(c), '.mat');
                DDMOutputfile(1) = obj.raw.fullPath(1);
                DDMOutputfile = append(obj.raw.movInfo.Path,filesep, 'calibrated1', filesep, 'DDMOutput1.mat');
                if ~exist(DDMOutputfile)
                    run = 1;
                else
                    if strcmp(obj.info.runMethod, 'run')
                        run = 1;
                        delete(DDMOutputfile);
                    else
                        if strcmp(obj.info.ddmParam.Scanning, 'on')
                            run = 1;
                        else
                            run = 0;
                        end
                    end
                end

                if run == 1
                    Padsize = zeros(1,3);
                    CriticalAngle = 0;
                    FrameSize = size(obj.AllFrames);
                    NumBins = round(FrameSize(1)./4);
                    nFrames = FrameSize(3);
                    FrameSize(3) = 1;
                    ExpTime = obj.info.ddmParam.ExpTime; 
                    DDMOutput = [];
                    AnisotropyOutput = [];

                    f = waitbar(0, 'initializing');
                    for dt=1:nFrames-1
                        AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3));
                        [AvgFFT] = obj.CalculateDelta(AvgFFT,FrameSize,dt, c); 
                        
                        [RadiallyAveragedDDMSignal, AnisotropyValues]=  obj.AverageRadialy3D(AvgFFT,FrameSize, CriticalAngle, NumBins);
                        DDMOutput(:,1) =[NaN ; RadiallyAveragedDDMSignal(:,1)];
                        DDMOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)];
    
                        if strcmp(obj.info.ddmParam.AngularAnisotropy, 'on')
                            AnisotropyOutput(:,1) = [NaN ; RadiallyAveragedDDMSignal(:,1)];
                            AnisotropyOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)]; 
                        end
                        waitbar(dt./nFrames,f,append('Calculating AvgFFT frame by frame - Timelag ', num2str(dt),...
                            ' out of ', num2str(nFrames-1)));
                    end
                    close all
                    close(f)
                    OldOutput = DDMOutput;
                    DDMOutput =  obj.ConvertOutput(DDMOutput, ExpTime);
                    if strcmp(obj.info.ddmParam.AngularAnisotropy, 'on')
                        AnisotropyOutput = obj.ConvertOutput(AnisotropyOutput, ExpTime);
                    end
                    obj.DDMOutput{c,1} = DDMOutput;

                    if strcmp(obj.info.ddmParam.Scanning, 'off')
                        filename = append(obj.calibrated{1, 1}.mainPath, filesep, 'DDMOutput', num2str(c), '.mat');
                        save(filename, 'DDMOutput')
                        if strcmp(obj.info.ddmParam.AngularAnisotropy, 'on')
                            filename = append(obj.calibrated{1, 1}.mainPath, filesep, 'AnisotropyOutput', num2str(c), '.mat');
                            save(filename, 'AnisotropyOutput')
                        end
                    end
                else
                    disp('Found DDMOuput file - loading it');
                    DDMOutput = load(DDMOutputfile);
                    DDMOutput = DDMOutput.DDMOutput;
                    obj.DDMOutput{c,1} = DDMOutput;
                    disp('Found DDMOuput file - Done');
                end
            end
        end

        function fitDDM(obj)
            for c = 1:size(obj.DDMOutput, 1)
                for i = 1:size(obj.DDMOutput{c,1},1)-1
                    D(i,:) = obj.DDMOutput{c,1}.DDMSignalValue{i,1};
                    q(i) = obj.DDMOutput{c,1}.QVector(i);
                    try
                        A0(i,1) = D(i, find(diff(smooth(D(i,:))) < 0, 1, 'first'));
                    catch
                        A0(i,1) = D(i, end);
                    end
                    B0(i,1) = D(i,1);
                end
                t = obj.DDMOutput{c,1}.Time{1,1} ;
                q1 = obj.info.ddmParam.Qmin;
                if strcmp(obj.info.ddmParam.Scanning, 'on')
                    q2 = 2*pi./(obj.info.PxSize*10^(-3));
                else
                    q2 = obj.info.ddmParam.Qmax;
                end
    
                q_idx = (q >= q1) & (q <= q2);
                q_sel = q(q_idx);
                D_sel = D(q_idx, :);

                params0 = [A0; B0];
                
                options = optimset('MaxIter', 10000, 'Display', 'iter');
                [params_opt, fval, exitflag] = fminsearch(@cost_function, params0, options);
    
                % Extract optimal parameters
                n_q = length(q_sel);
                A_opt = params_opt(1:n_q);
                B_opt = params_opt(n_q+1:end);
                
                obj.FitResults{c,1}.ParamA = A_opt;
                obj.FitResults{c,1}.ParamB = B_opt;
    
                [~, MSD_final] = cost_function(params_opt);
    
                obj.FitResults{c,1}.MSD = MSD_final;
                obj.FitResults{c,1}.tau = [1:size(MSD_final, 2)]*obj.info.ddmParam.ExpTime;
                obj.FitResults{c,1}.MSDstddev = fval;
                obj.FitResults{c,1}.MSDExit = exitflag;
    
                Filename = append(obj.calibrated{1, 1}.mainPath, filesep, 'MSD_', num2str(c), '.mat');
                save(Filename, 'MSD_final');

                ValidPart = find(isnan(obj.FitResults{c,1}.MSD), 1, 'first');
                if ~isempty(ValidPart)
                    obj.FitResults{c,1}.MSD = obj.FitResults{c,1}.MSD(1:ValidPart-1);
                else
                    % obj.FitResults{c,1}.MSD = obj.FitResults{c,1}.MSD(1:(size(obj.FitResults{c,1}.MSD,2)/2));
                end
            end


            function [sigma2, MSD_avg] = cost_function(params)
                % COST_FUNCTION: Computes the dispersion \sigma^2 and MSD for given parameters.
                
                n_q = length(q_sel);
                A = params(1:n_q); % Amplitude A(q)
                B = params(n_q+1:end); % Noise baseline B(q)
                
                % Step 3: Calculate MSD(t|q)
                MSD_tq = zeros(n_q, length(t));
                for i = 1:n_q
                    q_i = q_sel(i);
                    % Ensure no invalid log values by checking (D(q,t) - B(q)) / A(q) < 1
                    valid_idx = (D_sel(i, :) - B(i))./A(i) < 1;
                    MSD_tq(i, valid_idx) = (-4 ./ (q_i^2)) .* log(1 - (D_sel(i, valid_idx) - B(i)) ./ A(i));
                    % For invalid values, set to NaN to exclude from averages
                    MSD_tq(i, ~valid_idx) = NaN;
                end
                
                % Step 4: Determine subset J(t) and calculate N(t)
                J_t = ones(size(MSD_tq));% > (-4 ./ q_sel.^2)'; %< (4 ./ q_sel.^2)';
                N_t = sum(J_t, 1);
                
                % Step 5: Calculate MSD(t)
                MSD_avg = zeros(1, length(t));
                for j = 1:length(t)
                    valid_idx = J_t(:, j);
                    if N_t(j) > 0
                        MSD_avg(1,j) = nanmean(MSD_tq(valid_idx, j));
                    else
                        MSD_avg(j) = NaN; % No valid q values for this time
                    end
                end
                
                % Step 6: Calculate dispersion \sigma^2_t and \sigma^2
                sigma2_t = zeros(1, length(t));
                for j = 1:length(t)
                    valid_idx = J_t(:, j);
                    if N_t(j) > 1
                        log_ratios = log(MSD_tq(valid_idx, j) ./ MSD_avg(j));
                        sigma2_t(j) = sum(log_ratios.^2) / (N_t(j) - 1);
                    else
                        sigma2_t(j) = NaN;
                    end
                end
                
                sigma2 = nansum(sigma2_t); % Total dispersion
            end
        end

        function fitDDMNoOptim(obj)
            for c = 1:size(obj.DDMOutput, 1)
                for i = 1:size(obj.DDMOutput{c,1},1)-1
                    D(i,:) = obj.DDMOutput{c,1}.DDMSignalValue{i,1};
                    q(i) = obj.DDMOutput{c,1}.QVector(i);
                    try
                        A0(i,1) = D(i, find(diff(smooth(D(i,:))) < 0, 1, 'first'));
                    catch
                        A0(i,1) = D(i, end);
                    end
                    B0(i,1) = D(i,1);
                end
                t = obj.DDMOutput{c,1}.Time{1,1} ;
                q1 = obj.info.ddmParam.Qmin;
                if strcmp(obj.info.ddmParam.Scanning, 'on')
                    q2 = 2*pi./(obj.info.PxSize*10^(-3));
                else
                    q2 = obj.info.ddmParam.Qmax;
                end
    
                q_idx = (q >= q1) & (q <= q2);
                q_sel = q(q_idx);
                D_sel = D(q_idx, :);
                A_sel = A0(q_idx);
                B_sel = B0(q_idx);
                q_sel = q(q_idx);

               for i = 1:size(D_sel, 1)
                   ISF = D_sel(i,:);
                   %model = append(num2str(A_sel(i)-B_sel(i)), '*(1-exp(-x*a*(', num2str(q_sel(i)), '^2)))+', num2str(B_sel(i)));
                   model = append('a*(1-exp(-x*b*(', num2str(q_sel(i)), '^2)))+c');
                   f = fit(t', ISF', model);
                   param = coeffvalues(f);
                   A_opt(i,1) = param(1);
                   B_opt(i,1) = param(3);
                   Diffusion(i,1) = param(2);
               end
                
                obj.FitResults{c,1}.ParamA = A_opt;
                obj.FitResults{c,1}.ParamB = B_opt;

                MSD_final = 8*nanmean(rmoutliers(Diffusion))*t;
    
                obj.FitResults{c,1}.MSD = MSD_final;
                obj.FitResults{c,1}.tau = [1:size(MSD_final, 2)]*obj.info.ddmParam.ExpTime;
                obj.FitResults{c,1}.MSDstddev = NaN;
                obj.FitResults{c,1}.MSDExit = NaN;
    
                Filename = append(obj.calibrated{1, 1}.mainPath, filesep, 'MSD_', num2str(c), '.mat');
                Filename(1) = obj.raw.fullPath(1);
                Filename = append(obj.raw.movInfo.Path, filesep, 'calibrated1', filesep, 'MSD_', num2str(c), '.mat');
                save(Filename, 'MSD_final');

                ValidPart = find(isnan(obj.FitResults{c,1}.MSD), 1, 'first');
                if ~isempty(ValidPart)
                    obj.FitResults{c,1}.MSD = obj.FitResults{c,1}.MSD(1:ValidPart-1);
                else
                    % obj.FitResults{c,1}.MSD = obj.FitResults{c,1}.MSD(1:(size(obj.FitResults{c,1}.MSD,2)/2));
                end
            end
        end

        function [AverageDDMValueAtR, anisotropy_values] = AverageRadialy3D(obj, AvgFFT,FrameSize,critangle,NumBins)
            % Averages the scattering function in 3d radially. Returns the
            % average values in function of time and q-vector.

            [RadialValueInQSpace, ValidRange]=obj.Get3DGrid(FrameSize,critangle);
            MinQ = min(RadialValueInQSpace(:));
            MaxQ = max(RadialValueInQSpace(:));
            BinSize = (MaxQ - MinQ)./NumBins;
            FoundRadii = 0;
            R = MinQ;
            next = 1;
            AverageDDMValueAtR = nan(NumBins,2);
            AvgFFT = gather(AvgFFT);

            [N, M] = size(AvgFFT);
            cx = floor(N/2) + 1;
            cy = floor(M/2) + 1;
            num_theta = 360;
            theta_vals = linspace(0, 2*pi, num_theta);

            %%% Calculate Average
            while ~isempty(FoundRadii)
                R = R+BinSize;
                if R <= MaxQ
                    if size(AvgFFT) == size(RadialValueInQSpace)
                        if strcmp(obj.info.ddmParam.AngularAnisotropy, 'on')
                            RadiusQSpace = zeros(size(AvgFFT));
                            RadiusQSpace(RadialValueInQSpace>=R-BinSize & RadialValueInQSpace<R &  ValidRange) = AvgFFT(RadialValueInQSpace>=R-BinSize & RadialValueInQSpace<R &  ValidRange);
                            AngularResults = obj.angularIntensityProfile(RadiusQSpace);
    
                            I_max = max(AngularResults.intensity);
                            I_min = min(AngularResults.intensity);
                            anisotropy_values(next,:) = [R, (I_max - I_min) / (I_max + I_min)];
                        else
                            anisotropy_values = [];
                        end

                        FoundRadii = AvgFFT(RadialValueInQSpace>=R-BinSize & RadialValueInQSpace<R &  ValidRange );
                        AverageDDMValueAtR(next,:) = [R ,  nanmean(FoundRadii, 'all')];
                        next = next+1;

                    elseif size(AvgFFT)
                        error('Differential image and reference frame q vector do not have same dimensions')
                    end
                else
                    break
                end
            end
        end

        function result = angularIntensityProfile(obj, A)
           
            nbins = 360; % default: 1° angular resolution
            [nrows, ncols] = size(A);
            cx = (ncols + 1)/2;
            cy = (nrows + 1)/2;
            [x, y] = meshgrid(1:ncols, 1:nrows);
            vals = A(:);
            x = x(:);
            y = y(:);

            theta = atan2(y - cy, x - cx);  % radians [-pi, pi]
            theta = mod(theta, 2*pi);       % convert to [0, 2π)

            nonzero = vals ~= 0;
            theta = theta(nonzero);
            vals = vals(nonzero);

            edges = linspace(0, 2*pi, nbins+1);
            [~, ~, bin] = histcounts(theta, edges);

            angular_mean = accumarray(bin, vals, [nbins, 1], @mean, NaN);

            theta_centers = (edges(1:end-1) + edges(2:end)) / 2;

            valid = ~isnan(angular_mean);
            theta_centers = theta_centers(valid);
            angular_mean = angular_mean(valid);

            result.theta = theta_centers;
            result.angle_deg = rad2deg(theta_centers);
            result.intensity = angular_mean;         
        end


         function  [RadialValueInQSpace,ValidIndeces]=Get3DGrid(obj,sizes,critangle)

            z = 0;
            for i = 1:2
                if rem(sizes(i), 2) == 0
                    x{i} = 2*pi*(-round(sizes(i)/2):1:round(sizes(i)/2)-1)./((obj.info.PxSize*10^(-3)).*sizes(i));
                elseif rem(size(i), 2) == 1
                    x{i} = 2*pi*(-round(sizes(i)/2)+1:1:round(sizes(i)/2)-1)./((obj.info.PxSize*10^(-3)).*sizes(i));
                end
            end

            [qx, qy, qz] = meshgrid(x{2}, x{1}, z);
            RadialValueInQSpace = sqrt(qx.^2 + qy.^2 + qz.^2);
            ValidIndeces =  qz./sqrt(qy.^2 + qx.^2 + qz.^2)<abs(cosd(critangle));
        end

         function [AvgFFT] =  CalculateDelta(obj,AvgFFT,FrameSize,dt,c)
            % DiffFrames = {};
            counts = 0;
            nFrames = size(obj.AllFrames, 3);
            for t=1:nFrames-dt
                FrameDelta = obj.AllFrames(:,:,t+dt,c)-obj.AllFrames(:,:,t,c);  
                % DiffFrames{1, end+1} = fftshift(fft2(FrameDelta-mean(FrameDelta,'all')));
                % DiffFrames{2, end} = FrameDelta-mean(FrameDelta,'all');
                FrameDelta = abs( fftshift(fftn(FrameDelta-mean(FrameDelta,'all'),[FrameSize(1),FrameSize(2),FrameSize(3)]))).^2;                  
                AvgFFT = AvgFFT + FrameDelta;  
                counts = counts + 1;
            end
            AvgFFT = AvgFFT./counts;
         end 

         function DDMOutput = ConvertOutput(obj, DDMOutList,ExpTime)
            DDMOutput = table(0,{[]},{[]},'VariableNames',{'QVector','Time','DDMSignalValue'});
            warning('off','all')
            for i=2:size(DDMOutList,1)
                DDMOutput.QVector(i-1) = DDMOutList(i,1);
                DDMOutput.Time(i-1) = {DDMOutList(1,2:end)*ExpTime};
                DDMOutput.DDMSignalValue(i-1) ={DDMOutList(i,2:end)};
    
            end
            warning('on','all')
         end

         function getParams(obj)
             for c = 1:size(obj.FitResults,1)
                 Results{1,c}.Diff = MSD.getDiffCoeff(obj.FitResults{1, 1}.MSD,obj.FitResults{1, 1}.tau,...
                     obj.info.ddmParam.FitRDiff,obj.info.Dimension);
                
                 switch obj.info.Dimension
                    case '1D'
                        equation = '2*a*x+b';
                    case '2D'
                        equation = '4*a*x+b';
                    case '3D'
                        equation = '6*a*x+b';
                    otherwise
                        error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
                 end
            
                 assert(min(size(obj.FitResults{1, 1}.MSD))==1,'MSD needs to be provided as a vector')
                 %  assert(and(fitRange<=1,isnumeric(fitRange)),'fit Range needs to be numerical between 0 and 1');

                 tofit = obj.FitResults{1, 1}.MSD(1:obj.info.ddmParam.FitRDiff);
                 tau   = obj.FitResults{1, 1}.tau(1:obj.info.ddmParam.FitRDiff);
                 Lower = [obj.info.FitMinDiffLim, min(tofit)];
                 Upper = [obj.info.FitMaxDiffLim, max(tofit)];
                 Start = [obj.info.FitDiffEstim, tofit(1) - (tofit(2) - tofit(1))];
                 % [f, gov]     = fit(tau(:),tofit(:),equation, 'Lower', Lower, 'Upper', Upper, 'StartPoint', Start);
                 [f, gov]     = fit(tau(:),tofit(:),'4*a*x');

                 % if gov.rsquare > 0
                     g = coeffvalues(f);
                     D = g(1);
                     Results{1,c}.Diff = D;
                     try
                        Results{1,c}.alpha = MSD.getDiffTypeAlpha2(obj.FitResults{1, 1}.MSD',obj.info.ddmParam.ExpTime,sqrt(Results{1,c}.Diff)*obj.info.ddmParam.ExpTime);
                     catch
                         Results{1,c}.alpha = NaN;
                     end
                     Results{1,c}.n = MSD.getViscosity(Results{1,c}.Diff,obj.info.ddmParam.ParticleSize*10^(-3),obj.info.ddmParam.Temp);
                 % else
                 %     Results{1,c}.Diff = NaN;
                 %     Results{1,c}.alpha = NaN;
                 %     Results{1,c}.n = NaN;
                 % end
             end

             obj.MSDResults = Results;
             Filename = append(obj.calibrated{1, 1}.mainPath, filesep, 'MSDResults.mat');
             Filename(1) = obj.raw.fullPath(1);
             Filename = append(obj.raw.movInfo.Path, filesep, 'calibrated1', filesep, 'MSDResults.mat');
             save(Filename, "Results");
         end
    end
end

