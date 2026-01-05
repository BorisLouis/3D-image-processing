classdef TranslationalTracking < handle

    properties
        raw
        info
        Traces
        Results
        ResultsStepsize
    end
    
    methods
        function obj = TranslationalTracking(Path, info)
            obj.raw.Path = Path;
            obj.info = info;
        end

        function LoadTraces(obj, Filename)
            f2Load = append(obj.raw.Path, filesep, Filename, '.mat');
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});

            if ~isnan(obj.info.CutTraces)
                if isfield(data, 'traces')
                    data = data.traces;
                end

                if size(data, 2) < size(data, 1)
                    data = data';
                end
                mSizes = cellfun(@(t) size(t,1), data(1,:));
                mask = mSizes > obj.info.MinSize;
                data = data(:,mask)';

                Length = obj.info.CutTraces;
                f = waitbar(0, 'initializing');
                for step = 1:max([data{end,1}.t])-Length
                    waitbar(step./(max([data{end,1}.t])-Length), f, append('Cutting traces - Frame ', num2str(step), '/', num2str(max([data{end,1}.t])-Length)))
                    CurrTraceCell = {};
                    for i = 1:size(data,1)
                        CurrTrace = data{i,1};
                        CurrTraceCutted = CurrTrace(ismember(CurrTrace.t, (step:(step+Length))), :);    
                        if ~isempty(CurrTraceCutted)
                            if size(CurrTraceCutted, 1) > obj.info.MinSize
                                CurrTraceCell{end+1,1} = CurrTraceCutted;
                            end
                        end
                    end
                    dataMatrix{step, 1} = CurrTraceCell;
                end
                close(f)
            else
                try
                    Data = data.traces;
                catch
                    Data = data;
                end
                mSizes = cellfun(@(t) size(t,1), Data(1,:));
                mask = mSizes > obj.info.MinSize;
                Data = Data(:,mask)';
                dataMatrix{1, 1} = Data;
            end
            obj.Traces = dataMatrix;
        end

        function Analysis(obj, Radius, Loop)
            % if Loop == 1
            %     if exist(append(obj.raw.Path, filesep, 'msdRes', num2str(1), '.mat'))
            %         Loaded = load(append(obj.raw.Path, filesep, 'msdRes', num2str(1), '.mat'));
            %         Results{1, 1} = Loaded.Results;
            %         run = 0;
            %     else
            %         run = 1;
            %     end
            % elseif Loop == 2
            %     if exist(append(obj.raw.Path, filesep, 'msdRes', num2str(2), '.mat'))
            %         Loaded2 = load(append(obj.raw.Path, filesep, 'msdRes', num2str(2), '.mat'));
            %         Results{2, 1} = Loaded2.Results;
            %         run = 0;
            %     else
            %         run = 1;
            %     end
            % end
                
            run = 1;
            if run == 0 
                obj.Results = Results;
            elseif run == 1
                nRows = size(obj.Traces, 1);

                try
                    MinIdx = strfind(obj.raw.Path, 'min');
                    TimeStamp = obj.raw.Path(MinIdx+3:MinIdx+4);
                    TimeStamp = str2num(erase(TimeStamp, '_'));
                catch
                    TimeStamp = 0;
                end
    
                Time     = nan(nRows, 1);
                DiffMean = nan(nRows, 1);
                DiffStd  = nan(nRows, 1);
                ViscMean = nan(nRows, 1);
                ViscStd  = nan(nRows, 1);
                AnExpMean = nan(nRows, 1);
                AnExpStd  = nan(nRows, 1);
                DiffAll = cell(nRows, 1);
                ViscAll = cell(nRows, 1);
                AnExpAll = cell(nRows, 1);
                TimeResults = table(Time, DiffMean,DiffStd,DiffAll,ViscMean, ViscStd, ViscAll,AnExpMean,AnExpStd,AnExpAll);
    
                f = waitbar(0, 'initializing');
                %for k = 1:nRows
                for k = 1:20
                    currMov = obj.Traces{k, 1};
                    AllStepSizes = [];
                    if ~isempty(currMov)
                        allRes = struct('msdx',0,'msdy',0,'msdz',0,'msdr',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
                            'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
                            'vZ',0,'vR',0,'Gcomplex_x', 0, 'Gcomplex_y', 0, 'Gcomplex_z', 0, 'Gcomplex_r', 0,...
                            'Gloss_x', 0, 'Gloss_y', 0, 'Gloss_z', 0, 'Gloss_r', 0,...
                            'Gstorage_x', 0, 'Gstorage_y', 0, 'Gstorage_z', 0, 'Gstorage_r', 0, 'RoG_x', 0,...
                            'RoG_y', 0, 'RoG_z', 0', 'RoG_r', 0,...
                            'RoGTensor_x', 0, 'RoGTensor_y', 0, 'RoGTensor_z', 0, 'RoGTensor_r', 0,...
                            'Lp_x', 0, 'Lp_y', 0, 'Lp_z', 0, 'Lp_r', 0,...
                            'EndToEnd_x', 0, 'EndToEnd_y', 0, 'EndToEnd_z', 0, 'EndToEnd_r', 0,...
                            'num', 0);
                        if strcmp(obj.info.Experiment, 'Tracking-Segmentation')
                            allRes.Mask = 0;
                        elseif strcmp(obj.info.Experiment, 'Tracking-Phase')
                            allRes.Phase = 0;
                        end
                        allRes(length(currMov)).msdX = [];
                        maxLength = max(cellfun(@height, currMov(:,1)));
    
                        for i = 1:length(currMov)
                            waitbar(i./length(currMov), f, append('Microrheology analysis - part ', num2str(k), '/', num2str(size(obj.Traces, 1))));
                            currPart = currMov{i};
                        
                            coordinates = [currPart.col, currPart.row, currPart.z];
                            CM = mean(coordinates,1);
                            coordinates = coordinates-CM;
                        
                            %in X
                            coord = coordinates(:,1)/10^3;
                            Dimension = '1D';
                            % [~, allRes(i).msdx, allRes(i).tau, allRes(i).DX, allRes(i).nX, allRes(i).aX, allRes(i).vX] = obj.TraceAnalysis(coord, Dimension, Radius);
                            [allRes(i).StepSizes_x] = obj.GetStepSizes(coord);
                            % [allRes(i).Gcomplex_x, allRes(i).Gstorage_x, allRes(i).Gloss_x] = obj.PassiveMicrorheology(allRes(i).msdx, allRes(i).tau, Radius, 0.5);
                            % [allRes(i).RoGTensor_x, allRes(i).RoG_x] = obj.GyrationTensor(coord);
                            % [allRes(i).Lp_x, allRes(i).EndToEnd_x] = obj.PersistenceLength(coord);
                            % 
                            %in Y
                            coord = coordinates(:,2)/10^3;
                            Dimension = '1D';
                            % [~, allRes(i).msdy, ~, allRes(i).DY, allRes(i).nY, allRes(i).aY, allRes(i).vY] = obj.TraceAnalysis(coord, Dimension, Radius);
                            [allRes(i).StepSizes_y] = obj.GetStepSizes(coord);
                            % [allRes(i).Gcomplex_y, allRes(i).Gstorage_y, allRes(i).Gloss_y] = obj.PassiveMicrorheology(allRes(i).msdy, allRes(i).tau, Radius, 0.5);
                            % [allRes(i).RoGTensor_y, allRes(i).RoG_y] = obj.GyrationTensor(coord);
                            % [allRes(i).Lp_y, allRes(i).EndToEnd_y] = obj.PersistenceLength(coord);

                            %inZ
                            if strcmp(obj.info.Dimension, '3D')
                                coord = coordinates(:,3)/10^3;
                                Dimension = '1D';
                                % [~, allRes(i).msdz, ~, allRes(i).DZ, allRes(i).nZ, allRes(i).aZ, allRes(i).vZ] = obj.TraceAnalysis(coord, Dimension, Radius);
                                [allRes(i).StepSizes_z] = obj.GetStepSizes(coord);
                                % [allRes(i).Gcomplex_z, allRes(i).Gstorage_z, allRes(i).Gloss_z] = obj.PassiveMicrorheology(allRes(i).msdz, allRes(i).tau, Radius, 0.5);
                                % [allRes(i).RoGTensor_z, allRes(i).RoG_z] = obj.GyrationTensor(coord);
                                % [allRes(i).Lp_z, allRes(i).EndToEnd_z] = obj.PersistenceLength(coord);
                            else
                                allRes(i).msdz = [];
                                allRes(i).DZ = NaN;
                                allRes(i).nZ = NaN;
                                allRes(i).aZ = NaN;
                                allRes(i).vZ = NaN;
                                allRes(i).Gcomplex_z = [];
                                allRes(i).Gloss_z = [];
                                allRes(i).Gstorage_z = [];
                                allRes(i).RoGTensor = [];
                                allRes(i).RoG = NaN;
                                allRes(i).Lp_z = NaN;
                                allRes(i).EndToEnd_z = NaN;
                            end
    
                            %inR
                            if strcmp(obj.info.Dimension, '3D')
                                coord = coordinates(:,1:3)/10^3;
                                Dimension = '3D';
                            elseif strcmp(obj.info.Dimension, '2D')
                                coord = coordinates(:,1:2)/10^3;
                                Dimension = '2D';
                            end          
                            % [~, allRes(i).msdr, ~, allRes(i).DR, allRes(i).nR, allRes(i).aR, allRes(i).vR] = obj.TraceAnalysis(coord, Dimension, Radius);
                            [allRes(i).StepSizes_r] = obj.GetStepSizes(coord);
                            % [allRes(i).Gcomplex_r, allRes(i).Gstorage_r, allRes(i).Gloss_r] = obj.PassiveMicrorheology(allRes(i).msdr, allRes(i).tau, Radius, 0.5);
                            % [allRes(i).RoGTensor_r, allRes(i).RoG_r] = obj.GyrationTensor(coord);
                            % [allRes(i).Lp_r, allRes(i).EndToEnd_r] = obj.PersistenceLength(coord);
    
                            allRes(i).num  = length(allRes(i).msdr);
                    
                            % if strcmp(obj.info.Experiment, 'Tracking-Segmentation')
                            %     allRes(i).Mask = round(mean(currPart.InSegment));
                            % elseif strcmp(obj.info.Experiment, 'Tracking-Phase')
                            %     try
                            %         allRes(i).Phase = nanmean(currPart.Phase);
                            %         allRes(i).IntPhaseCh = nanmean(currPart.IntPhaseCh);
                            %         allRes(i).GradientMagnitude = nanmean(currPart.GradientMagnitude);
                            %         allRes(i).LocalVariance = nanmean(currPart.LocalVariance);
                            %         allRes(i).SharpnessLaplacian = nanmean(currPart.SharpnessLaplacian);
                            %     catch
                            %         allRes(i).Phase = NaN;
                            %         allRes(i).IntPhaseCh = NaN;
                            %         allRes(i).GradientMagnitude = NaN;
                            %         allRes(i).LocalVariance = NaN;
                            %         allRes(i).SharpnessLaplacian = NaN;
                            %     end
                            % end
                            if all(size(AllStepSizes) == [0 0])
                                AllStepSizes = [AllStepSizes; allRes(i).StepSizes_r];
                            elseif size(AllStepSizes, 2) > size(allRes(i).StepSizes_r, 2)
                                ToAdd = [allRes(i).StepSizes_r, nan(size(allRes(i).StepSizes_r,1), size(AllStepSizes, 2) - size(allRes(i).StepSizes_r, 2))];
                                AllStepSizes = [AllStepSizes; ToAdd];
                            elseif size(AllStepSizes, 2) < size(allRes(i).StepSizes_r, 2)
                                ToAdd = [AllStepSizes, nan(size(AllStepSizes,1), size(allRes(i).StepSizes_r, 2) - size(AllStepSizes, 2))];
                                AllStepSizes = [ToAdd; allRes(i).StepSizes_r];
                            else
                                AllStepSizes = [AllStepSizes; allRes(i).StepSizes_r];
                            end
                        end
    
                        disp(['The diffusion coefficient is ', num2str(mean([allRes.DR], 'omitnan')), ' \mum^2/s and the viscosity is ', num2str(mean([allRes.nR], 'omitnan')), ' cp']);
    
                        TimeResults.Time(k) = TimeStamp + k*obj.info.expTime;
                        TimeResults.DiffMean(k) = mean([allRes.DR]); 
                        TimeResults.DiffStd(k)  = std([allRes.DR]);
                        TimeResults.ViscMean(k) = mean([allRes.nR]);
                        TimeResults.ViscStd(k)  = std([allRes.nR]);
                        TimeResults.AnExpMean(k) = mean([allRes.aR]);
                        TimeResults.AnExpStd(k)  = std([allRes.aR]);
                        TimeResults.DiffAll{k} = [allRes.DR];
                        TimeResults.ViscAll{k} = [allRes.nR];
                        TimeResults.AnExpAll{k} = [allRes.aR];  
                        TimeResults.AllStepSizes{k} = AllStepSizes;
                    else
                        %disp(append('No traces found that are longer than MinSize (', num2str(MinSize), ' datapoints)'))
                        TimeResults.Time(k) = k;
                    end
    
                    Results{k,1} = allRes;
                    %Results{k,2} = TimeResults;
                end
                obj.Results{Loop, 1} = Results;
                filename = append(obj.raw.Path, filesep, 'msdRes', num2str(Loop), '.mat');
                save(filename, "Results");
                Results{end,2} = TimeResults;
                obj.Results{Loop, 1} = Results;
                filename = append(obj.raw.Path, filesep, 'msd_TimeResults', num2str(Loop), '.mat');
                save(filename, "TimeResults");
            end
        end

        function PlotTrends(obj)
            fig = figure;
            baseColors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
            for i = 1:size(obj.Results, 1)
                try  
                    Data = obj.Results{i,1}{end,2};
                    Data(isnan([Data.Time]),:) = [];
                    Time = Data.Time;
                    DiffMean = Data.DiffMean;
                    DiffStd = Data.DiffStd; 

                    baseColor = baseColors(i, :);  % MATLAB default blue
                    lighterColor = baseColor + 0.5 * (1 - baseColor); % lighter for shading

                    upper = DiffMean + DiffStd;
                    lower = DiffMean - DiffStd;

                    hold on;
                    fill([Time; flipud(Time)], ...
                         [upper; flipud(lower)], ...
                         lighterColor, ...
                         'FaceAlpha', 0.8, 'EdgeColor', 'none');
                    hold on
                    plot(Time, DiffMean, 'Color', baseColor, 'LineWidth', 2);
                catch
                end
            end
            xlabel('Time (s)', 'FontSize', 12);
            ylabel('Diffusion (µm^2/s)', 'FontSize', 12);
            grid on; box on; axis tight;
            set(gca, 'FontSize', 10);
            legend({'', append(num2str(obj.info.Radius1*100), ' nm'), '', append(num2str(obj.info.Radius2*100), ' nm')}, 'Location', 'best');
            title('Diffusion over Time');
            saveas(fig, append(obj.raw.Path, filesep, 'DiffusionTrend.png'));


            fig = figure;
            baseColors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
            for i = 1:size(obj.Results, 1)
                try  
                    Data = obj.Results{i,1}{end,2};
                    Data(isnan([Data.Time]),:) = [];
                    Time = Data.Time;
                    ViscMean = Data.ViscMean;
                    ViscStd = Data.ViscStd; 

                    baseColor = baseColors(i, :);  % MATLAB default blue
                    lighterColor = baseColor + 0.5 * (1 - baseColor); % lighter for shading

                    upper = ViscMean + ViscStd;
                    lower = ViscMean - ViscStd;

                    hold on;
                    fill([Time; flipud(Time)], ...
                         [upper; flipud(lower)], ...
                         lighterColor, ...
                         'FaceAlpha', 0.8, 'EdgeColor', 'none');
                    hold on
                    plot(Time, ViscMean, 'Color', baseColor, 'LineWidth', 2);
                catch
                end
            end
            xlabel('Time (s)', 'FontSize', 12);
            ylabel('Viscosity (cP)', 'FontSize', 12);
            grid on; box on; axis tight;
            set(gca, 'FontSize', 10);
            set(gca, 'YScale', 'log');
            legend({'', append(num2str(obj.info.Radius1*1000), ' nm'), '', append(num2str(obj.info.Radius2*1000), ' nm')}, 'Location', 'best');
            title('Viscosity over Time');
            ylim([0 50]);
            saveas(fig, append(obj.raw.Path, filesep, 'ViscosityTrend.png'));
        end

        function [AvStep, msdx, tau, D, n, a, v] = TraceAnalysis(obj, coord, Dimension, Radius)
            AvStep = MSD.getAvStepSize(coord); 
            msdx = MSD.calc(coord);%convert to um;
            tau = (1:length(msdx))'*obj.info.expTime;
            D  = MSD.getDiffCoeff(msdx,tau,obj.info.DiffFit,Dimension);
            n  = MSD.getViscosity(D,Radius,obj.info.Temp);
            a  = MSD.getDiffTypeAlpha2(msdx,obj.info.expTime, AvStep);
            if strcmp(Dimension, '1D')
                v  = abs(coord(1,1) - coord(end,1)/10^3/(length(coord)*obj.info.expTime)); %um/s
            elseif strcmp(Dimension, '3D')
                d   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2 + (coord(1,3)-coord(end,3))^2);
                v = d/(length(coord)*obj.info.expTime); %um/s
            elseif strcmp(Dimension, '2D')
                d   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2);
                v = d/(length(coord)*obj.info.expTime); %um/s
            end
        end

        function PlotDistributions(obj, Loop)
            m = size(obj.Results{Loop,1},1);
            DR_all = cell(m,1);
            allDR = []; 
            for i = 1:m
                if strcmp(obj.info.StepsizeAnalysis, 'on')
                    StepMatrix = (obj.Results{Loop,1}{end,2}.AllStepSizes{i,1}(:,1)).^2./obj.info.expTime;
                    DR_all{i} = StepMatrix(:);
                    allDR = [allDR; StepMatrix];
                else
                    DR = [obj.Results{Loop,1}{i,1}.DR];
                    DR_all{i} = DR(:);      % ensure column vector
                    allDR = [allDR, DR];
                end
            end
            lowEdge  = 0;     % lower bound: 1st percentile
            highEdge = 10;    % upper bound: 99th percentile

            numBins = 50;                     % adjust as needed
            edges = linspace(lowEdge, highEdge, numBins+1);
            binCenters = edges(1:end-1) + diff(edges)/2;
            H = zeros(m, numBins);

            for t = 1:m
                H(t,:) = histcounts(DR_all{t}, edges);
            end

            Fig = figure;
            surf(binCenters, [1:m].*obj.info.expTime, H, 'EdgeColor', 'none');
            view(2); % optional: use view(3) for perspective
            
            xlabel('Diffusion coefficient (µm^2/s)');
            ylabel('Time (s)');
            zlabel('Counts (abs)');
            a=colorbar;
            a.Label.String = 'Counts (abs)';
            title(append('Diffusion coefficient over time - ', num2str(1000*obj.info.(append('Radius', num2str(Loop)))), ' nm'));

            if strcmp(obj.info.StepsizeAnalysis, 'on')
                title(append('Diffusion coefficient over time (1st order stepsize) - ', num2str(1000*obj.info.(append('Radius', num2str(Loop)))), ' nm'));
                saveas(Fig, append(obj.raw.Path, filesep, 'DistributionPlot_Diffusion_Stepsize', num2str(Loop), '.png'));
            else
                saveas(Fig, append(obj.raw.Path, filesep, 'DistributionPlot_Diffusion', num2str(Loop), '.png'));
            end
        end

        function [Gmag, Gp, Gpp] = PassiveMicrorheology(obj, msd, tau, Radius, fracKeep)
            nPts = numel(msd);
            nKeep = max(5, round(fracKeep*nPts));
            msd = msd(1:nKeep);
            tau = tau(1:nKeep);

            logTau = log10(tau);
            logMSD = log10(msd);
            windowLength = min(11, numel(tau)); % must be odd, <= number of points
            if mod(windowLength,2)==0
                windowLength = windowLength-1; % make odd
            end
            polyOrder = 2;
            smoothLogMSD = sgolayfilt(logMSD, polyOrder, windowLength);
            alpha = gradient(smoothLogMSD) ./ gradient(logTau);
            alpha(alpha < 0) = 0;
            alpha(alpha > 1) = 1;
            
            % complex modulus magnitude (Mason–Weitz GSER)
            Gmag = (1.3806498*10^(-23)*obj.info.Temp) ./ (pi*Radius*10^(-9)*msd*10^-(12) .* gamma(1+alpha));
            
            % storage and loss moduli
            Gp  = Gmag .* cos(pi*alpha/2);
            Gpp = Gmag .* sin(pi*alpha/2);
        end

        function [G, Rg] = GyrationTensor(obj, coords)
            [N, d] = size(coords);
            if d < 1 || d > 3
                error('coords must be an n×1, n×2, or n×3 matrix.');
            end

            r_mean = mean(coords, 1);
            R = coords - r_mean;
            G = (R' * R) / N;
            Rg = sqrt(trace(G));
        end

        function [Lp, Re] = PersistenceLength(obj, coords)
            % PersistenceLength estimates the persistence length of a polymer-like
            % trace using the tangent–tangent correlation method -
            % Worm-like chain model
     
            [N, d] = size(coords);
            if d < 1 || d > 3
                error('coords must be n×1, n×2, or n×3');
            end

            % --- Compute tangent vectors ---
            dr = diff(coords, 1, 1);            % segment vectors
            segLengths = sqrt(sum(dr.^2, 2));   % segment lengths
            t = dr ./ segLengths;               % normalized tangents (N-1 x d)
        
            % arc-length step (mean step length)
            ds = mean(segLengths);
        
            % --- Tangent-tangent correlation function ---
            maxLag = floor((N - 1) / 2);   % enough statistics
        
            C = zeros(maxLag,1);
            s = (1:maxLag)' * ds;          % separation distances
        
            for k = 1:maxLag
                dots = sum(t(1:end-k,:) .* t(1+k:end, :), 2);
                C(k) = mean(dots);
            end
        
            % --- Fit C(s) = exp(-s / Lp) using linear regression on ln(C) ---
            % Only use positive values of C for log
            valid = C > 0.05;  % threshold removes noise
        
            if sum(valid) < 3
                Lp = NaN;
            else
                y = log(C(valid));
                x = s(valid);
            
                % linear fit: ln(C) = -s/Lp
                xm = mean(x);
                ym = mean(y);
                slope = sum((x - xm) .* (y - ym)) / sum((x - xm).^2);
                Lp = -1 / slope;
            
                Lp = -1 / slope;    % persistence length
            end
        
            % --- End-to-end distance ---
            Re = norm(coords(end,:) - coords(1,:));
        
        end

        function [StepSizes] = GetStepSizes(obj, coords)
            [N, d] = size(coords);
            stepVectors = diff(coords, 1, 1);   % (N-1 × d)
            StepSizes = NaN(N, N-1);
            for tau = 1:N-1
                StepSizeList = [];
                for start = 1:tau
                    IdxList = [1:N]+(start-1);
                    idx = mod(IdxList, tau) == 0;
                    SelectedCoords = coords(idx, :);
                    stepVectors = diff(SelectedCoords,1,1);
                    StepSizeList = [StepSizeList; sqrt(sum(stepVectors.^2, 2))];
                end
                StepSizes(1:size(StepSizeList, 1), tau) = StepSizeList;
            end
        end

        function StepSizeAnalysis(obj,Loop)
            cdf1 = @(r2, m1) 1 - exp(-r2./m1);
            cdf2 = @(r2, a, m1, m2) 1 - ( a .* exp(-r2./m1) + (1-a).*exp(-r2./m2) );
            cdf3 = @(r2, a1, a2, m1, m2, m3) 1 - ( a1.*exp(-r2./m1) + a2.*exp(-r2./m2) + (1-a1-a2).*exp(-r2./m3));
            TableOfDiffusionCoeff = [];
            TableOfFractions = [];

            for Frame = 1:size(obj.Results{Loop,1},1)
                steps = obj.Results{Loop,1}{end,2}.AllStepSizes{Frame,1};

                numLag = size(steps,2);
                fitResults = struct;

                meanR = [0, 0, 0];
                
                for k = 1
                    r = steps(:,k);
                    r2 = r.^2;
                    r2 = r2(~isnan(r2) & r2>0);
                    [F_emp, r2_sorted] = ecdf(r2);
                
                    % Initial guesses
                    m1_0 = median(r2);
                    m2_0 = 2*m1_0;
                    m3_0 = 3*m1_0;
                
                    % Fit one-component
                    try
                        model1 = @(p,x) cdf1(x,p(1));
                        p1{k} = nlinfit(r2_sorted, F_emp, model1, m1_0);
                    catch
                        p1{k} = [nan];
                    end
                
                    % % Fit two-component
                    try
                        model2 = @(p,x) cdf2(x,p(1),p(2),p(3));
                        p2{k} = nlinfit(r2_sorted, F_emp, model2, [0.5, m1_0, m2_0]);
                    catch
                        p2{k} = [nan nan];
                    end
        
                    % % Fit three-component
                    try
                        model3 = @(p,x) cdf3(x,p(1),p(2),p(3), p(4), p(5));
                        p3{k} = nlinfit(r2_sorted, F_emp, model3, [0.33, 0.33, m1_0, m2_0, m3_0]);
                    catch
                        p3{k} = [nan nan nan];
                    end
    
                    r2 = (steps(:,k)).^2;
                    if and(all(p3{k}(1:2) > 0.05),  all(p3{k}(1:2) < 1))
                        bestModel{Frame, k} = 'three-component';
                        m(1) = p3{k}(3);
                        m(2) = p3{k}(4);
                        m(3) = p3{k}(5);
                        a(1) = p3{k}(1);
                        a(2) = p3{k}(2);
                        a(3) = 1-a(1)-a(2);
                        [~, minIdx] = min([abs(r2 - m(1)*k), abs(r2 - m(2)*k), abs(r2 - m(3)*k)], [], 2);
                        pop{k,1} = r2(minIdx == 1);
                        pop{k,2} = r2(minIdx == 2);
                        pop{k,3} = r2(minIdx == 3); 
                        maxPops = 3;
                    elseif abs(min(p2{k}(2:3))./max(p2{k}(2:3))) < 0.90
                        bestModel{Frame, k} = 'two-component';
                        m(1) = p2{k}(2);
                        m(2) = p2{k}(3);
                        a(1) = p2{k}(1);
                        a(2) = 1 - a(1);
                        [~, minIdx] = min([abs(r2 - m(1)*k), abs(r2 - m(2)*k)],[],2);
                        pop{k,1} = r2(minIdx == 1);
                        pop{k,2} = r2(minIdx == 2);
                        maxPops = 2;
                    else
                        try
                            bestModel{Frame, k} = 'one-component';
                            m(1) = p1{k};
                            a(1) = 1;
                            pop{k,1} = r2;
                            maxPops = 1;
                        catch
                            m = NaN;
                            a = NaN;
                            pop{k,1} = r2;
                            maxPops = 1;
                        end
                    end
                end  

                if strcmp(obj.info.Dimension, '1D')
                    n = 1;
                elseif strcmp(obj.info.Dimension, '2D')
                    n = 2;
                elseif strcmp(obj.info.Dimension, '3D')
                    n = 3;
                end
                m = [m, nan(1, 3 - size(m, 2))];
                % TableOfDiffusionCoeff = [TableOfDiffusionCoeff; m];
                TableOfDiffusionCoeff = [TableOfDiffusionCoeff; (m)./(2*n*obj.info.expTime)];
                a = [a, nan(1, 3 - size(a, 2))];
                TableOfFractions = [TableOfFractions; a];
                m = [];
                a = [];
                Time(Frame, 1) = Frame.*obj.info.expTime;
            end

            Results.Diff = TableOfDiffusionCoeff;
            Results.Fraction = TableOfFractions;
            obj.ResultsStepsize = Results;

            save(append(obj.raw.Path, filesep, 'ResultsStepsizeAnalysis_', num2str(Loop), '.mat'), "Results");

            Fig = figure; hold on
            colors = lines(3);  % 3 distinct color         
            for k = 1:3
                scatter(Time, TableOfDiffusionCoeff(:,k), 30, colors(k,:),'filled', ...
                    'MarkerFaceAlpha','flat', 'AlphaData', TableOfFractions(:,k) );
                LegendNames{1,k} = append('Population ', num2str(k));
            end
            
            xlabel('Polymerization time (s)');
            ylabel('Diffusion coefficient (µm^2/s)');
            ylim([0 5])
            legend(LegendNames);
            grid on
            saveas(Fig, append(obj.raw.Path, filesep, 'DiffusionTrend_', num2str(Loop), '.png'));
        end
    end
end

