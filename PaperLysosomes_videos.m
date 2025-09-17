close all; clc; clear all

HeadFolder = 'D:\MultiColor - lysosome tracking\raw data for presentation';
HeadFolders = dir(HeadFolder);
MaxFrame = [500, 500, 500, 200, 300, 300, 300, 300, 300, 300];
ExpTimes = [30, 30, 30, 200, 100, 100, 100, 100, 100, 100];
CellNames = {'A549 - mSiPEI', 'A549 - mSi', 'HeLa - mSiPEI', 'HeLa - mSi', 'HepG2 - mSiPEI', 'HepG2 - mSi',...
    'KM12C - mSiPEI', 'KM12C - mSi', 'MCF - mSiPEI', 'MCF - mSi'};
imgH = 512;
imgW = 512;

for s = 4:size(HeadFolders,1)
    try
        mainFolder1 = append(HeadFolders(s).folder, filesep, HeadFolders(s).name);
    
        %% Parameters
        CellName = CellNames{s-2};
        framesToLoad = 1:MaxFrame(s-2);
        frameRate = round(1000./ExpTimes(s-2));     % frames per second
        pxSize = 81;       % nm per pixel
        TimeRange = {'0','30'};
        stepThreshold = 500;
        DMin = 0;
        DMax = 0.25;
        Amin = 0;
        Amax = 2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Background Particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
        try
            %% Paths for dataset 1 (left)
            trackFile1  = fullfile(mainFolder1,'TracesWMask.mat');
        
            %% Video output names
            mp4File_seg = fullfile(mainFolder1,'TraceMovie_Segment_bgParticles.mp4');
            mp4File_diff = fullfile(mainFolder1,'TraceMovie_Diffusion_bgParticles.mp4');
            mp4File_alpha = fullfile(mainFolder1,'TraceMovie_Alpha_bgParticles.mp4');
            
            %% Load track data
            load(trackFile1)%,'traces3D'); 
            traces1 = Traces3D;
            
            %% Load raw videos
            Folder = dir(mainFolder1);
            idx = find(and(contains({Folder.name}, '.HIS'), contains({Folder.name}, 'Cam2')) == 1);
            videoFile1 = fullfile(mainFolder1, Folder(idx).name);
            mov1 = Load.Movie.his.getFrame(videoFile1, framesToLoad);
            [imgH, imgW, nFrames] = size(mov1);
            nY = imgH; nX = imgW;
            
            % Parameters for scalebar
            scaleBarLength_nm = 10000; % 10 µm
            scaleBarLength_px = round(scaleBarLength_nm / pxSize);
            barHeight = 10;   % pixels thick
            margin = 10;     % pixels from bottom-left
            
            % add scalebar to each frame (bright line)
            for i = 1:nFrames
                frameImg = mov1(:,:,i);  % uint16 frame
                brightValue = max(frameImg(:)); % maximum for uint16 image
                % Bottom-right coordinates (y increases downward; top-left is (1,1))
                rowStart = nY - margin - barHeight;
                rowEnd   = rowStart + barHeight - 1;
                colEnd   = nX - margin;
                colStart = colEnd - scaleBarLength_px + 1;
                if rowStart < 1, rowStart = 1; end
                if colStart < 1, colStart = 1; end
            
                frameImg(rowStart:rowEnd, colStart:colEnd) = brightValue;
            
                % Insert '10 µm' text above bar — use simple approach with text drawn later:
                mov1(:,:,i) = frameImg;
            end
            
            %% Custom colormap: #5770FF → #7F32BD → #FB8808 → #D84420 → #690002
            c1 = [ 87, 112, 255]/255;  % #5770FF
            c2 = [127,  50, 189]/255;  % #7F32BD
            c3 = [251, 136,   8]/255;  % #FB8808
            c4 = [216,  68,  32]/255;  % #D84420
            c5 = [105,   0,   2]/255;  % #690002
            
            nSteps = 500; nSeg = 4; stepsPerSeg = round(nSteps/nSeg);
            r = [linspace(c1(1),c2(1),stepsPerSeg), linspace(c2(1),c3(1),stepsPerSeg), linspace(c3(1),c4(1),stepsPerSeg), linspace(c4(1),c5(1),stepsPerSeg)];
            g = [linspace(c1(2),c2(2),stepsPerSeg), linspace(c2(2),c3(2),stepsPerSeg), linspace(c3(2),c4(2),stepsPerSeg), linspace(c4(2),c5(2),stepsPerSeg)];
            b = [linspace(c1(3),c2(3),stepsPerSeg), linspace(c2(3),c3(3),stepsPerSeg), linspace(c3(3),c4(3),stepsPerSeg), linspace(c4(3),c5(3),stepsPerSeg)];
            cmap = [r(:), g(:), b(:)];
            cmap = cmap(1:nSteps,:);
            
            %% Precompute per-trace metrics: InSegment majority, D (um^2/s) and alpha
            nTraces1 = size(traces1,1);
            isInSegment = false(nTraces1,1);
            Dvals = NaN(nTraces1,1);
            alphas = NaN(nTraces1,1);
            
            % parameters for MSD calculation and fit
            expTime = 1 / frameRate;   % seconds per frame (exposure time between frames)
            minPtsForFit = 3;          % need at least 3 lag points to fit
            maxLagFraction = 0.5;      % use up to 50% of trace length for MSD lags
            
            for j = 1:nTraces1
                currTrace = traces1{j,1}; % table expected with columns 'row','col','z', ... 'InSegment'
                % --- determine InSegment majority ---
                if ismember('InSegment', currTrace.Properties.VariableNames)
                    segVec = currTrace.InSegment;
                    % Safety: ensure numeric and logical values (0/1)
                    segNumeric = double(segVec);
                    % handle NaNs by treating them as 0
                    segNumeric(isnan(segNumeric)) = 0;
                    isInSegment(j) = mean(segNumeric) > 0.5;
                else
                    isInSegment(j) = false;
                end
            
                % --- compute MSD (2D) and fit for D and alpha ---
                % Build coordinate matrix: [x, y] in nanometers (as in your other script)
                if all(ismember({'col','row'}, currTrace.Properties.VariableNames))
                    coord_nm = [currTrace.col, currTrace.row];
                else
                    % fallback if differently named
                    try
                        coord_nm = [currTrace.x, currTrace.y];
                    catch
                        coord_nm = [];
                    end
                end
            
                if isempty(coord_nm) || size(coord_nm,1) < 4
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % center coordinates (like your original code)
                CM = mean(coord_nm,1);
                coord_nm_centered = coord_nm - CM;   % still in nm
            
                % convert to um for MSD (your other script used /1e3)
                coord_um = coord_nm_centered / 1e3;  % um
            
                % compute MSD for lags 1..maxLag (2D)
                nPts = size(coord_um,1);
                maxLag = max(1, floor(nPts * maxLagFraction));
                msd = NaN(maxLag,1);
                for lag = 1:maxLag
                    diffs = coord_um(1+lag:end,:) - coord_um(1:end-lag,:);
                    sqdist = sum(diffs.^2,2); % 2D
                    msd(lag) = mean(sqdist);
                end
            
                % prepare tau (s)
                tau = (1:maxLag)' * expTime;
            
                % require sufficient valid msd points
                validIdx = ~isnan(msd) & msd>0;
                if sum(validIdx) < minPtsForFit
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % Use first part of the msd (valid points) for fit; prefer not to include highest noisy lags:
                useIdx = find(validIdx);
                % choose up to first Nfit points (or all valid if small):
                Nfit = min(length(useIdx), max(3, min(10, length(useIdx))));
                fitIdx = useIdx(1:Nfit);
            
                % log-log fit: log(MSD) = log(4 D) + alpha * log(t)
                x = log(tau(fitIdx));
                y = log(msd(fitIdx));
                % linear fit
                p = polyfit(x,y,1);
                alpha_est = p(1);
                logC = p(2);
                % C = 4 D  => D = exp(logC)/4
                D_est = exp(logC) / 4;
            
                % Validate results (positive D)
                if D_est <= 0 || ~isfinite(alpha_est)
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                else
                    Dvals(j) = D_est;     % um^2 / s
                    alphas(j) = alpha_est;
                end
            end
            
            % Some summary stats (used later for scaling)
            D_min = nanmin(Dvals); D_max = nanmax(Dvals);
            alpha_min = nanmin(alphas); alpha_max = nanmax(alphas);
            
            % avoid degenerate ranges
            if isempty(D_min) || ~isfinite(D_min), D_min = 0; end
            if isempty(D_max) || ~isfinite(D_max), D_max = 1; end
            if D_min == D_max, D_max = D_min + eps; end
            if isempty(alpha_min) || ~isfinite(alpha_min), alpha_min = 0; end
            if isempty(alpha_max) || ~isfinite(alpha_max), alpha_max = 1; end
            if alpha_min == alpha_max, alpha_max = alpha_min + eps; end
            
            % helper mapping functions
            ts2color = @(val, vmin, vmax) max(1, min(nSteps, round((val - vmin)/(vmax - vmin)*(nSteps-1))+1 ));
            
            %% Common figure setup function (we will reuse for each video)
            createFigureAndHandles = @(cbarLabel) deal(...
                figure('Position',[100 100 1400 600],'Renderer','opengl'), ...
                tiledlayout(1,21,'TileSpacing','compact','Padding','compact'));
        
            %% Pre-create line objects per trace (use one set per video to avoid flicker)
            % We'll create 3 separate figures and sets of line handles.
            
            % ------------------ 1) SEGMENT video ------------------
            [FigSeg, tLayout] = createFigureAndHandles('Segment');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax2 = nexttile(tLayout, [1 1]); hold(ax2,'on'); % colorbar/legend
            set(ax1, 'YDir', 'reverse')
            img_h1 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            % create line handles
            lineH_seg = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_seg(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % create a small legend in ax2: two colors (in segment / out)
            cla(ax2);
            legendColors = [c1; c4]; % in-seg = c1 (blue), out = gray
            legendImgs = reshape(legendColors, [1, 2, 3]); % small image
            imagesc(ax2, [0 1], [0 1], permute(repmat(legendColors,1,1),[1 3 2])); % quick hack
            axis(ax2,'off');
            text(ax2,0.5,0.7,'In segment','HorizontalAlignment','center','FontWeight','bold');
            text(ax2,0.5,0.3,'Out of segment','HorizontalAlignment','center');
            
            % timer text top-right
            timeText_seg = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            % Video writer
            v_seg = VideoWriter(mp4File_seg,'MPEG-4'); v_seg.FrameRate = frameRate; open(v_seg);
            
            % Loop frames and update
            for i = 1:nFrames
                set(img_h1,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_seg(j).XData = x;
                            lineH_seg(j).YData = y;
                            if isInSegment(j)
                                lineH_seg(j).Color = legendColors(1,:);
                            else
                                lineH_seg(j).Color = legendColors(2,:);
                            end
                        end
                    end
                end
                % update timer (seconds)
                secNow = round((i-1) * expTime);
                set(timeText_seg, 'String', sprintf('%d s', secNow));
                % draw scalebar text '10 µm' above bar:
                % place relative to scalebar location
                txtX = colStart + 15; txtY = rowStart - 15;
                % remove previous text if exists (simple approach: store in appdata or create once outside loop)
                % For simplicity create a persistent text object:
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigSeg);
                writeVideo(v_seg, frame);
            end
            close(v_seg);
            close(FigSeg);
            
            % ------------------ 2) DIFFUSION video ------------------
            [FigDiff, tLayout] = createFigureAndHandles('Diffusion');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax3 = nexttile(tLayout, [1 1]); hold(ax3,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h2 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            lineH_diff = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_diff(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Draw colorbar for diffusion (map D_min -> D_max onto cmap)
            imagesc(ax3, [0 1], linspace(DMin,DMax,nSteps), permute(cmap,[1 3 2]));
            set(ax3,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax3,[DMin, DMax]); yticklabels(ax3,{sprintf('%.2g',DMin), sprintf('%.2g',DMax)});
            ylabel(ax3,'Diffusion D (µm^2/s)','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_diff = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_diff = VideoWriter(mp4File_diff,'MPEG-4'); v_diff.FrameRate = frameRate; open(v_diff);
            
            for i = 1:nFrames
                set(img_h2,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_diff(j).XData = x;
                            lineH_diff(j).YData = y;
                            Dval = Dvals(j);
                            if isnan(Dval)
                                lineH_diff(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(Dval, DMin, DMax);
                                lineH_diff(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_diff, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigDiff);
                writeVideo(v_diff, frame);
            end
            close(v_diff);
            close(FigDiff);
            
            % ------------------ 3) ALPHA video ------------------
            [FigAlpha, tLayout] = createFigureAndHandles('Alpha');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax4 = nexttile(tLayout, [1 1]); hold(ax4,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h3 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            lineH_alpha = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_alpha(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Colorbar for alpha
            imagesc(ax4, [0 1], linspace(Amin,Amax,nSteps), permute(cmap,[1 3 2]));
            set(ax4,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax4,[Amin, Amax]); yticklabels(ax4,{sprintf('%.2g',Amin), sprintf('%.2g',Amax)});
            ylabel(ax4,'Anomalous exponent \alpha','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_alpha = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_alpha = VideoWriter(mp4File_alpha,'MPEG-4'); v_alpha.FrameRate = frameRate; open(v_alpha);
            
            for i = 1:nFrames
                set(img_h3,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_alpha(j).XData = x;
                            lineH_alpha(j).YData = y;
                            aval = alphas(j);
                            if isnan(aval)
                                lineH_alpha(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(aval, Amin, Amax);
                                lineH_alpha(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_alpha, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigAlpha);
                writeVideo(v_alpha, frame);
            end
            close(v_alpha);
            close(FigAlpha);
            
            disp('All three MP4s saved:');
            disp(mp4File_seg);
            disp(mp4File_diff);
            disp(mp4File_alpha);
        catch
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Background Lysosomes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
        try
            %% Paths for dataset 1 (left)
            trackFile1  = fullfile(mainFolder1,'TracesWMask.mat');
        
            %% Video output names
            mp4File_seg = fullfile(mainFolder1,'TraceMovie_Segment_bgLysosomes.mp4');
            mp4File_diff = fullfile(mainFolder1,'TraceMovie_Diffusion_bgLysosomes.mp4');
            mp4File_alpha = fullfile(mainFolder1,'TraceMovie_Alpha_bgLysosomes.mp4');
            
            %% Load track data
            load(trackFile1)%,'traces3D'); 
            traces1 = Traces3D;
            
            %% Load raw videos
            Folder = dir(mainFolder1);
            idx = find(and(contains({Folder.name}, '.HIS'), contains({Folder.name}, 'Cam1')) == 1);
            videoFile1 = fullfile(mainFolder1, Folder(idx).name);
            mov1 = Load.Movie.his.getFrame(videoFile1, framesToLoad);
            [imgH, imgW, nFrames] = size(mov1);
            nY = imgH; nX = imgW;
            
            % Parameters for scalebar
            scaleBarLength_nm = 10000; % 10 µm
            scaleBarLength_px = round(scaleBarLength_nm / pxSize);
            barHeight = 10;   % pixels thick
            margin = 10;     % pixels from bottom-left
            
            % add scalebar to each frame (bright line)
            for i = 1:nFrames
                frameImg = mov1(:,:,i);  % uint16 frame
                brightValue = max(frameImg(:)); % maximum for uint16 image
                % Bottom-right coordinates (y increases downward; top-left is (1,1))
                rowStart = nY - margin - barHeight;
                rowEnd   = rowStart + barHeight - 1;
                colEnd   = nX - margin;
                colStart = colEnd - scaleBarLength_px + 1;
                if rowStart < 1, rowStart = 1; end
                if colStart < 1, colStart = 1; end
            
                frameImg(rowStart:rowEnd, colStart:colEnd) = brightValue;
            
                % Insert '10 µm' text above bar — use simple approach with text drawn later:
                mov1(:,:,i) = frameImg;
            end
            
            %% Custom colormap: #5770FF → #7F32BD → #FB8808 → #D84420 → #690002
            c1 = [ 87, 112, 255]/255;  % #5770FF
            c2 = [127,  50, 189]/255;  % #7F32BD
            c3 = [251, 136,   8]/255;  % #FB8808
            c4 = [216,  68,  32]/255;  % #D84420
            c5 = [105,   0,   2]/255;  % #690002
            
            nSteps = 500; nSeg = 4; stepsPerSeg = round(nSteps/nSeg);
            r = [linspace(c1(1),c2(1),stepsPerSeg), linspace(c2(1),c3(1),stepsPerSeg), linspace(c3(1),c4(1),stepsPerSeg), linspace(c4(1),c5(1),stepsPerSeg)];
            g = [linspace(c1(2),c2(2),stepsPerSeg), linspace(c2(2),c3(2),stepsPerSeg), linspace(c3(2),c4(2),stepsPerSeg), linspace(c4(2),c5(2),stepsPerSeg)];
            b = [linspace(c1(3),c2(3),stepsPerSeg), linspace(c2(3),c3(3),stepsPerSeg), linspace(c3(3),c4(3),stepsPerSeg), linspace(c4(3),c5(3),stepsPerSeg)];
            cmap = [r(:), g(:), b(:)];
            cmap = cmap(1:nSteps,:);
            
            %% Precompute per-trace metrics: InSegment majority, D (um^2/s) and alpha
            nTraces1 = size(traces1,1);
            isInSegment = false(nTraces1,1);
            Dvals = NaN(nTraces1,1);
            alphas = NaN(nTraces1,1);
            
            % parameters for MSD calculation and fit
            expTime = 1 / frameRate;   % seconds per frame (exposure time between frames)
            minPtsForFit = 3;          % need at least 3 lag points to fit
            maxLagFraction = 0.5;      % use up to 50% of trace length for MSD lags
            
            for j = 1:nTraces1
                currTrace = traces1{j,1}; % table expected with columns 'row','col','z', ... 'InSegment'
                % --- determine InSegment majority ---
                if ismember('InSegment', currTrace.Properties.VariableNames)
                    segVec = currTrace.InSegment;
                    % Safety: ensure numeric and logical values (0/1)
                    segNumeric = double(segVec);
                    % handle NaNs by treating them as 0
                    segNumeric(isnan(segNumeric)) = 0;
                    isInSegment(j) = mean(segNumeric) > 0.5;
                else
                    isInSegment(j) = false;
                end
            
                % --- compute MSD (2D) and fit for D and alpha ---
                % Build coordinate matrix: [x, y] in nanometers (as in your other script)
                if all(ismember({'col','row'}, currTrace.Properties.VariableNames))
                    coord_nm = [currTrace.col, currTrace.row];
                else
                    % fallback if differently named
                    try
                        coord_nm = [currTrace.x, currTrace.y];
                    catch
                        coord_nm = [];
                    end
                end
            
                if isempty(coord_nm) || size(coord_nm,1) < 4
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % center coordinates (like your original code)
                CM = mean(coord_nm,1);
                coord_nm_centered = coord_nm - CM;   % still in nm
            
                % convert to um for MSD (your other script used /1e3)
                coord_um = coord_nm_centered / 1e3;  % um
            
                % compute MSD for lags 1..maxLag (2D)
                nPts = size(coord_um,1);
                maxLag = max(1, floor(nPts * maxLagFraction));
                msd = NaN(maxLag,1);
                for lag = 1:maxLag
                    diffs = coord_um(1+lag:end,:) - coord_um(1:end-lag,:);
                    sqdist = sum(diffs.^2,2); % 2D
                    msd(lag) = mean(sqdist);
                end
            
                % prepare tau (s)
                tau = (1:maxLag)' * expTime;
            
                % require sufficient valid msd points
                validIdx = ~isnan(msd) & msd>0;
                if sum(validIdx) < minPtsForFit
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % Use first part of the msd (valid points) for fit; prefer not to include highest noisy lags:
                useIdx = find(validIdx);
                % choose up to first Nfit points (or all valid if small):
                Nfit = min(length(useIdx), max(3, min(10, length(useIdx))));
                fitIdx = useIdx(1:Nfit);
            
                % log-log fit: log(MSD) = log(4 D) + alpha * log(t)
                x = log(tau(fitIdx));
                y = log(msd(fitIdx));
                % linear fit
                p = polyfit(x,y,1);
                alpha_est = p(1);
                logC = p(2);
                % C = 4 D  => D = exp(logC)/4
                D_est = exp(logC) / 4;
            
                % Validate results (positive D)
                if D_est <= 0 || ~isfinite(alpha_est)
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                else
                    Dvals(j) = D_est;     % um^2 / s
                    alphas(j) = alpha_est;
                end
            end
            
            % Some summary stats (used later for scaling)
            D_min = nanmin(Dvals); D_max = nanmax(Dvals);
            alpha_min = nanmin(alphas); alpha_max = nanmax(alphas);
            
            % avoid degenerate ranges
            if isempty(D_min) || ~isfinite(D_min), D_min = 0; end
            if isempty(D_max) || ~isfinite(D_max), D_max = 1; end
            if D_min == D_max, D_max = D_min + eps; end
            if isempty(alpha_min) || ~isfinite(alpha_min), alpha_min = 0; end
            if isempty(alpha_max) || ~isfinite(alpha_max), alpha_max = 1; end
            if alpha_min == alpha_max, alpha_max = alpha_min + eps; end
            
            % helper mapping functions
            ts2color = @(val, vmin, vmax) max(1, min(nSteps, round((val - vmin)/(vmax - vmin)*(nSteps-1))+1 ));
            
            %% Common figure setup function (we will reuse for each video)
            createFigureAndHandles = @(cbarLabel) deal(...
                figure('Position',[100 100 1400 600],'Renderer','opengl'), ...
                tiledlayout(1,21,'TileSpacing','compact','Padding','compact'));
        
            %% Pre-create line objects per trace (use one set per video to avoid flicker)
            % We'll create 3 separate figures and sets of line handles.
            
            % ------------------ 1) SEGMENT video ------------------
            [FigSeg, tLayout] = createFigureAndHandles('Segment');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax2 = nexttile(tLayout, [1 1]); hold(ax2,'on'); % colorbar/legend
            set(ax1, 'YDir', 'reverse')
            img_h1 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            % create line handles
            lineH_seg = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_seg(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % create a small legend in ax2: two colors (in segment / out)
            cla(ax2);
            legendColors = [c1; c4]; % in-seg = c1 (blue), out = gray
            legendImgs = reshape(legendColors, [1, 2, 3]); % small image
            imagesc(ax2, [0 1], [0 1], permute(repmat(legendColors,1,1),[1 3 2])); % quick hack
            axis(ax2,'off');
            text(ax2,0.5,0.7,'In segment','HorizontalAlignment','center','FontWeight','bold');
            text(ax2,0.5,0.3,'Out of segment','HorizontalAlignment','center');
            
            % timer text top-right
            timeText_seg = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            % Video writer
            v_seg = VideoWriter(mp4File_seg,'MPEG-4'); v_seg.FrameRate = frameRate; open(v_seg);
            
            % Loop frames and update
            for i = 1:nFrames
                set(img_h1,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_seg(j).XData = x;
                            lineH_seg(j).YData = y;
                            if isInSegment(j)
                                lineH_seg(j).Color = legendColors(1,:);
                            else
                                lineH_seg(j).Color = legendColors(2,:);
                            end
                        end
                    end
                end
                % update timer (seconds)
                secNow = round((i-1) * expTime);
                set(timeText_seg, 'String', sprintf('%d s', secNow));
                % draw scalebar text '10 µm' above bar:
                % place relative to scalebar location
                txtX = colStart + 15; txtY = rowStart - 15;
                % remove previous text if exists (simple approach: store in appdata or create once outside loop)
                % For simplicity create a persistent text object:
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigSeg);
                writeVideo(v_seg, frame);
            end
            close(v_seg);
            close(FigSeg);
            
            % ------------------ 2) DIFFUSION video ------------------
            [FigDiff, tLayout] = createFigureAndHandles('Diffusion');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax3 = nexttile(tLayout, [1 1]); hold(ax3,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h2 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            lineH_diff = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_diff(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Draw colorbar for diffusion (map D_min -> D_max onto cmap)
            imagesc(ax3, [0 1], linspace(DMin,DMax,nSteps), permute(cmap,[1 3 2]));
            set(ax3,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax3,[DMin, DMax]); yticklabels(ax3,{sprintf('%.2g',DMin), sprintf('%.2g',DMax)});
            ylabel(ax3,'Diffusion D (µm^2/s)','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_diff = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_diff = VideoWriter(mp4File_diff,'MPEG-4'); v_diff.FrameRate = frameRate; open(v_diff);
            
            for i = 1:nFrames
                set(img_h2,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_diff(j).XData = x;
                            lineH_diff(j).YData = y;
                            Dval = Dvals(j);
                            if isnan(Dval)
                                lineH_diff(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(Dval, DMin, DMax);
                                lineH_diff(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_diff, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigDiff);
                writeVideo(v_diff, frame);
            end
            close(v_diff);
            close(FigDiff);
            
            % ------------------ 3) ALPHA video ------------------
            [FigAlpha, tLayout] = createFigureAndHandles('Alpha');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax4 = nexttile(tLayout, [1 1]); hold(ax4,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h3 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            
            lineH_alpha = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_alpha(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Colorbar for alpha
            imagesc(ax4, [0 1], linspace(Amin,Amax,nSteps), permute(cmap,[1 3 2]));
            set(ax4,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax4,[Amin, Amax]); yticklabels(ax4,{sprintf('%.2g',Amin), sprintf('%.2g',Amax)});
            ylabel(ax4,'Anomalous exponent \alpha','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_alpha = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_alpha = VideoWriter(mp4File_alpha,'MPEG-4'); v_alpha.FrameRate = frameRate; open(v_alpha);
            
            for i = 1:nFrames
                set(img_h3,'CData', mov1(:,:,i));
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_alpha(j).XData = x;
                            lineH_alpha(j).YData = y;
                            aval = alphas(j);
                            if isnan(aval)
                                lineH_alpha(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(aval, Amin, Amax);
                                lineH_alpha(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_alpha, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigAlpha);
                writeVideo(v_alpha, frame);
            end
            close(v_alpha);
            close(FigAlpha);
            
            disp('All three MP4s saved:');
            disp(mp4File_seg);
            disp(mp4File_diff);
            disp(mp4File_alpha);
        catch
        end
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Background Particles +
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% transparent lysosomes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
    
        try
    
            maskFile = fullfile(mainFolder1,'SegmentMovie','SegmentMovie.mat');
            load(maskFile,'Mask'); % loads cell array Mask{nFrames}, each 512x512 double (0/1)
        
            % Transparency and color for mask overlay
            maskColor = cat(3, ones(imgH,imgW), zeros(imgH,imgW), zeros(imgH,imgW)); % red
            alphaVal = 0.3;
        
            %% Paths for dataset 1 (left)
            trackFile1  = fullfile(mainFolder1,'TracesWMask.mat');
        
            %% Video output names
            mp4File_seg = fullfile(mainFolder1,'TraceMovie_Segment_bgParticles_maskLysosomes.mp4');
            mp4File_diff = fullfile(mainFolder1,'TraceMovie_Diffusion_bgParticles_maskLysosomes.mp4');
            mp4File_alpha = fullfile(mainFolder1,'TraceMovie_Alpha_bgParticles_maskLysosomes.mp4');
            
            %% Load track data
            load(trackFile1)%,'traces3D'); 
            traces1 = Traces3D;
            
            %% Load raw videos
            Folder = dir(mainFolder1);
            idx = find(and(contains({Folder.name}, '.HIS'), contains({Folder.name}, 'Cam1')) == 1);
            videoFile1 = fullfile(mainFolder1, Folder(idx).name);
            mov1 = Load.Movie.his.getFrame(videoFile1, framesToLoad);
            [imgH, imgW, nFrames] = size(mov1);
            nY = imgH; nX = imgW;
            
            % Parameters for scalebar
            scaleBarLength_nm = 10000; % 10 µm
            scaleBarLength_px = round(scaleBarLength_nm / pxSize);
            barHeight = 10;   % pixels thick
            margin = 10;     % pixels from bottom-left
            
            % add scalebar to each frame (bright line)
            for i = 1:nFrames
                frameImg = mov1(:,:,i);  % uint16 frame
                brightValue = max(frameImg(:)); % maximum for uint16 image
                % Bottom-right coordinates (y increases downward; top-left is (1,1))
                rowStart = nY - margin - barHeight;
                rowEnd   = rowStart + barHeight - 1;
                colEnd   = nX - margin;
                colStart = colEnd - scaleBarLength_px + 1;
                if rowStart < 1, rowStart = 1; end
                if colStart < 1, colStart = 1; end
            
                frameImg(rowStart:rowEnd, colStart:colEnd) = brightValue;
            
                % Insert '10 µm' text above bar — use simple approach with text drawn later:
                mov1(:,:,i) = frameImg;
            end
            
            %% Custom colormap: #5770FF → #7F32BD → #FB8808 → #D84420 → #690002
            c1 = [ 87, 112, 255]/255;  % #5770FF
            c2 = [127,  50, 189]/255;  % #7F32BD
            c3 = [251, 136,   8]/255;  % #FB8808
            c4 = [216,  68,  32]/255;  % #D84420
            c5 = [105,   0,   2]/255;  % #690002
            
            nSteps = 500; nSeg = 4; stepsPerSeg = round(nSteps/nSeg);
            r = [linspace(c1(1),c2(1),stepsPerSeg), linspace(c2(1),c3(1),stepsPerSeg), linspace(c3(1),c4(1),stepsPerSeg), linspace(c4(1),c5(1),stepsPerSeg)];
            g = [linspace(c1(2),c2(2),stepsPerSeg), linspace(c2(2),c3(2),stepsPerSeg), linspace(c3(2),c4(2),stepsPerSeg), linspace(c4(2),c5(2),stepsPerSeg)];
            b = [linspace(c1(3),c2(3),stepsPerSeg), linspace(c2(3),c3(3),stepsPerSeg), linspace(c3(3),c4(3),stepsPerSeg), linspace(c4(3),c5(3),stepsPerSeg)];
            cmap = [r(:), g(:), b(:)];
            cmap = cmap(1:nSteps,:);
            
            %% Precompute per-trace metrics: InSegment majority, D (um^2/s) and alpha
            nTraces1 = size(traces1,1);
            isInSegment = false(nTraces1,1);
            Dvals = NaN(nTraces1,1);
            alphas = NaN(nTraces1,1);
            
            % parameters for MSD calculation and fit
            expTime = 1 / frameRate;   % seconds per frame (exposure time between frames)
            minPtsForFit = 3;          % need at least 3 lag points to fit
            maxLagFraction = 0.5;      % use up to 50% of trace length for MSD lags
            
            for j = 1:nTraces1
                currTrace = traces1{j,1}; % table expected with columns 'row','col','z', ... 'InSegment'
                % --- determine InSegment majority ---
                if ismember('InSegment', currTrace.Properties.VariableNames)
                    segVec = currTrace.InSegment;
                    % Safety: ensure numeric and logical values (0/1)
                    segNumeric = double(segVec);
                    % handle NaNs by treating them as 0
                    segNumeric(isnan(segNumeric)) = 0;
                    isInSegment(j) = mean(segNumeric) > 0.5;
                else
                    isInSegment(j) = false;
                end
            
                % --- compute MSD (2D) and fit for D and alpha ---
                % Build coordinate matrix: [x, y] in nanometers (as in your other script)
                if all(ismember({'col','row'}, currTrace.Properties.VariableNames))
                    coord_nm = [currTrace.col, currTrace.row];
                else
                    % fallback if differently named
                    try
                        coord_nm = [currTrace.x, currTrace.y];
                    catch
                        coord_nm = [];
                    end
                end
            
                if isempty(coord_nm) || size(coord_nm,1) < 4
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % center coordinates (like your original code)
                CM = mean(coord_nm,1);
                coord_nm_centered = coord_nm - CM;   % still in nm
            
                % convert to um for MSD (your other script used /1e3)
                coord_um = coord_nm_centered / 1e3;  % um
            
                % compute MSD for lags 1..maxLag (2D)
                nPts = size(coord_um,1);
                maxLag = max(1, floor(nPts * maxLagFraction));
                msd = NaN(maxLag,1);
                for lag = 1:maxLag
                    diffs = coord_um(1+lag:end,:) - coord_um(1:end-lag,:);
                    sqdist = sum(diffs.^2,2); % 2D
                    msd(lag) = mean(sqdist);
                end
            
                % prepare tau (s)
                tau = (1:maxLag)' * expTime;
            
                % require sufficient valid msd points
                validIdx = ~isnan(msd) & msd>0;
                if sum(validIdx) < minPtsForFit
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                    continue
                end
            
                % Use first part of the msd (valid points) for fit; prefer not to include highest noisy lags:
                useIdx = find(validIdx);
                % choose up to first Nfit points (or all valid if small):
                Nfit = min(length(useIdx), max(3, min(10, length(useIdx))));
                fitIdx = useIdx(1:Nfit);
            
                % log-log fit: log(MSD) = log(4 D) + alpha * log(t)
                x = log(tau(fitIdx));
                y = log(msd(fitIdx));
                % linear fit
                p = polyfit(x,y,1);
                alpha_est = p(1);
                logC = p(2);
                % C = 4 D  => D = exp(logC)/4
                D_est = exp(logC) / 4;
            
                % Validate results (positive D)
                if D_est <= 0 || ~isfinite(alpha_est)
                    Dvals(j) = NaN;
                    alphas(j) = NaN;
                else
                    Dvals(j) = D_est;     % um^2 / s
                    alphas(j) = alpha_est;
                end
            end
            
            % Some summary stats (used later for scaling)
            D_min = nanmin(Dvals); D_max = nanmax(Dvals);
            alpha_min = nanmin(alphas); alpha_max = nanmax(alphas);
            
            % avoid degenerate ranges
            if isempty(D_min) || ~isfinite(D_min), D_min = 0; end
            if isempty(D_max) || ~isfinite(D_max), D_max = 1; end
            if D_min == D_max, D_max = D_min + eps; end
            if isempty(alpha_min) || ~isfinite(alpha_min), alpha_min = 0; end
            if isempty(alpha_max) || ~isfinite(alpha_max), alpha_max = 1; end
            if alpha_min == alpha_max, alpha_max = alpha_min + eps; end
            
            % helper mapping functions
            ts2color = @(val, vmin, vmax) max(1, min(nSteps, round((val - vmin)/(vmax - vmin)*(nSteps-1))+1 ));
            
            %% Common figure setup function (we will reuse for each video)
            createFigureAndHandles = @(cbarLabel) deal(...
                figure('Position',[100 100 1400 600],'Renderer','opengl'), ...
                tiledlayout(1,21,'TileSpacing','compact','Padding','compact'));
        
            %% Pre-create line objects per trace (use one set per video to avoid flicker)
            % We'll create 3 separate figures and sets of line handles.
            
            % ------------------ 1) SEGMENT video ------------------
            [FigSeg, tLayout] = createFigureAndHandles('Segment');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax2 = nexttile(tLayout, [1 1]); hold(ax2,'on'); % colorbar/legend
            set(ax1, 'YDir', 'reverse')
            img_h1 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            hold(ax1,'on');
            mask_h1 = imshow(maskColor, 'Parent', ax1);
            set(mask_h1, 'AlphaData', alphaVal * Mask{1});
            
            % create line handles
            lineH_seg = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_seg(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % create a small legend in ax2: two colors (in segment / out)
            cla(ax2);
            legendColors = [c1; c4]; % in-seg = c1 (blue), out = gray
            legendImgs = reshape(legendColors, [1, 2, 3]); % small image
            imagesc(ax2, [0 1], [0 1], permute(repmat(legendColors,1,1),[1 3 2])); % quick hack
            axis(ax2,'off');
            text(ax2,0.5,0.7,'In segment','HorizontalAlignment','center','FontWeight','bold');
            text(ax2,0.5,0.3,'Out of segment','HorizontalAlignment','center');
            
            % timer text top-right
            timeText_seg = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            % Video writer
            v_seg = VideoWriter(mp4File_seg,'MPEG-4'); v_seg.FrameRate = frameRate; open(v_seg);
            
            % Loop frames and update
            for i = 1:nFrames
                set(img_h1,'CData', mov1(:,:,i));
                set(mask_h1,'AlphaData', alphaVal * Mask{i});
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_seg(j).XData = x;
                            lineH_seg(j).YData = y;
                            if isInSegment(j)
                                lineH_seg(j).Color = legendColors(1,:);
                            else
                                lineH_seg(j).Color = legendColors(2,:);
                            end
                        end
                    end
                end
                % update timer (seconds)
                secNow = round((i-1) * expTime);
                set(timeText_seg, 'String', sprintf('%d s', secNow));
                % draw scalebar text '10 µm' above bar:
                % place relative to scalebar location
                txtX = colStart + 15; txtY = rowStart - 15;
                % remove previous text if exists (simple approach: store in appdata or create once outside loop)
                % For simplicity create a persistent text object:
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigSeg);
                writeVideo(v_seg, frame);
            end
            close(v_seg);
            close(FigSeg);
            
            % ------------------ 2) DIFFUSION video ------------------
            [FigDiff, tLayout] = createFigureAndHandles('Diffusion');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax3 = nexttile(tLayout, [1 1]); hold(ax3,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h2 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            hold(ax1,'on');
            mask_h2 = imshow(maskColor, 'Parent', ax1);
            set(mask_h2, 'AlphaData', alphaVal * Mask{1});
            
            lineH_diff = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_diff(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Draw colorbar for diffusion (map D_min -> D_max onto cmap)
            imagesc(ax3, [0 1], linspace(DMin,DMax,nSteps), permute(cmap,[1 3 2]));
            set(ax3,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax3,[DMin, DMax]); yticklabels(ax3,{sprintf('%.2g',DMin), sprintf('%.2g',DMax)});
            ylabel(ax3,'Diffusion D (µm^2/s)','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_diff = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_diff = VideoWriter(mp4File_diff,'MPEG-4'); v_diff.FrameRate = frameRate; open(v_diff);
            
            for i = 1:nFrames
                set(img_h2,'CData', mov1(:,:,i));
                set(mask_h2,'AlphaData', alphaVal * Mask{i});
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_diff(j).XData = x;
                            lineH_diff(j).YData = y;
                            Dval = Dvals(j);
                            if isnan(Dval)
                                lineH_diff(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(Dval, DMin, DMax);
                                lineH_diff(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_diff, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigDiff);
                writeVideo(v_diff, frame);
            end
            close(v_diff);
            close(FigDiff);
            
            % ------------------ 3) ALPHA video ------------------
            [FigAlpha, tLayout] = createFigureAndHandles('Alpha');
            ax1 = nexttile(tLayout, [1 20]); hold(ax1,'on');
            ax4 = nexttile(tLayout, [1 1]); hold(ax4,'on'); % colorbar
            set(ax1, 'YDir', 'reverse')
            img_h3 = imagesc(ax1, mov1(:,:,1)); colormap(ax1, gray); axis(ax1,'image','off');
            hold(ax1,'on');
            mask_h3 = imshow(maskColor, 'Parent', ax1);
            set(mask_h3, 'AlphaData', alphaVal * Mask{1});
            
            lineH_alpha = gobjects(nTraces1,1);
            for j = 1:nTraces1
                lineH_alpha(j) = plot(ax1, NaN, NaN, 'LineWidth', 2);
            end
            title(ax1, CellName, 'FontSize', 12);
            
            % Colorbar for alpha
            imagesc(ax4, [0 1], linspace(Amin,Amax,nSteps), permute(cmap,[1 3 2]));
            set(ax4,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',12,'FontWeight','bold');
            yticks(ax4,[Amin, Amax]); yticklabels(ax4,{sprintf('%.2g',Amin), sprintf('%.2g',Amax)});
            ylabel(ax4,'Anomalous exponent \alpha','FontSize',12,'FontWeight','bold');
            
            % timer text top-right
            timeText_alpha = text(ax1, nX - margin - 40, margin + 10, '0 s', 'Color','white','FontSize',14,'FontWeight','bold', 'HorizontalAlignment','right');
            
            v_alpha = VideoWriter(mp4File_alpha,'MPEG-4'); v_alpha.FrameRate = frameRate; open(v_alpha);
            
            for i = 1:nFrames
                set(img_h3,'CData', mov1(:,:,i));
                set(mask_h3,'AlphaData', alphaVal * Mask{i});
                for j = 1:nTraces1
                    currTrace = traces1{j,1};
                    if height(currTrace) > 20
                        idx = find(currTrace.t <= i);
                        if numel(idx) >= 2
                            x = currTrace.col(idx) / pxSize;
                            y = currTrace.row(idx) / pxSize;
                            
                            % compute step distances in nm
                            dx_nm = diff(currTrace.col(idx));
                            dy_nm = diff(currTrace.row(idx));
                            stepDist_nm = sqrt(dx_nm.^2 + dy_nm.^2);
                            
                            % insert NaN for steps >500 nm
                            tooBig = stepDist_nm > stepThreshold;
                            x([false; tooBig]) = NaN;
                            y([false; tooBig]) = NaN;
                            lineH_alpha(j).XData = x;
                            lineH_alpha(j).YData = y;
                            aval = alphas(j);
                            if isnan(aval)
                                lineH_alpha(j).Color = [0.7 0.7 0.7]; % fallback gray
                            else
                                cidx = ts2color(aval, Amin, Amax);
                                lineH_alpha(j).Color = cmap(cidx,:);
                            end
                        end
                    end
                end
                % timer
                secNow = round((i-1) * expTime);
                set(timeText_alpha, 'String', sprintf('%d s', secNow));
                if i==1
                    scaletextH = text(ax1, colStart+40, rowStart-5, '10 µm', 'Color','white','FontSize',12,'FontWeight','bold', 'VerticalAlignment','bottom');
                end
            
                frame = captureFrameForVideo(FigAlpha);
                writeVideo(v_alpha, frame);
            end
            close(v_alpha);
            close(FigAlpha);
            
            disp('All three MP4s saved:');
            disp(mp4File_seg);
            disp(mp4File_diff);
            disp(mp4File_alpha);
        catch
        end
    catch
    end
end

%% function to draw a frame (common)
function frame = captureFrameForVideo(figHandle)
    drawnow;
    frame = getframe(figHandle);
end

% %% --- Save last frame grayscale only (unchanged) ---
% FigLastGray = figure('Position',[100 100 1400 600]);
% tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% 
% % LEFT
% axL = nexttile(1);
% imshow(mov1(:,:,end),[]); % grayscale
% title(CellName);
% 
% % Save figure
% saveas(FigLastGray, fullfile(mainFolder1,'LastFrame_Grayscale.png'));
% saveas(FigLastGray, fullfile(mainFolder1,'LastFrame_Grayscale.svg'));
% close(FigLastGray);
% 
% %% --- Save last frame with traces (example using alpha coloring) ---
% FigLastTraces = figure('Position',[100 100 1400 600]);
% tiledlayout(1,21,'TileSpacing','compact','Padding','compact');
% ax1 = nexttile([1 20]); hold on;
% ax3 = nexttile([1 1]); hold on; % colorbar
% set(ax1, 'YDir', 'reverse')
% 
% % LEFT video
% imagesc(ax1, mov1(:,:,end)); colormap(ax1,gray); axis(ax1,'image','off');
% for j = 1:nTraces1
%     currTrace = traces1{j,1};
%     if height(currTrace) > 20
%         idx = find(currTrace.t <= nFrames);
%         if numel(idx)>=2
%             x = currTrace.col(idx)/pxSize;
%             y = currTrace.row(idx)/pxSize;
%             % color by alpha for last-frame preview (fallback to gray)
%             aval = alphas(j);
%             if isnan(aval)
%                 clr = [0.7 0.7 0.7];
%             else
%                 clr = cmap(ts2color(aval,Amin,Amax),:);
%             end
%             plot(ax1, x, y, 'LineWidth', 1, 'Color', clr);
%         end
%     end
% end
% title(ax1,CellName);
% 
% % Colorbar
% imagesc(ax3,[0 1], linspace(alpha_min,alpha_max,nSteps), permute(cmap,[1 3 2]));
% set(ax3,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',14,'FontWeight','bold');
% yticks(ax3,[Amin,Amax]); yticklabels(ax3,{sprintf('%.2g',alpha_min), sprintf('%.2g',alpha_max)});
% ylabel(ax3,'Anomalous exponent \alpha','FontSize',14,'FontWeight','bold');
% 
% % Save figure
% saveas(FigLastTraces, fullfile(mainFolder1,'LastFrame_Traces_Alpha.png'));
% saveas(FigLastTraces, fullfile(mainFolder1,'LastFrame_Traces_Alpha.svg'));
% close(FigLastTraces);