classdef MPPhaseMovie < Core.MPMovie
    %MPPHASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        QPmap
        Cropped
    end
    
    methods
        function obj = MPPhaseMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function getPhaseMovie(obj, q)
            if strcmp(obj.info.runMethod, 'load')
                if exist(append(obj.raw.movInfo.Path, filesep, 'PhaseMovie'))
                    run = 0;
                else
                    run = 1;
                end
            else
                run = 1;
            end

            if run == 1
                f = waitbar(0,'Initializing');
                if strcmp(obj.info.frame2Load, 'all')
                    nFrames = obj.calibrated{1, 1}.nFrames; 
                elseif isa(obj.info.frame2Load, 'double')
                    nFrames = max(obj.info.frame2Load);
                end
    
                s.optics = obj.info.optics;
                s.proc = obj.info.proc;
    
                mkdir(append(obj.raw.movInfo.Path, filesep, 'PhaseMovie'));
    
                n = 1;
                Step = 0;
                ChunkSize = 100;
                for k = 1:ChunkSize:nFrames
                    Step = Step + 1;
                    idx = k:min(k+ChunkSize-1, nFrames);
                    Startidx = idx(1)-1;
                    for i = idx
                        waitbar(n./nFrames,f,append('Calculating phase map ', num2str(n),' out of ', num2str(nFrames)));
                        Stack = obj.getFrame(n, q);
                        [Stack, StartX, StartY] = QP_package.cropXY(Stack);
                        % PadValues = 134;
                        % for l = 1:size(Stack, 3)
                        %     StackPadded(:,:,l) = padarray(Stack(:,:,l), [50 50], PadValues);
                        % end
                        % QPmap(:,:,:,i-Startidx) = QP_package.getQP(StackPadded,s);
                        QPmap(:,:,:,n) = QP_package.getQP(Stack,s);
                        n = n+1;
                    end
                    % QPmap = QPmap(51:end-50, 51:end-50, :,:);
                end
                close(f)
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat');
                save(Filename, 'QPmap');
                obj.QPmap = QPmap;
                QPmap = [];
               
                obj.Cropped.StartX = StartX;
                obj.Cropped.StartY = StartY;
            else
                load(append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat'))
                Stack = obj.getFrame(1, q);
                obj.Cropped.StartX = floor((size(Stack,1) - size(QPmap, 1))./2);
                obj.Cropped.StartY = floor((size(Stack,2) - size(QPmap, 2))./2);
                obj.QPmap = QPmap;
            end
            

            % if strcmp(obj.info.runMethod, 'load')
            %     if exist(append(obj.raw.movInfo.Path, filesep, 'MovieRawData'))
            %         run2 = 0;
            %     else
            %         run2 = 1;
            %     end
            % else
            %     run2 = 1;
            % end
            % 
            % if run2 == 1
            %     f = waitbar(0,'Initializing');
            %     if strcmp(obj.info.frame2Load, 'all')
            %         nFrames = obj.calibrated{1, 1}.nFrames; 
            %     elseif isa(obj.info.frame2Load, 'double')
            %         nFrames = max(obj.info.frame2Load);
            %     end
            % 
            %     mkdir(append(obj.raw.movInfo.Path, filesep, 'MovieRawData'));
            % 
            %     n = 1;
            %     Step = 0;
            %     ChunkSize = 100;
            %     for k = 1:ChunkSize:nFrames
            %         Step = Step + 1;
            %         idx = k:min(k+ChunkSize-1, nFrames);
            %         Startidx = idx(1)-1;
            %         for i = idx
            %             waitbar(n./nFrames,f,append('Calculating raw map ', num2str(n),' out of ', num2str(nFrames)));
            %             StackRaw(:,:,:,i-Startidx) = obj.getFrame(n, q);
            %             n = n+1;
            %         end
            %         Filename = append(obj.raw.movInfo.Path, filesep, 'MovieRawData', filesep, 'MovieRawData', num2str(Step), '.mat');
            %         save(Filename, 'StackRaw');
            %         StackRaw = [];
            %     end
            %     close(f)
            % else
            %     disp('Raw intensity data already saved next to phase')
            % end
        end

        function [MeanPhasePart] = calibrateAlpha(obj, q)
            QPmap = obj.QPmap(:, 25:360,:,:);
            MeanQP = min(QPmap, 4);
            BigFig = figure();
            PlaneInFocus = 1;
            for plane = 1:size(obj.QPmap, 3)
                try
                    PlaneImage = squeeze(QPmap(:,:,plane, :));
                    PlaneImage(PlaneImage > 0) = NaN;
                    MinProjection = min(PlaneImage, [], 3, 'omitnan');
                    [N, edges] = histcounts(MinProjection(:));
    
                    x = MinProjection(:);
                    [f,xi] = ksdensity(x);
                    StartIdx = find(xi > -1.6, 1, 'first');
                    EndIdx = find(xi < 0, 1, 'last');
                    [~,locs] = min(abs(diff(smooth(f(StartIdx:EndIdx)))));   % minima of f
                    threshold = xi(locs(1) + StartIdx);
                    x_small = x(x < threshold);
                    x_big   = x(x >= threshold);
                    [N, edges] = histcounts(x(:));
                    [~, Idx] = min(abs(edges - threshold));
                    SecondDeriv = diff(smooth(smooth(N(1:Idx-1))));
                    figure();
                    subplot(1,2,1)
                    fitty = fit(edges(1:Idx-2)', SecondDeriv, 'gauss1');
                    plot(fitty, edges(1:Idx-2)', SecondDeriv)
                    A = fitty.a1;
                    mu = fitty.b1;
                    sigma = fitty.c1;
                    FWHM = 2*sqrt(2*log(2)) * sigma;
                    title('Gaussian fit on derivative of histogram of psf pixels')
                    subplot(1,2,2)
                    plot(edges(1:Idx-1), smooth(smooth(N(1:Idx-1))))
                    hold on
                    plot(edges(Idx:end-1), smooth(smooth(N(Idx:end))))
                    hold on
                    xline(mu + FWHM./2)
                    title('histogram of minprojection image')
                    xlabel('pixel values minprojection')
                    ylabel('Counts')
    
                    MeanPhasePart(plane,1) = mu + FWHM./2;
                catch
                    MeanPhasePart(plane,1) = NaN;
                end
                %%% find particle centers

                MeanImage = MeanQP(:,:,plane);
                [ParticlePixels, radii, metric] = imfindcircles(-imgaussfilt(MeanImage, 5), [5 10], 'Sensitivity', 0.93);
                fig = figure;
                imagesc(MeanImage);      % show the image
                hold on;
                plot(ParticlePixels(:,1), ParticlePixels(:,2), 'r.', 'MarkerSize', 20);   % red dots on centers
                viscircles(ParticlePixels, radii, 'EdgeColor', 'b');               % optional: draw circles
                ParticlePixels(sum(or(ParticlePixels-5 < 1, ParticlePixels+5 > min(size(MeanImage))),2) == 1, :) = [];
                saveas(fig, append(obj.raw.movInfo.Path, filesep, 'PhaseIm_minProjection_plane', num2str(plane), '.png'));
                save(append(obj.raw.movInfo.Path, filesep, 'PhaseIm_minProjection_plane', num2str(plane), '.mat'), "MeanImage");
                %%% find background pixels
                I_filt = imgaussfilt(MeanImage, 10);
                blobMaskmin = imregionalmin(I_filt,4);
                blobMaskmax = imregionalmax(I_filt,4);
                blobMask = blobMaskmax + blobMaskmin;
                blobMaskDilated = imdilate(blobMask, strel('disk', 25));
                % blobMaskDilated
                % [y, x] = find(blobMaskDilated == 0);
                % idx = randperm(numel(x), size(ParticlePixels, 1));
                % backgroundPoints = [y(idx), x(idx)];

                for frame = 1:size(QPmap,4)
                    Frame = imgaussfilt(QPmap(:,:,plane, frame),5);
                    for pixel = 1:size(ParticlePixels, 1)
                        ParticleROI = Frame(round(ParticlePixels(pixel,1)-5:ParticlePixels(pixel,1)+5),...
                            round(ParticlePixels(pixel,2)-5:ParticlePixels(pixel,2)+5));
                        %ParticleROI = imgaussfilt(ParticleROI, 5);
                        Particle(frame, pixel) = min(ParticleROI, [], 'all'); 
                    end 
                    Bg = QPmap(:,:,plane, frame);
                    Bg(blobMaskDilated == 1) = NaN;
                    Background(frame, 1) = nanmedian(Bg, 'all');
                end

                Fig2 = figure();
                for i = 1:size(Particle,2)
                    plot(Particle(:,i), 'g')
                    hold on
                end
                plot(Background, 'r')
                ylim([-1.6 1.6])

                % subplot(2, size(QPmap,3)./2, plane)
                for i = 1:size(Particle, 2)
                    try
                        h1 =  plot(Particle(:,i), 'g');
                        hold on
                        h2 = plot(Background(:,i), 'r');
                        hold on
                        if i == 1
                            firstGreen = h1;
                            firstRed   = h2;
                        end
                    catch
                    end
                end
                legend([firstGreen firstRed], {'Particle','Background'});

                Particle = [];
                Background = [];
                saveas(Fig2, append(obj.raw.movInfo.Path, filesep, 'ParticlePhase_minProjection_plane', num2str(plane), '.png'));
                % MeanImage = MeanQP(:,:,plane);
                % [ParticlePixels, radii, metric] = imfindcircles(-imgaussfilt(MeanImage, 5), [5 10], 'Sensitivity', 0.93);
                % fig = figure;
                % imagesc(MeanImage);      % show the image
                % hold on;
                % plot(ParticlePixels(:,1), ParticlePixels(:,2), 'r.', 'MarkerSize', 20);   % red dots on centers
                % viscircles(ParticlePixels, radii, 'EdgeColor', 'b');               % optional: draw circles
                % ParticlePixels(sum(or(ParticlePixels-5 < 1, ParticlePixels+5 > min(size(MeanImage))),2) == 1, :) = [];
                % 
                % %%% find background pixels
                % I_filt = imgaussfilt(MeanImage, 10);
                % blobMaskmin = imregionalmin(I_filt,4);
                % blobMaskmax = imregionalmax(I_filt,4);
                % blobMask = blobMaskmax + blobMaskmin;
                % blobMaskDilated = imdilate(blobMask, strel('disk', 25));
                % % blobMaskDilated
                % % [y, x] = find(blobMaskDilated == 0);
                % % idx = randperm(numel(x), size(ParticlePixels, 1));
                % % backgroundPoints = [y(idx), x(idx)];
                % 
                % for frame = 1:size(QPmap,4)
                %     Frame = imgaussfilt(QPmap(:,:,plane, frame),5);
                %     for pixel = 1:size(ParticlePixels, 1)
                %         ParticleROI = Frame(round(ParticlePixels(pixel,1)-5:ParticlePixels(pixel,1)+5),...
                %             round(ParticlePixels(pixel,2)-5:ParticlePixels(pixel,2)+5));
                %         ParticleROI = imgaussfilt(ParticleROI, 5);
                %         Particle(frame, pixel) = min(ParticleROI, [], 'all'); 
                %     end 
                %     Bg = QPmap(:,:,plane, frame);
                % end
                % 
                % % subplot(2, size(QPmap,3)./2, plane)
                % for i = 1:size(Particle, 2)
                %     try
                %         h1 =  plot(Particle(:,i), 'g')
                %         hold on
                %         h2 = plot(Background(:,i), 'r')
                %         hold on
                %         if i == 1
                %             firstGreen = h1;
                %             firstRed   = h2;
                %         end
                %     catch
                %     end
                % end
                % legend([firstGreen firstRed], {'Particle','Background'});
                % 
                % Particle = [];
                % Background = [];
            end
        end
    end
end

