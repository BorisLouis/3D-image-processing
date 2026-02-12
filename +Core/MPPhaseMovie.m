classdef MPPhaseMovie < Core.MPMovie
    %MPPHASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        QPmap
        QPmap_fullPhase
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
                s.optics.dz = mean(abs(diff(obj.calibrated{1, 1}.oRelZPos)));
    
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
                        [QPmap_fullPhase(:,:,n) ,~] = QP_package.getQPFullPhase(Stack,s);
                        [QPmap(:,:,:,n), ~] = QP_package.getQP(Stack,s);
                        n = n+1;
                    end
                    % QPmap = QPmap(51:end-50, 51:end-50, :,:);
                end
                close(f)
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat');
                save(Filename, 'QPmap');
                obj.QPmap = QPmap;
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovieFullPhase.mat');
                save(Filename, 'QPmap_fullPhase');
                obj.QPmap_fullPhase = QPmap_fullPhase;
               
                obj.Cropped.StartX = StartX;
                obj.Cropped.StartY = StartY;
            else
                load(append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat'))
                Stack = obj.getFrame(1, q);
                obj.Cropped.StartX = floor((size(Stack,1) - size(QPmap, 1))./2);
                obj.Cropped.StartY = floor((size(Stack,2) - size(QPmap, 2))./2);
                obj.QPmap = QPmap;
            end
        end

        function calibrateAlpha2(obj, q)
            f = waitbar(0,'Initializing');
            if strcmp(obj.info.frame2Load, 'all')
                nFrames = obj.calibrated{1, 1}.nFrames; 
            elseif isa(obj.info.frame2Load, 'double')
                nFrames = max(obj.info.frame2Load);
            end

            s.optics = obj.info.optics;
            s.optics.dz = 0.1;
            s.optics.alpha = 1.5;
            s.optics.lambda = 0.698;
            s.optics.n = 1.0;
            s.proc = obj.info.proc;

            for Frame = 1:nFrames
                waitbar(Frame./nFrames, f, 'Loading Frames');
                Fullvid(:,:,:,Frame) = obj.getFrame(Frame, q);
            end

            for Plane = 1:size(Fullvid, 3)
                waitbar(Plane./size(Fullvid, 3), f, 'Calculating phasemaps');
                Stack = squeeze(Fullvid(:,:,Plane,:));
                [Stack, StartX, StartY] = QP_package.cropXY(Stack);
                QPmap(:,:,:,Plane) = QP_package.getQP(Stack,s);
            end

            TotQPmap = sum(QPmap,4);
            Fig = figure();
            imagesc(mean(TotQPmap, 3))
            axis image;
            title('Click points, then press OK');
            [x, y] = ginput;
            coords = double([x y]);
            close(Fig) 

            figure()
            for i = 1:8%size(coords)
                Profile = squeeze(QPmap(round(coords(1,2)), round(coords(1,1)) , :, i));     
                plot(Profile)
                hold on
            end

            for Plane = 1:size(Fullvid, 3)
                Fig = figure();
                imagesc(mean(QPmap(:,:,:,Plane),3))
                axis image;
                colormap gray;
                title('Click points, then press OK');
                [x, y] = ginput;
                coords = double([x y]);
                close(Fig)  

                CoordsPlanes{Plane} = coords;
            end

            figure()
            for Plane = 1:size(Fullvid, 3)
                subplot(2, ceil(size(Fullvid, 3)./2), Plane)
                PlaneCoords = CoordsPlanes{1};
                for i = 1:size(PlaneCoords, 1)
                    c = round(PlaneCoords(i,:));
                    % ROI = squeeze(QPmap(c(2), c(1), :, Plane));
                    ROI = QPmap(c(2)-10:c(2)+10, c(1)-10:c(1)+10, :, Plane);
                    ROIs = squeeze(sum(sum(ROI,1),2));
                    plot(ROIs);
                    hold on
                    List(:, Plane) = ROIs;
                end
            end
            TotPhase = sum(List, 2);
                
        end

        function [Results] = calibrateAlpha(obj, q)
            QPmap = obj.QPmap(:, 25:400, :,:);
            MinQP = min(QPmap,[], 4);
            PlaneInFocus = 1;
            OutputFolder = append(obj.raw.movInfo.Path, filesep, 'Output_CalibrateAlpha');
            mkdir(append(obj.raw.movInfo.Path, filesep, 'Output_CalibrateAlpha'));
            zpos = [1:size(QPmap, 4)];

            for plane = 1:size(MinQP, 3)
                I = MinQP(:,:,plane); 
                h = 1;             
                I2 = imhmin(I, h);
                bw = imregionalmin(I2);
                bw = bw & (I < prctile(I(:), 0.25)); 
                PartCoord = regionprops(bw, 'Centroid');
                CoordsOverPlanes{plane} = PartCoord;
            end

            numPlanes  = numel(CoordsOverPlanes);
            distThresh = 10;
            centroids = cell(1, numPlanes);
            for p = 1:numPlanes
                centroids{p} = cat(1, CoordsOverPlanes{p}.Centroid);
            end
            refCoords = centroids{1};
            numRef    = size(refCoords, 1);
            CleanCoordsOverPlanes = [];
            for i = 1:numRef
                refPt = refCoords(i, :);
                CleanParticle{1,1} = refPt;
                Fail = 0;
                for p = 2:numPlanes
                    dists = sqrt(sum((centroids{p} - refPt).^2, 2));
                    if all(dists >= distThresh)
                        CleanParticle{1,p} = nan;
                        Fail = 1;
                    else
                        idx = find(dists < distThresh, 1, 'first');
                        CleanParticle{1,p} = centroids{p}(idx,:);
                    end
                end

                if Fail == 0
                    CleanCoordsOverPlanes = [CleanCoordsOverPlanes;CleanParticle];
                end
            end

            Fig1 = figure();
            for plane = 1:size(MinQP, 3)
                subplot(2,4,plane)
                imagesc(MinQP(:,:,plane));
                hold on
                for j = 1:size(CleanCoordsOverPlanes,1)
                    c = CleanCoordsOverPlanes{j,plane};   % [x y]
                    plot(c(1), c(2), 'r.', 'MarkerSize', 10);
                end
                title(append('Plane ', num2str(plane)))
                clim([-3.15 3.15])
                sgtitle('Particle Pixels')
            end
            filename = append(OutputFolder, filesep, 'ParticlePositions.png');
            saveas(Fig1, filename);

            % Fig2 = figure();
            % for plane = 1:numPlanes
            %     PartCoord = CleanCoordsOverPlanes(:,plane);
            %     subplot(2,4,plane)
            %     for i = 1:numel(PartCoord)
            %         c = round(PartCoord{i,1});
            %         List = squeeze(QPmap(c(2), c(1), plane, :));
            %         for Frame = 1:size(QPmap, 4)
            %             ROI = imgaussfilt(QPmap(c(2)-6:c(2)+6, c(1)-6:c(1)+6, plane, Frame), 3);
            %             List(Frame,1) = max(ROI, [], 'all');
            %         end
            %         plot(List);
            %         hold on
            %     end
            % end

            figure()
            for i = 1:size(CleanCoordsOverPlanes, 1)
                PartCoord = CleanCoordsOverPlanes(i,:);
                for plane = 1:numPlanes
                    c = round(PartCoord{:,plane});
                    ROI(:,:,plane,:) = QPmap(c(2)-10:c(2)+10, c(1)-10:c(1)+10, plane, :);
                end
                ROIs = squeeze(sum(sum(sum(ROI,1),2),3));
                plot(ROIs);
                hold on
            end
                    % CorrectAngle = 0;
                    % Angle(1) = List(1);
                    % for j = 2:size(List)
                    %     if List(j) - List(j-1) < -2
                    %         CorrectAngle = CorrectAngle + pi*2;
                    %     elseif List(j) - List(j-1) > 2
                    %         CorrectAngle = CorrectAngle - pi*2;
                    %     end
                    %     Angle(j) = List(j) + CorrectAngle;
                    % end
                    % y = Angle(:);                 % column vector
                    % x = zpos';            % x-axis
                    % plot(curve, zpos, Angle);
                    
                    % % Initial guesses
                    % A0 = max(y) - min(y);
                    % [~, idx] = max(y);
                    % mu0 = x(idx);
                    % sigma0 = numel(y)/10;
                    % C0 = min(y);
                    % 
                    % % Fit
                    % ft = fittype('A*exp(-(x-mu)^2/(2*sigma^2)) + C', ...
                    %              'independent','x','coefficients',{'A','mu','sigma','C'});
                    % 
                    % opts = fitoptions(ft);
                    % opts.StartPoint = [A0 mu0 sigma0 C0];
                    % opts.Lower = [0 0 0 -Inf];     % constrain amplitude & width
                    % 
                    % [curve, gof] = fit(x, y, ft, opts);
                    % 
                    % if gof.rsquare > 0.70
                    %     height = curve.A;
                    %     width  = curve.sigma;
                    %     center = curve.mu;
                    %     baseline = curve.C;
                    % 
                    %     plot(curve, zpos, Angle);
                    %     legend off
                    %     Results.height(i, plane) = height;
                    %     Results.width(i,plane) = width;
                    %     Results.center(i,plane) = center;
                    %     Results.baseline(i,plane) = baseline;
                    % 
                    %     ParticleNames{i} = append('Particle ', num2str(i));
                    %     hold on
                    % else
                    %     Results.height(i, plane) = nan;
                    %     Results.width(i,plane) = nan;
                    %     Results.center(i,plane) = nan;
                    %     Results.baseline(i,plane) = nan;
                    % end
            %     end
            %     xlabel('z-position in stack');
            %     ylabel('Phase (°)')
            % end
            % filename = append(OutputFolder, filesep, 'PhaseProfiles.png');
            % saveas(Fig2, filename);
            % save(append(OutputFolder, filesep, 'CalcAlphaFitResults.png'), "Results");
            Results = [];
        end
    end
end

