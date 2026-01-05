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
        end

        function [MeanPhasePart] = calibrateAlpha(obj, q)
            QPmap = obj.QPmap(:, 25:400, :,:);
            MinQP = min(QPmap,[], 4);
            PlaneInFocus = 1;
            OutputFolder = append(obj.raw.movInfo.Path, filesep, 'Output_CalibrateAlpha');
            mkdir(append(obj.raw.movInfo.Path, filesep, 'Output_CalibrateAlpha'));
            zpos = [1:size(QPmap, 4)];
            Fig1 = figure()
            for plane = 1:size(MinQP, 3)
                I = MinQP(:,:,plane); 
                h = 1;             
                I2 = imhmin(I, h);
                bw = imregionalmin(I2);
                bw = bw & (I < prctile(I(:), 0.1)); 
                D = bwdist(bw);
                PartCoord = regionprops(bw, 'Centroid');
                subplot(2,4,plane)
                imagesc(MinQP(:,:,plane));
                hold on
                for i = 1:numel(PartCoord)
                    c = PartCoord(i).Centroid;   % [x y]
                    plot(c(1), c(2), 'r.', 'MarkerSize', 10);
                end
                title(append('Plane ', num2str(plane)))
                clim([-3.15 3.15])
                sgtitle('Particle Pixels')
                CoordsOverPlanes{plane} = PartCoord;
            end
            filename = append(OutputFolder, filesep, 'ParticlePositions.png');
            saveas(Fig1, filename);

            Fig2 = figure();
            for plane = 1:size(MinQP, 3)
                PartCoord = CoordsOverPlanes{plane};
                subplot(2,4,plane)
                for i = 1:numel(PartCoord)
                    c = round(PartCoord(i).Centroid);
                    plot(zpos, squeeze(QPmap(c(2), c(1), plane, :)));
                    Values(i, plane) = min(squeeze(QPmap(c(2), c(1), plane, :)));
                    ParticleNames{i} = append('Particle ', num2str(i));
                    hold on
                end
                title(append('Plane ', num2str(plane)))
                ylim([-3.15 3.15])
                xlabel('z-position in stack');
                ylabel('Phase (°)')
            end
            filename = append(OutputFolder, filesep, 'PhaseProfiles.png');
            saveas(Fig2, filename);

            Values(Values == 0) = NaN;
            save(append(OutputFolder, filesep, 'MinPhases.png'), Values);
        end
    end
end

