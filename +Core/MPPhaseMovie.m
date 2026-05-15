classdef MPPhaseMovie < Core.MPMovie
    %MPPHASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        QPmap
        QPmapPerPlane
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
                        [QPmap(:,:,:,n), ~] = QP_package.getQP(Stack,s, obj.calibrated{1,1}.oRelZPos);
                        n = n+1;
                    end
                end
                close(f)
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat');
                save(Filename, 'QPmap', '-v7.3');
                obj.QPmap = QPmap;
               
                obj.Cropped.StartX = StartX;
                obj.Cropped.StartY = StartY;
            else
                disp('Found Phasemap - loading it')
                load(append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie.mat'))
                Stack = obj.getFrame(1, q);
                obj.Cropped.StartX = floor((size(Stack,1) - size(QPmap, 1))./2);
                obj.Cropped.StartY = floor((size(Stack,2) - size(QPmap, 2))./2);
                obj.QPmap = QPmap;
                disp('Found Phasemap - loaded')
            end
        end

        function getPhaseMapPerPlane(obj, q)
            if strcmp(obj.info.runMethod, 'load')
                if exist(append(obj.raw.movInfo.Path, filesep, 'PhaseMapsPerPlane'))
                    run = 0;
                else
                    run = 1;
                end
            else
                run = 1;
            end

            run = 1;
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
    
                mkdir(append(obj.raw.movInfo.Path, filesep, 'PhaseMapsPerPlane'));
                
                for k = 1:nFrames 
                    waitbar(k./nFrames,f,append('Loading frame ', num2str(k)));
                    Stack(:,:,:,k) = obj.getFrame(k, q);
                    MotorPos(k,:) = obj.raw.frameInfo(2*k).Pos;
                end

                RelzPos = MotorPos(:,3) - mean(MotorPos(:,3));
                for k = 1:size(Stack, 3)
                    waitbar(k./size(Stack, 3),f,append('Calculating phase map - plane ', num2str(k)));
                    PlaneStack = squeeze(Stack(:,:,k,:));
                    [PlaneStack, StartX, StartY] = QP_package.cropXY(PlaneStack);
                    [QPmapPerPlane(:,:,:,k), ~] = QP_package.getQP(PlaneStack,s, RelzPos);             
                end
                close(f)
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMapsPerPlane', filesep, 'QPmapPerPlane.mat');
                save(Filename, 'QPmapPerPlane');
                obj.QPmapPerPlane = QPmapPerPlane;
               
                obj.Cropped.StartX = StartX;
                obj.Cropped.StartY = StartY;
            else
                disp('Found Phasemap per plane- loading it')
                load(append(obj.raw.movInfo.Path, filesep, 'PhaseMapsPerPlane', filesep, 'QPmapPerPlane.mat'))
                Stack = obj.getFrame(1, q);
                obj.Cropped.StartX = floor((size(Stack,1) - size(QPmapPerPlane, 1))./2);
                obj.Cropped.StartY = floor((size(Stack,2) - size(QPmapPerPlane, 2))./2);
                obj.QPmapPerPlane = QPmapPerPlane;
                disp('Found Phasemap per plane - loaded')
            end
        end
    end
end

