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
                    QPMapFrame = QP_package.getQP(Stack,s);
                    QPmap(:,:,:,i-Startidx) = QP_package.getQP(Stack, s);
                    n = n+1;
                end
                Filename = append(obj.raw.movInfo.Path, filesep, 'PhaseMovie', filesep, 'PhaseMovie', num2str(Step), '.mat');
                save(Filename, 'QPmap');
                QPmap = [];
            end
            close(f)

            obj.Cropped.StartX = StartX;
            obj.Cropped.StartY = StartY;
        end

    end
end

