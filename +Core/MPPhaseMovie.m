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


            for i = 1:nFrames
                waitbar(i./nFrames,f,append('Calculating phase map ', num2str(i),' out of ', num2str(nFrames)));
                Stack = obj.getFrame(i, q);
                [Stack, StartX, StartY] = QP_package.cropXY(Stack);
                QPmap(:,:,:,i) = QP_package.getQP(Stack, s);
            end

            obj.QPmap = QPmap;
            obj.Cropped.StartX = StartX;
            obj.Cropped.StartY = StartY;
        end

    end
end

