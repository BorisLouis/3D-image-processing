classdef MPSegmentMovie < Core.MPMovie

    properties
    end
    
    methods
        function obj = MPSegmentMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function getSegmentMovie(obj, q)
            if strcmp(obj.info.frame2Load, 'all')
                nFrames = obj.calibrated{1, 1}.nFrames; 
            elseif isa(obj.info.frame2Load, 'double')
                nFrames = max(obj.info.frame2Load);
            end

            for i = 1:nFrames
                Stack = obj.getFrame(i, q);
                QP = QP_package.getQP(Stack, s);
            end
        end
    end
end

