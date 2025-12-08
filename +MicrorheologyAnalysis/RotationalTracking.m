classdef RotationalTracking < handle
    
    properties
        raw
        info
    end
    
    methods
        function obj = RotationalTracking(Path, info)
            obj.raw.Path = Path;
            obj.info = info;
        end
    end
end

