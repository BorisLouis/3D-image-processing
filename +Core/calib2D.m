classdef calib2D < Core.Movie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        cal
    end
    
    methods
        function obj = calib2D(raw,cal)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Core.Movie(raw);
            
            switch nargin
                case 1
                case 2
                    obj.cal = cal;
            end 
        end
        
        function set.cal(obj,cal)
            
            assert(isstruct(cal), 'Calibration is expected to be a struct');
            assert(numel(fieldnames(cal))==3, 'Calibration is expected to have 3 Fields');
            assert(isfield(cal,'file'),'One of the field should be "file" ');
            
            obj.cal = cal;
                
        end
        
        function calc(obj)
            
            path = obj.raw.movInfo.Path;
            [file2Analyze] = obj.getFileInPath( path, '.mat');
            
            if (~isempty(file2Analyze))
                
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                calibration = tmp.calibration;

            else
            
                path = obj.raw.fullPath;
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(path);
                assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
                assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');

                [calib, inform] = mpSetup.cali.calculate(fullPath);

                calibration.info = inform;
                calibration.file = calib;

                filename = [file2Analyze.folder filesep 'calibration.mat'];
                calibration.fullPath = filename;
                save(filename,'calibration');
                
            end
            
            obj.cal = calibration;
            disp('Done');
            
        end
        
        function cal = getCal(obj)
            cal = obj.cal;
        end
        
    end
end

