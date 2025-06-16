classdef ChannelCalibration < Core.Movie
    %CHANNELCALIBRATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cal
    end
    
    methods
        function obj = ChannelCalibration(raw,info)
            %In this case raw path is used for both the loading of the
            %movie and the calibration calculation
            obj = obj@Core.Movie(raw,info);
          
        end
        
        function set.cal(obj,cal)
            assert(isstruct(cal),'2DCal expected to be a struct');
            obj.cal = cal;
                
        end
        
        function calc(obj, nChan)
            
            switch nargin
                case 1
                    nChan = 1;
                case 2
                otherwise
                    error('too many input arguments')
            end
            
            assert(~isempty(obj.raw),'Please provide a path');
            %Calculate from the raw path stored in the movie
            path = obj.raw.movInfo.Path;
            [file2Analyze] = obj.getFileInPath( path, 'calibration.mat');
            
            if (~isempty(file2Analyze))
                %If calibration already exist we load it
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                calibration = tmp.calibration;

            else
                %Otherwise we Calculate it
                f = waitbar(0, 'initializing');

                if and(contains(obj.raw.fullPath1, 'Cam1'), contains(obj.raw.fullPath2, 'Cam2'))
                    path1 = obj.raw.fullPath1;
                    path2 = obj.raw.fullPath2;
                else and(contains(obj.raw.fullPath1, 'Cam2'), contains(obj.raw.fullPath2, 'Cam1'))
                    path1 = obj.raw.fullPath2;
                    path2 = obj.raw.fullPath1;
                end
                [frameInfo1, movInfo1] = Load.Movie.(erase(obj.raw.ext, '.')).getInfo(path1);
                [frameInfo2, movInfo2] = Load.Movie.(erase(obj.raw.ext, '.')).getInfo(path2);

                [optimizer,metric] = imregconfig("multimodal");
                
                for i = 1:movInfo1.maxFrame
                    waitbar(i./movInfo1.maxFrame,f,'Loading frames');
                    [mov1(:,:,i)] = Load.Movie.(erase(obj.raw.ext, '.')).getFrame(path1, i);
                    [mov2(:,:,i)] = Load.Movie.(erase(obj.raw.ext, '.')).getFrame(path2, i);
                end

                for i = 1:movInfo1.maxFrame
                    waitbar(i./movInfo1.maxFrame,f,'Calculating transformation');
                    tform = imregcorr(mov1(:,:,i), mov2(:,:,i), 'similarity');
                    Scale(i) =  tform.Scale;
                    RotationAngle(i) = tform.RotationAngle;
                    Translation(:,:,i) = tform.Translation;
                end

                waitbar(0.99,f,'Saving in file');

                calibration.file.Scale = Scale;
                calibration.file.RotationAngle = RotationAngle;
                calibration.file.Translation = Translation;

                filename = [path filesep 'calibration.mat'];
                calibration.fullPath = filename;
                save(filename,'calibration');
                
                close(f)
            end
            
            obj.cal = calibration;
            disp('Done');
            
        end
        
        function cal = getCal(obj)
            
            cal = obj.cal;
            
        end
    end
end

