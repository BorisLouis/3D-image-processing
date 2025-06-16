classdef Channel2DCalibration < handle
    %CHANNEL2DCALIBRATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        ext
        Calibrations
        info
        allCal
        cal    
    end
    
    methods
        function obj = Channel2DCalibration(path2MPCal,info)
            obj.path = path2MPCal.path;
            obj.ext  = path2MPCal.ext;
            obj.info = info;
        end
        
        function [cal] = getCal(obj)
            cal = obj.cal;
        end

        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            folderContent = dir(path);
            %Get how many folder are in the main folder
            idx = sum(cellfun(@sum,{folderContent.isdir}));
            %Matlab always store ., .. as folder for relative path so we
            %want to find more than 2 folder in folderContent.
            % assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            
            obj.path = path;
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
                path = obj.raw.fullPath;
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(path);
                assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
                assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');
                if strcmp(obj.info.method, 'Darkfield Phase')
                    DarkfieldPhase = 1;
                else
                    DarkfieldPhase = 0;
                end
                [calib, inform, transformations] = mpSetup.cali.calculate(path,nChan,DarkfieldPhase);

                calibration.info = inform;
                calibration.file = calib;
                if calibration.file.multiModal == true
                    calibration.file.transformations = transformations;
                else
                end
                [folder,~] = fileparts(path);
                filename = [folder filesep 'calibration.mat'];
                calibration.fullPath = filename;
                save(filename,'calibration');
                
            end
            
            obj.cal = calibration;
            disp('Done');
            
        end

        function retrieveMovies(obj)
            disp('Retrieving movies from indicated folder...')
            %we get the MPCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                if ~isempty(file2Analyze)
                    
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    %we extract z motor position to check if the movie
                    %is indeed a zCalibration (expect zStack)
                    tmp = Core.ChannelCalibration(file, obj.info);
                    obj.Calibrations.(['MPCal' num2str(i-2)]) = tmp;

                else

                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);

                end

            end
            disp('=======> DONE ! <========')
        end

        function calcIndivCal(obj)
            
            fieldN = fieldnames(obj.Calibrations);
            nfields = numel(fieldN);
      
            for i = 1: nfields
                
                if isfield(obj.info,'nChan')
                    nChan = obj.info.nChan;
                else
                    nChan = 1;
                end
                
                disp(['Retrieving data from MPCalibration file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                
                obj.Calibrations.(fieldN{i}).calc(nChan);
                currentCal = obj.Calibrations.(fieldN{i}).getCal;
                
                obj.allCal(i).file = currentCal.file;
                
            end
        end
        
        function calcCombinedCal(obj)
            disp('Combining data from different calibration...');
            allData = obj.allCal;
            nFiles = length(allData);  

            for i = 1:nFiles
                allScale(i,:) = allData(i).file.Scale; 
                allRotational(i,:)  = allData(i).file.RotationAngle;   
                allTranslational(i, :, :) = allData(i).file.Translation;
            end

            Scale = mean(allScale, "all");
            Rotational = mean(allRotational, "all");
            Translational(1,1) = mean(allTranslational(:,1,:), 'all');
            Translational(1,2) = mean(allTranslational(:,2,:), 'all');

            tform = simtform2d(Scale, Rotational, Translational);

            fileName = [obj.path filesep '2DCal.mat'];
            save(fileName,'tform')
            
            disp('================>DONE<====================');           
        end
    end
end

