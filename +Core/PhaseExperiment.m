classdef PhaseExperiment < handle 
    %PHASEEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        ext
        PhaseMovies
        info
        cal2D
    end
    
    methods
        function obj = PhaseExperiment(folder2Data,cal2D,info,infoChannel)
            
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.info = info;
            if ~isempty(infoChannel)
                obj.info = cell2struct([struct2cell(info); struct2cell(infoChannel)], ...
                         [fieldnames(info); fieldnames(infoChannel)], 1);
            end
            obj.info.optics.dx = obj.info.PxSize*10^(-3);
            obj.cal2D = cal2D;
        end
        
        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            folderContent = dir(path);
            %Get how many folder are in the main folder
            idx = sum(cellfun(@sum,{folderContent.isdir}));
            %Matlab always store ., .. as folder for relative path so we
            %want to find more than 2 folder in folderContent.
            assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            
            obj.path = path;
            
        end

        function set.cal2D(obj,cal2D)
            if isempty(cal2D)
                obj.cal2D = [];
            elseif isstruct(cal2D)
                obj.cal2D = cal2D;
            else
                assert(ischar(cal2D), 'Path should be given as a string');
                assert(isfolder(cal2D), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')

                [file2Analyze] = Core.Movie.getFileInPath(cal2D,'2DCal.mat');

                if strcmp(obj.info.Dimension, '2D')
                    obj.info.runCal = false;
                end

                if isempty(file2Analyze)
                    if strcmp(obj.info.Dimension, '2D')
                        obj.info.runCal = true;
                        obj.cal2D = cal2D;
                    else
                        error('No 2D calibration file found in the given folder');
                    end
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(and(isstruct(cal), and(isfield(cal,'camConfig'),isfield(cal,'file'))),...
                        '2D calibration is supposed to be a struct with 4 fields');
                    obj.cal2D = cal;

                end
            end
        end

        function retrievePhaseMask(obj, q)
            fieldsN = fieldnames(obj.PhaseMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                try
                    disp(['Retrieving data from phase file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.PhaseMovies.(fieldsN{i});
    
                    currentTrackMov.getPhaseMovie(q);
                    [MinPhasePart(:,i)] = currentTrackMov.calibrateAlpha(q);
                catch
                    disp(append('Failed phase masking - Movie ', num2str(i), ' / ', num2str(nfields), ' ...'));
                end
            end
            % save(append(obj.path, filesep, "CalcAlpha_alpha1_70.mat"), "MinPhasePart");
            % disp(append('saved min projection values to calc alpha - file ', obj.path))
        end
    end
end

