classdef TICSExperiment < handle
    %TICSEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        ext
        TICSMovies
        info
        cal2D
    end
    
    methods
        function obj = TICSExperiment(folder2Data,cal2D,info,infoChannel)
            
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.info = info;
            if ~isempty(infoChannel)
                obj.info = cell2struct([struct2cell(info); struct2cell(infoChannel)], ...
                         [fieldnames(info); fieldnames(infoChannel)], 1);
            end
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
                    if strcmp(obj.info.Dimension, '2D')
                        obj.cal2D = cal2D;
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
        end

        function retrieveTICSData(obj, q)
            fieldsN = fieldnames(obj.TICSMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            ViscName = cell([nfields, 1]);
            ViscValuesMean = nan([nfields, 1]);
            ViscValuesStd = nan([nfields, 1]);
            for i = 1: nfields
                try
                    disp(['Retrieving data from TICS file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.TICSMovies.(fieldsN{i});
                    currentTrackMov.LoadAllFrames;
            
                    currentTrackMov.calculateOmega;
                    currentTrackMov.getAutocorrmap;
                    currentTrackMov.getDiffusionmap;
      
                catch
                    disp(append('Failed TICS movie ', num2str(i), ' / ', num2str(nfields), ' ...'));
                end

                ViscName(i,1) = {currentTrackMov.raw.fullPath};
                ViscValuesMean(i, 1) = currentTrackMov.Results{1, 1}.ViscMean;
                ViscValuesStd(i, 1) = currentTrackMov.Results{1, 1}.ViscStd;

                ViscAnalysis = MicrorheologyAnalysis.analyzeViscosityMap(currentTrackMov.ViscosityMap{1, 1}, 'PixelSize', currentTrackMov.info.PxSize,...
                    'ClipNegative', 1, 'PatchPercentile', 75, 'ThresholdSweep', [50:5:95], 'NumGrayLevels', 16,...
                    'NVariogramPairs', 200000, 'PlotResults', 0, 'Title', currentTrackMov.raw.movInfo.Path);
                save(append(currentTrackMov.raw.movInfo.Path, filesep, 'ViscosityMapAnalysis.mat'), "ViscAnalysis");
            end
            
            ViscResults = table(ViscName, ViscValuesMean, ViscValuesStd, 'VariableNames', {'file', 'mean visc (cP)', 'std visc (cP)'});
            save(append(obj.path, filesep, 'TICSResults_viscosity.mat'), "ViscResults");
        end
    end
end

