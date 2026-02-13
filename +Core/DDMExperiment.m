classdef DDMExperiment < handle
    %DDMEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        ext
        DDMMovies
        info
        cal2D
    end
    
    methods
        function obj = DDMExperiment(folder2Data,cal2D,info,infoChannel)
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

        function retrieveDDMMovies(obj,q)
            fieldsN = fieldnames(obj.DDMMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];

            Cell = [];
            Diffusion = [];
            AnExp = [];
            Viscosity = [];
            IntMean = [];
            IntVar = [];
            Size = [];

            for i = 1: nfields
                try
                    disp(['Retrieving data from DDM file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.DDMMovies.(fieldsN{i});
    
                    if strcmp(currentTrackMov.info.ddmParam.Scanning, 'off')
                        currentTrackMov.getFullFrames;
                        currentTrackMov.mainDDM('NumBins',30);
                        currentTrackMov.fitDDMNoOptim;
                        % currentTrackMov.fitDDM;
                        currentTrackMov.getParams;
                        disp(['DDM analysis completed for file ' num2str(i) ' / ' num2str(nfields) ' ...']);   
                    elseif strcmp(currentTrackMov.info.ddmParam.Scanning, 'on')
                        [PxIdx] = currentTrackMov.getROI;
                        for z = 1:8
                            ViscosityMap{z,1} = nan(currentTrackMov.raw.movInfo.Length, currentTrackMov.raw.movInfo.Width);
                            DiffusionMap{z,1} = nan(currentTrackMov.raw.movInfo.Length, currentTrackMov.raw.movInfo.Width);
                            AnExpMap{z,1} = nan(currentTrackMov.raw.movInfo.Length, currentTrackMov.raw.movInfo.Width);
                        end
                        for CurrentPx = PxIdx'
                            currentTrackMov.getFrameParts(CurrentPx);
                            currentTrackMov.mainDDM('NumBins',30);
                            currentTrackMov.fitDDM;
                            currentTrackMov.getParams;
                            for j = 1:size(currentTrackMov.MSDResults, 1)
                                ViscosityMap{j,1}(CurrentPx) = currentTrackMov.MSDResults{j,1}.n;
                                AnExpMap{j,1}(CurrentPx) = currentTrackMov.MSDResults{j,1}.alpha;
                                DiffusionMap{j,1}(CurrentPx) = currentTrackMov.MSDResults{j,1}.Diff;
                            end
                        end
                        currentTrackMov.MSDResults{1, 1}.Diff = nanmean(DiffusionMap{1,1});
                        currentTrackMov.MSDResults{1, 1}.alpha = nanmean(AnExpMap{1,1});
                        currentTrackMov.MSDResults{1, 1}.n = nanmean(ViscosityMap{1,1});

                        FileName = append(currentTrackMov.calibrated{1, 1}.mainPath, filesep, 'ViscosityMap.mat');
                        save(FileName, 'ViscosityMap');
                        FileName = append(currentTrackMov.calibrated{1, 1}.mainPath, filesep, 'AnExpMap.mat');
                        save(FileName, 'AnExpMap');
                        FileName = append(currentTrackMov.calibrated{1, 1}.mainPath, filesep, 'DiffusionMap.mat');
                        save(FileName, 'DiffusionMap');
                        disp(['DDM analysis completed for file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    end
                    Cell{end+1, 1} = currentTrackMov.raw.movInfo.Path;
                    Diffusion = [Diffusion; currentTrackMov.MSDResults{1, 1}.Diff];
                    AnExp = [AnExp; currentTrackMov.MSDResults{1, 1}.alpha];
                    Viscosity = [Viscosity; currentTrackMov.MSDResults{1, 1}.n];
                    IntMean = [IntMean; currentTrackMov.IntResults.MeanInt];
                    IntVar = [IntVar; currentTrackMov.IntResults.VarInt];
                    Size = [Size; currentTrackMov.SizeResults];
                catch
                    disp(append('Failed DDM movie ', num2str(i), ' / ', num2str(nfields), ' ...'));
                end
            end

            DDMResults = table(Cell, Diffusion, AnExp, Viscosity, IntMean, IntVar, Size, 'VariableNames', {'Sample', 'Diff (µm^2/s)', 'Anom. exp',...
                'Viscosity (cP)', 'Intensity mean', 'Intensity var', 'Size (µm^2)'});
            Filename = append(obj.path, filesep, 'DDMResults.mat');
            save(Filename, 'DDMResults')
        end
    end
end

