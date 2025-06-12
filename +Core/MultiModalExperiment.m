classdef MultiModalExperiment < handle
    %MULTIMODALEXPERIMENT 
    %Class made to combine multiple channels. Take the info, calibrate both
    %channels, create right objects for every channel and later on combine
    %info and process data when needed
    
    properties
        path
        ext
        MoviesCh1
        MoviesCh2
        cal2D
        info
        info1
        info2
        ZCal
        SRCal
        SRCal2
        traces3D
        traces3Dcommon
        MSD
    end
    
    methods
        function obj = MultiModalExperiment(folder2Data,cal2D,info, info1, info2,SRCalPath,zCalPath)
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.cal2D = cal2D;
            obj.info = info;
            obj.info1 = info1;
            obj.info2 = info2;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath; 

            if strcmp(obj.info.Channel1, 'Translational Tracking')
            elseif strcmp(obj.info.Channel1, 'Rotational Tracking')
                obj.MoviesCh1 = Core.TrackingExperimentRotational(folder2Data, obj.cal2D, obj.info, obj.SRCal.path, obj.ZCal.path, 1, obj.info1);
            elseif strcmp(obj.info.Channel1, 'Phase')
            elseif strcmp(obj.info.Channel1, 'Segmentation')
            end

            if strcmp(obj.info.Channel2, 'Translational Tracking')
            elseif strcmp(obj.info.Channel2, 'Rotational Tracking')
                obj.MoviesCh2 = Core.TrackingExperimentRotational(folder2Data, obj.cal2D, obj.info, obj.SRCal.path, obj.ZCal.path, 2, obj.info2);
            elseif strcmp(obj.info.Channel2, 'Phase')
            elseif strcmp(obj.info.Channel2, 'Segmentation')
            end

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
            else
                assert(ischar(cal2D), 'Path should be given as a string');
                assert(isfolder(cal2D), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')

                [file2Analyze] = Core.Movie.getFileInPath(cal2D,'2DCal.mat');

                if isempty(file2Analyze)
                    error('No 2D calibration file found in the given folder');
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

        function set.SRCal(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal.cal = SRCal;
                obj.SRCal.path = SRCal;
            else
                assert(isfolder(SRCal), 'The given path is not a folder');
               
                for q = 1:2
                %Check Given path
                    [file2Analyze] = Core.Movie.getFileInPath(SRCal,append('SRCalibration', num2str(q), '.mat'));
    
                    if isempty(file2Analyze)
                        error('No SR calibration file found in the given folder');
                    else
                        fileName = [file2Analyze.folder filesep file2Analyze.name];
                        cal = load(fileName);
                        field = fieldnames(cal);
                        cal = cal.(field{1});
                        assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                            'SR calibration is supposed to be a struct with 2 fields');
    
                        if q == 1
                            obj.SRCal.cal = cal;
                            obj.SRCal.path = SRCal;
                        elseif q == 2
                            obj.SRCal2.cal = cal;
                            obj.SRCal2.path = SRCal;
                        end
                    end
                end
            end
            
        end

         function set.ZCal(obj,zCal)
            if isempty(zCal)
                obj.ZCal.cal = zCal;
                obj.ZCal.path = zCal;
            else
                assert(isfolder(zCal), 'The given path is not a folder');

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(zCal,'zCalibration.mat');

                if isempty(file2Analyze)
                    error('No z calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(isstruct(cal),'zCalibration is supposed to be in cells format');
                    assert(and(isfield(cal,'fitZParam'),isfield(cal,'calib')),...
                        'Something is wrong in the fields of your Z calibration');

                    obj.ZCal.cal = cal;
                    obj.ZCal.path = zCal;
                end
            end    
         end


         function RetrieveMovies(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory

            %Load the movies for planes 1-8
            for i = 3:size(folder2Mov,1)
                %Check if the directory
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                if ~isempty(file2Analyze)
                    
                    obj.info.multiModal = 1;
                    count = 0;
                    count = count+1;
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;

                    tmp = Core.MPMovie(file , obj.cal2D, obj.info);
                    
                    if count == 1
                        tmp.giveInfo;
                    else
                        tmp.info = obj.trackMovies.(['mov' num2str(1)]).getInfo; 
                    end
                    tmp.calibrate;

                    if strcmp(obj.info.Channel1, 'Rotational Tracking')
                        Movie1 = Core.MPTrackingMovieRotational(file , obj.MoviesCh1.cal2D, obj.MoviesCh1.info, obj.MoviesCh1.SRCal.path, obj.MoviesCh1.ZCal.path, 1);
                        Movie1.calibrated = tmp.calibrated{1,1};    
                        if count == 1
                            Movie1.giveInfo;
                        else
                            Movie1.info = obj.MoviesCh1.trackMovies.(['mov' num2str(1)]).getInfo; 
                        end
                        obj.MoviesCh1.trackMovies.(['mov' num2str(((i-2)*2)-1)]) = Movie1;
                    elseif strcmp(obj.info.Channel1, 'Translational Tracking')
                    elseif strcmp(obj.info.Channel1, 'Phase')
                    elseif strcmp(obj.info.Channel1, 'Segmentation')
                    end

                    if strcmp(obj.info.Channel2, 'Rotational Tracking')
                        Movie2 = Core.MPTrackingMovieRotational(file , obj.MoviesCh2.cal2D, obj.MoviesCh2.info, obj.MoviesCh2.SRCal.path, obj.MoviesCh2.ZCal.path, 1);
                        Movie2.calibrated = tmp.calibrated{1,2};
                        if count == 1
                            Movie2.giveInfo;
                        else
                            Movie2.info = obj.MoviesCh1.trackMovies.(['mov' num2str(1)]).getInfo; 
                        end
                        obj.MoviesCh2.trackMovies.(['mov' num2str(((i-2)*2)-1)]) = Movie2;
                    elseif strcmp(obj.info.Channel2, 'Translational Tracking')
                    elseif strcmp(obj.info.Channel2, 'Phase')
                    elseif strcmp(obj.info.Channel2, 'Segmentation')
                    end

                else
                    
                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ' obj.ext ' file and is therefore ignored']);
                    
                end
                
            end       
            
            if isempty(obj.MoviesCh1)              
               error(['No %s was found for planes 9-16. Please check:\n',...
                   '1) that you gave the correct file extension.\n',...
                   '2) that you gave the path of a folder containing folders containing movies with the given extension'],obj.ext);        

            end
            disp('=======> DONE ! <========')

          end

          function RunAnalysis(obj)
              %%% first run channel 1 analysis
              if strcmp(obj.info.Channel1, 'Segmentation')
              elseif strcmp(obj.info.Channel1, 'Phase')
              elseif strcmp(obj.info.Channel1, 'Translational Tracking')
              elseif strcmp(obj.info.Channel1, 'Rotational Tracking')
                    frame = obj.info1.testFrame;
                    % testMov = obj.MoviesCh1.trackMovies.mov1;
                    % testMov2 = obj.MoviesCh2.trackMovies.mov1;
                    % testMov.getTransformation(testMov2, frame);
                    % testMov.getROIs;
                    % testMov2.getROIs;
                    % testMov.showCandidate(testMov2, frame);

                    val2Use = 'bestFocus';
                    obj.MoviesCh1.retrieveTrackData(obj.MoviesCh1.info.detectParam,obj.MoviesCh1.info.trackParam, 1);
                    traces = obj.getTraces3D;
                    
                    trackingExp.ConsolidateChannels3;
                    
                    %% save Data
                    trackingExp.saveData;
              end
          end

          function retrieveTrackData(obj,detectParam, trackParam)
                %Checking user input
                assert(nargin==3, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter');
                %assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
                assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                    'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
                fieldsN = fieldnames(obj.trackMovies);
                %Extraction of Data
                nfields = numel(fieldsN);
                allTraces = [];
                for i = 1: nfields
                    
                    disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.trackMovies.(fieldsN{i});
                    
                    %Molecule detection
                    % if currentTrackMov.info.rotationalCalib == 0
                        currentTrackMov.findCandidatePos(detectParam);
                    % else
                    % end
                    
                    %SR fitting
                    currentTrackMov.SRLocalizeCandidate(detectParam);
                    refPlane = round(currentTrackMov.calibrated{1,1}.nPlanes/2);
                    rot = true;
                    %apply SRCal
                    currentTrackMov.applySRCal(rot,refPlane);
                    
                    %apply ZCal
                    currentTrackMov.applyZCal;
                    
                    %Plane consolidation
                    frames = 1:currentTrackMov.calibrated{1,1}.nFrames;
                    currentTrackMov.consolidatePlanes(frames,detectParam)
                    
                    %superResolve
                    currentTrackMov.superResolve;
                    
                    %tracking occurs here
                    currentTrackMov.trackParticle(trackParam);
                    
                    [traces] = currentTrackMov.getTraces;
                    for q = 1:length(traces)
                        allTraces = [];
                        fileN = cell(length(traces{q,1}),1);
                        fileN(:,1) = {i};
                   
                        [xStep,xMotor] = currentTrackMov.getXPosMotor;
                        [yStep,yMotor] = currentTrackMov.getYPosMotor;
                        [zSt,zMotor]   = currentTrackMov.getZPosMotor;
        
                        colMot = cell(length(traces{q,1}),1);
                        colMot(:,1) = {xMotor};
                        colStep = cell(length(traces{q,1}),1);
                        colStep(:,1) = {xStep};
        
                        rowMot = cell(length(traces{q,1}),1);
                        rowMot(:,1) = {yMotor};
                        rowStep = cell(length(traces{q,1}),1);
                        rowStep(:,1) = {yStep};
        
                        zMot = cell(length(traces{q,1}),1);
                        zMot(:,1) = {zMotor};
                        zStep = cell(length(traces{q,1}),1);
                        zStep(:,1) = {zSt};
        
                        allTraces = [allTraces; traces{q,1}(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot ];
                        obj.traces3D{q,1} = allTraces;
                        obj.traces3D{q,2} = q; 
                    end
                end
                
                
    %             filename = [obj.path filesep 'traces3D.mat'];
    %             save(filename,'allTraces');
                
                
                disp('=================> DONE <===================');
            end
    end
end

