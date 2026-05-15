classdef TrackingExperimentRotational < handle

    properties
    
    path
    ext
    trackMovies
    cal2D
    info
    ZCal
    SRCal
    SRCal2
    traces3D
    traces3Dcommon
    MSD
        
    end

     methods
         function obj = TrackingExperimentRotational(folder2Data,cal2D,info,SRCalPath,zCalPath, q, infoChannel)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.info = info;
            if ~isempty(infoChannel)
                obj.info = cell2struct([struct2cell(info); struct2cell(infoChannel)], ...
                         [fieldnames(info); fieldnames(infoChannel)], 1);
            end
            obj.cal2D = cal2D;
            obj.ZCal = zCalPath;
            SRCalPath = {SRCalPath, q};
            obj.SRCal = SRCalPath; 
            % obj.SRCal2 = SRCalPath;
        end

        %Set function
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

        function set.SRCal(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal.cal = SRCal;
                obj.SRCal.path = SRCal;
            elseif isempty(SRCal{1,1})
                obj.SRCal.cal = SRCal;
                obj.SRCal.path = SRCal;
            else
                assert(isfolder(SRCal{1}), 'The given path is not a folder');
               

                q = SRCal{2};
                    %Check Given path
                    [file2Analyze] = Core.Movie.getFileInPath(SRCal{1},append('SRCalibration', num2str(q), '.mat'));
    
                    if isempty(file2Analyze)
                        error('No SR calibration file found in the given folder');
                    else
                        fileName = [file2Analyze.folder filesep file2Analyze.name];
                        cal = load(fileName);
                        field = fieldnames(cal);
                        cal = cal.(field{1});
                        assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                            'SR calibration is supposed to be a struct with 2 fields');
    
                        obj.SRCal.cal = cal;
                        obj.SRCal.path = SRCal{1};
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

        %get 3D traces
        
        function [traces3D] = getTraces3D(obj)
            traces3D = obj.traces3D;
        end


        function retrieveMovies(obj)
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
                    
                    count = 0;
                    count = count+1;
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    tmp = Core.MPTrackingMovieRotational(file , obj.cal2D, obj.info, obj.SRCal{1,1}.path, obj.ZCal.path);
                    
                    if count == 1
                        tmp.giveInfo;
                    else
                        tmp.info = obj.trackMovies.(['mov' num2str(1)]).getInfo; 
                    end
                    tmp.calibrate;
                    obj.trackMovies.(['mov' num2str(((i-2)*2)-1)]) = tmp;

                else
                    
                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ' obj.ext ' file and is therefore ignored']);
                    
                end
                
            end                
            
            if isempty(obj.trackMovies)              
               error(['No %s was found for planes 9-16. Please check:\n',...
                   '1) that you gave the correct file extension.\n',...
                   '2) that you gave the path of a folder containing folders containing movies with the given extension'],obj.ext);        
                
            end
            disp('=======> DONE ! <========')
        end

        function retrieveTrackDataPart1(obj,detectParam, trackParam, q)
            %Checking user input
            assert(nargin==4, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter and q');
            %assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
            fieldsN = fieldnames(obj.trackMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                try
                
                    disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.trackMovies.(fieldsN{i});
                    
                    %Molecule detection
                    % if currentTrackMov.info.rotationalCalib == 0
                        currentTrackMov.findCandidatePos(detectParam, q);
                    % else
                    % end
                    
                    %SR fitting
                    currentTrackMov.SRLocalizeCandidate(detectParam, q);
                    refPlane = round(currentTrackMov.calibrated{1,1}.nPlanes/2);
                    rot = true;
                    %apply SRCal
                    currentTrackMov.applySRCal(rot,refPlane);
                    
                    %apply ZCal
                    currentTrackMov.applyZCal;
                    
                    %Plane consolidation
                    if strcmp(obj.info.frame2Load, 'all')
                        frames = 1:currentTrackMov.calibrated{1,1}.nFrames;
                    elseif isa(obj.info.frame2Load, 'double')
                        frames = obj.info.frame2Load;
                    end
                    currentTrackMov.consolidatePlanes(frames,detectParam,q)
    
                    obj.trackMovies.(fieldsN{i}) = currentTrackMov;
                catch
                end
            end
        end

        function retrieveTrackDataPart2(obj,trackParam, q)
            %Checking user input
            assert(nargin==3, 'retrieveZCalData expects 2 inputs, 1)tracking parameter and q');
            %assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
            fieldsN = fieldnames(obj.trackMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                try

                    disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.trackMovies.(fieldsN{i});
                    
                    %superResolve
                    currentTrackMov.superResolve(q);
                    
                    currentTrackMov.correctDrift;
                    %tracking occurs here
                    currentTrackMov.trackParticle(trackParam,q);
                    
                    [traces] = currentTrackMov.getTraces;
    
                    allTraces = [];
                    fileN = cell(length(traces),1);
                    fileN(:,1) = {i};
               
                    [xStep,xMotor] = currentTrackMov.getXPosMotor;
                    [yStep,yMotor] = currentTrackMov.getYPosMotor;
                    [zSt,zMotor]   = currentTrackMov.getZPosMotor;
    
                    colMot = cell(length(traces),1);
                    colMot(:,1) = {xMotor};
                    colStep = cell(length(traces),1);
                    colStep(:,1) = {xStep};
    
                    rowMot = cell(length(traces),1);
                    rowMot(:,1) = {yMotor};
                    rowStep = cell(length(traces),1);
                    rowStep(:,1) = {yStep};
    
                    zMot = cell(length(traces),1);
                    zMot(:,1) = {zMotor};
                    zStep = cell(length(traces),1);
                    zStep(:,1) = {zSt};
    
                    allTraces = [allTraces; traces(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot ];
                    currentTrackMov.traces3D = [traces(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot];
                    obj.traces3D = allTraces;
                catch
                end
            end
            
            try
                filename = [obj.path filesep 'traces3D-', num2str(q), '.mat'];
                save(filename,'allTraces');
            catch
            end
            
            
            disp('=================> DONE <===================');
        end

        function retrieveTrackData(obj,detectParam, trackParam, q)
            %Checking user input
            assert(nargin==4, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter and q');
            %assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
            fieldsN = fieldnames(obj.trackMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                
                try
                    
                    disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                    currentTrackMov = obj.trackMovies.(fieldsN{i});

                    if strcmp(obj.info.runMethod, 'load')
                        filename1 = append(currentTrackMov.raw.movInfo.Path, filesep, 'calibrated', num2str(q), filesep,  'candidatePos.mat');
                        if exist(filename1)
                            run = 1;
                        else
                            run = 0;
                        end
                        filename2 = append(currentTrackMov.raw.movInfo.Path, filesep, 'calibrated', num2str(q), filesep,  'SRLocPos.mat');
                        if exist(filename2)
                            run  = 2;
                        end
                        filename3 = append(currentTrackMov.raw.movInfo.Path, filesep, 'calibrated', num2str(q), filesep,  'particle.mat');
                        if exist(filename3)
                            run  = 3;
                            load(filename3);
                            if isfield(particle, "SRList")
                                run = 4;
                            end
                        end
                        filename5 = append(currentTrackMov.raw.movInfo.Path, filesep, 'calibrated', num2str(q), filesep,  'Traces3D.mat');
                        if exist(filename5)
                            run  = 5;
                        end
                    else
                        run = 0;
                    end

                    if run < 1

                        if isempty(currentTrackMov.calibrated{1,1})
                            load(append(currentTrackMov.raw.movInfo.Path, filesep, 'calibrated', num2str(q), filesep,...
                                'calibrated', num2str(q), '.mat'));
                            currentTrackMov.calibrated{1,1} = calib;
                        end
                        currentTrackMov.findCandidatePos(detectParam, q);
                        currentTrackMov.showCandidateSingleChan(obj.info.TestFrame, q);
                    else
                        load(filename1)
                        currentTrackMov.candidatePos = candidate;
                        disp('Found candidate file - loaded it')
                    end

                    if strcmp(detectParam.fitting, 'on')
                        if run < 2
                            currentTrackMov.SRLocalizeCandidate(detectParam, q);
                            refPlane = round(currentTrackMov.calibrated{1,1}.nPlanes/2);
                            rot = true;
                            currentTrackMov.applySRCal(rot,refPlane);
                            currentTrackMov.applyZCal;
                        else
                            load(filename2)
                            currentTrackMov.corrLocPos = locPos;
                            currentTrackMov.unCorrLocPos = locPos;
                            disp('Found localisation positions file - loaded it')
                        end
    
                        if run < 3 
                            if strcmp(obj.info.frame2Load, 'all')
                                frames = 1:currentTrackMov.calibrated{1,1}.nFrames;
                            elseif isa(obj.info.frame2Load, 'double')
                                frames = obj.info.frame2Load;
                            end
                            currentTrackMov.consolidatePlanes(frames,detectParam,q)
                        else 
                            load(filename3)
                            currentTrackMov.particles = particle;
                            disp('Found particle file - loaded it')
                        end
    
                        if run < 4
                            currentTrackMov.superResolve(q);
                            currentTrackMov.correctDrift;
                        end
                    else
                        currentTrackMov.ConvCandToPart;
                        disp('Not tracking NPs but irregular shapes - no fitting')
                    end

                    if run < 5
                        currentTrackMov.trackParticle(trackParam,q);
                    else
                        load(filename5);
                        currentTrackMov.traces3D = TrackedData;
                        disp('Found trace file - loaded it')
                    end
                        
                    [traces] = currentTrackMov.getTraces;
    
                    allTraces = [];
                    fileN = cell(length(traces),1);
                    fileN(:,1) = {i};
               
                    [xStep,xMotor] = currentTrackMov.getXPosMotor;
                    [yStep,yMotor] = currentTrackMov.getYPosMotor;
                    [zSt,zMotor]   = currentTrackMov.getZPosMotor;
    
                    colMot = cell(length(traces),1);
                    colMot(:,1) = {xMotor};
                    colStep = cell(length(traces),1);
                    colStep(:,1) = {xStep};
    
                    rowMot = cell(length(traces),1);
                    rowMot(:,1) = {yMotor};
                    rowStep = cell(length(traces),1);
                    rowStep(:,1) = {yStep};
    
                    zMot = cell(length(traces),1);
                    zMot(:,1) = {zMotor};
                    zStep = cell(length(traces),1);
                    zStep(:,1) = {zSt};
    
                    allTraces = [allTraces; traces(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot ];
                    currentTrackMov.traces3D = [traces(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot];
                    obj.traces3D = allTraces;
                    obj.trackMovies.(fieldsN{i}) = currentTrackMov;
           
                    filename = [currentTrackMov.raw.movInfo.Path filesep 'traces3D_', num2str(q), '.mat'];
                    save(filename,'traces');
                       
                catch
                    disp(['Tracking in file ' num2str(i) ' / ' num2str(nfields) ' failed']);
                end

            end
            
            
            filename = [obj.path filesep 'traces3D_', num2str(q), '.mat'];
            save(filename,'allTraces');
            
            
            disp('=================> DONE <===================');
        end

        function MakeMovie(obj, sizeParticles, minSize, trailing, frameRate)
            MovNames = fieldnames(obj.trackMovies);

            for i = 1:numel(MovNames)
                CurrentMovie = obj.trackMovies.(MovNames{i});
                Traces = CurrentMovie.traces3D;
               
                TRAIL_COLOR  = [0 0 1];       % black traces
                DOT_COLOR    = [1 0 0];       % red dot at current position
                LINE_WIDTH   = 1.5;
                DOT_SIZE     = 0.5;
                FRAME_RATE   = 33;
                OUTPUT_FILE  = append(CurrentMovie.raw.movInfo.Path, filesep, 'TraceMovie');
                MIN_LENGTH   = 20;            % minimum track length to display
             
                trackLengths = cellfun(@(t) size(t, 1), Traces(:, 1));
                Traces = Traces(trackLengths > MIN_LENGTH, :);
                n_tracks = size(Traces, 1);
                
                % --- Load microscopy movie --------------------------------
                nFrames  = CurrentMovie.calibrated{1,1}.nFrames;
     
                startcal = strfind(CurrentMovie.calibrated{1,1}.filePath.plane2, 'calibrated');
                filePath = append(CurrentMovie.raw.movInfo.Path, filesep, CurrentMovie.calibrated{1,1}.filePath.plane2(startcal(1):end));
                
                % safer explicit load:
                frame1 = Load.Movie.tif.getFrame(filePath, 1);
                [img_h, img_w] = size(frame1);
                Movie = zeros(img_h, img_w, nFrames, 'like', frame1);
                Movie(:,:,1) = frame1;
                for f = 2:nFrames
                    Movie(:,:,f) = Load.Movie.tif.getFrame(filePath, f);
                end
                fprintf('Movie loaded.\n');
                
                % --- Extract (frame, x, y) per trace from Traces ----------
                % Traces(:,1) -> each cell is an n x 12 table
                % col 3 = rowM (x), col 4 = colM (y), col 9 = t (frame)
                track_data = cell(n_tracks, 1);
                f_min_all  = Inf;
                f_max_all  = -Inf;
                
                for k = 1:n_tracks
                    tbl   = Traces{k, 1};              % n x 12 table
                    frames = tbl{:, 9};               % column 't'
                    xp     = tbl{:, 3};               % column 'rowM'
                    yp     = tbl{:, 4};               % column 'colM'
                
                    [frames_s, ord] = sort(frames);
                    track_data{k}   = [frames_s, xp(ord), yp(ord)];
                
                    f_min_all = min(f_min_all, frames_s(1));
                    f_max_all = max(f_max_all, frames_s(end));
                end
                
                f_start    = max(1, f_min_all);
                f_end      = min(nFrames, f_max_all);
                n_vid_frames = f_end - f_start + 1;
                
                fprintf('Rendering frames %d to %d (%d total)\n', f_start, f_end, n_vid_frames);
                
                % --- Figure setup -----------------------------------------
                fig = figure('Name', 'Nanoparticle Traces', ...
                             'Color', 'k', ...
                             'Units', 'pixels', ...
                             'Position', [80, 80, img_w, img_h], ...
                             'MenuBar', 'none', ...
                             'ToolBar', 'none');
                
                ax = axes('Parent', fig, ...
                          'Units', 'pixels', ...
                          'Position', [0, 0, img_w, img_h], ...
                          'Color', 'k', ...
                          'XColor', 'none', ...
                          'YColor', 'none');
                axis(ax, 'image', 'off');
                hold(ax, 'on');
                
                % --- Normalise first frame for display --------------------
                img_norm = mat2gray(double(Movie(:,:,f_start)));
                h_img    = imagesc(ax, img_norm);
                colormap(ax, gray);
                axis(ax, 'image', 'off');
                hold(ax, 'on');
                
                % --- Pre-create one line + one dot handle per trace -------
                h_line = gobjects(n_tracks, 1);
                h_dot  = gobjects(n_tracks, 1);
                
                for k = 1:n_tracks
                    h_line(k) = plot(ax, NaN, NaN, '-', ...
                                     'Color',     TRAIL_COLOR, ...
                                     'LineWidth', LINE_WIDTH);
                    h_dot(k)  = plot(ax, NaN, NaN, 'o', ...
                                     'Color',          DOT_COLOR, ...
                                     'MarkerFaceColor', DOT_COLOR, ...
                                     'MarkerSize',      DOT_SIZE, ...
                                     'LineWidth', 0.5);
                end
                
                % Frame counter label
                h_txt = text(ax, img_w * 0.02, img_h * 0.04, '', ...
                             'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', ...
                             'VerticalAlignment', 'top', 'Interpreter', 'none');
                
                % --- VideoWriter setup ------------------------------------
                vid_out           = VideoWriter(OUTPUT_FILE, 'MPEG-4');
                vid_out.FrameRate = FRAME_RATE;
                vid_out.Quality   = 92;
                open(vid_out);
                
                % --- Render loop ------------------------------------------
                for f = f_start:f_end
                    vid_idx = f - f_start + 1;
                
                    % Update background (normalised to [0 1] per-frame for contrast)
                    frame_raw = double(Movie(:,:,f));
                    set(h_img, 'CData', mat2gray(frame_raw));
                
                    % Update each trace up to current frame
                    for k = 1:n_tracks
                        td   = track_data{k};          % [frame, x, y]
                        mask = td(:,1) <= f;
                
                        if sum(mask) < 1
                            set(h_line(k), 'XData', NaN, 'YData', NaN);
                            set(h_dot(k),  'XData', NaN, 'YData', NaN);
                        else
                            xp_tr = td(mask, 2);
                            yp_tr = td(mask, 3);
                            set(h_line(k), 'XData', yp_tr,      'YData', xp_tr);       % colM=x-axis, rowM=y-axis
                            set(h_dot(k),  'XData', yp_tr(end), 'YData', xp_tr(end));
                        end
                    end
                
                    set(h_txt, 'String', sprintf('frame %d', f));
                    drawnow limitrate;
                
                    writeVideo(vid_out, getframe(fig));
                
                    if mod(vid_idx, 50) == 0
                        fprintf('  ... %d / %d frames done\n', vid_idx, n_vid_frames);
                    end
                end
                
                close(vid_out);
                close(fig);
                fprintf('Done! Video saved to: %s\n', OUTPUT_FILE);
            end
        end
    
    
        function cleanedChannel = cleanTraces(obj, channel, threshold, ch)
                numTraces = size(channel, 1);
                merged = false(numTraces, 1); % Keep track of merged traces
                cleanedChannel = {};
                
                h = waitbar(0, 'initializing channel');
                for i = 1:numTraces
                    waitbar(i./numTraces,h, append('Cleaning traces channel ', num2str(ch)));
                    if merged(i)
                        continue
                    else
                        trace1 = channel{i, 1};
                        coords1 = table2array(trace1(:, 1:2)); 
                        time1 = table2array(trace1(:, 10));
                        intensity1 = trace1.intensity;
                        
                        for j = i+1:numTraces
                            if merged(j)
                                continue
                            else                    
                                trace2 = channel{j, 1};
                                coords2 = table2array(trace2(:, 1:2)); 
                                time2 = table2array(trace2(:, 10));
                                intensity2 = trace2.intensity;
    
                                common_time = intersect(time1, time2);
                                if ~isempty(common_time)
                                    minTimeGap = 0;
                                else
                                    minTimeGap = min(abs(time1 - min(time2)));
                                    minTimeGap = min(minTimeGap, min(abs(time2 - min(time1))));
                                end
                    
                                distances = sqrt(sum((mean(coords1) - mean(coords2)).^2, 2));
                        
                                % Check if they are within threshold distance
                                if distances < threshold && minTimeGap < 75
    
                                    % Merge traces
                                    merged(j) = true;
                                    common_time = intersect(trace1.t, trace2.t);
                                    merged_trace = [trace1; trace2];
                                    for o = 1:size(common_time)
                                        idx = find(merged_trace.t == common_time(o));
                                        [~, Idx] = max(merged_trace.intensity(idx));
                                        idx(Idx) = [];
                                        merged_trace(idx,:) = [];
                                    end      
                                    merged_trace = sortrows(merged_trace, 10);
                                    
                                    trace1 = merged_trace; % Update trace1 with new merged data
                                end
                            end
                        end
                        cleanedChannel{end+1, 1} = trace1; % Store merged trace
                    end
                end
                close(h)
                
        end


         %Plotting for individual movies
        function showLoc(obj,idx)
             fieldsN = fieldnames(obj.trackMovies);
             maxIdx = length(fieldsN);
             assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
             
             currentMov = obj.trackMovies.(fieldsN{idx});
             
             currentMov.showCorrLoc;
        end
        
        % function showTraces(obj,idx)
        %      fieldsN = fieldnames(obj.trackMovies);
        %      maxIdx = length(fieldsN);
        %      assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
        % 
        %      currentMov = obj.trackMovies.(fieldsN{idx});
        % 
        %      currentMov.showTraces;
        % end

        function showTraces(obj)
            figure()
            for j = 1:length(obj.traces3Dcommon)
                traces = obj.traces3Dcommon{j,1};
                
                subplot(round(j./2), 2, j)
                hold on 
                for i = 1: length(traces)
                    currentTrace = traces{i};
                    plot3(currentTrace.col, currentTrace.row, currentTrace.z)
                    
                end
                hold on
            end
        end
        
        function evalAccuracy(obj,dim,idx)
            
            fieldsN = fieldnames(obj.trackMovies);
            maxIdx = length(fieldsN);
            assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
            
            currentMov = obj.trackMovies.(fieldsN{idx});
            
            currentMov.evalAccuracy(dim);
            
        end
        
        function [int,SNR] = getAvgIntensity(obj)
            if ~isempty(obj.traces3Dcommon)
                assert(~isempty(obj.traces3Dcommon),'You need to extract 3D traces before extracting average intensity');
                for j = 1:length(obj.traces3Dcommon)
                    traces = obj.traces3Dcommon{j,1};
                    nTraces = length(traces);
                    int = zeros(nTraces,1);
                    SNR = zeros(nTraces,1);
                    for i = 1: length(traces)
                        currentTrace = traces{i};
                        intensity(i) = mean(currentTrace.intensity);
                        SNRatio(i) = mean(currentTrace.SNR);
                        
                    end
                    
                    int(j,1) = mean(intensity);
                    SNR(j,1) = mean(SNRatio);
                end
            elseif ~isempty(obj.traces3D)
                 for j = 1:length(obj.traces3D)
                    traces = obj.traces3D{j,1};
                    nTraces = length(traces);
                    int = zeros(nTraces,1);
                    SNR = zeros(nTraces,1);
                    for i = 1: length(traces)
                        currentTrace = traces{i};
                        intensity(i) = mean(currentTrace.intensity);
                        SNRatio(i) = mean(currentTrace.SNR);
                        
                    end
                    
                    int(j,1) = mean(intensity);
                    SNR(j,1) = mean(SNRatio);
                end
            else
                assert(~isempty(obj.traces3Dcommon),'You need to extract 3D traces before extracting average intensity');
            end
            
        end
        
        function [msd,traces] = getMSD(obj,dimension)
            
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting RMSD');
            
            traces = obj.traces3D;
            
            switch nargin
                case 1
                    
                    dimension = '3D';
                    
                case 2
                    
                otherwise
                    
                    error('too many input arguments');
                    
            end
            msd = cell(traces);
            MSDmat = zeros(obj.trackMovies.('mov1').raw.movInfo.maxFrame(1),length(traces));
            for i = 1 : size(traces,1)
                
                currentTrace = traces{i,:};
                
                if size(traces{i},1)>1
                    coord = [currentTrace.col,currentTrace.row,currentTrace.z];
                    [msdTmp,~] = MSD.calc(coord,dimension);
                    currentTrace.MSD = zeros(size(currentTrace.row,1),1);
                    currentTrace.MSD(1:end-1) = msdTmp;
                    traces{i} = currentTrace;
                    msd{i} = msdTmp;
                    MSDmat(1:length(msdTmp),i) = msdTmp(:);
                end
               
            end 
            sizes = cellfun(@size,traces,'UniformOutput',false);
            idxMat   = cellfun(@(x) x==11,sizes(:,1), 'UniformOutput', 0);
            idx = cellfun(@sum,idxMat,'UniformOutput',1);
            %delete traces where no MSD was calculated 
            traces(logical(~idx),:) = [];
            msd(logical(~idx),:) = [];
            obj.MSD = msd;
            obj.traces3D = traces;    
            
        end
        
        function saveData(obj,q)
            
            trackRes = struct; 
            disp('Saving Data');

            if ~isempty(obj.traces3D)
                trackData = obj.traces3D;
                trackRes.traces = trackData; 
                trackRes.info = obj.info;
                trackRes.path = obj.path;
                name = append('trackResults',num2str(q),'.mat');
                filename = [obj.path filesep name];
                save(filename,'trackRes');
                disp(append('Data of channel ', num2str(q), ' succesfully saved'));  
            else
                
                warning('No Data was saved because no traces or MSD could be found, please make sure you ran the analysis first');     
            end   
        end
        
        function showMSD(obj)
            MSD = obj.MSD;
            
            figure()
            hold on
            
            for i = 1:length(MSD)
                currentMSD = MSD{i};
                plot(currentMSD)
            end
                       
        end
     end
end
