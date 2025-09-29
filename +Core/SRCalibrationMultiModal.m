classdef SRCalibrationMultiModal < handle

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
        function obj = SRCalibrationMultiModal(folder2Data,cal2D,info, varargin)
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.info = info;
            obj.cal2D = cal2D;

            obj.MoviesCh1 = Core.SRCalibration(folder2Data, cal2D, obj.info);
            obj.MoviesCh2 = Core.SRCalibration(folder2Data, cal2D, obj.info);
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

        function retrieveSRCalMov(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            count = 0;
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                %If element i is a folder
                if folder2Mov(i).isdir
                    count = count+1;
                    folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                    file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                    if ~isempty(file2Analyze)    
                        file.path = file2Analyze.folder;
                        file.ext  = obj.ext;
                        
                        tmp = Core.MPSRCalMovie(file, obj.cal2D,obj.info);
                        tmp.calibrate;
                        
                        % if count == 1
                        %     tmp.giveInfo;
                        % else
                        %     tmp.info = obj.SRCalMovies.(['SRCal' num2str(1)]).getInfo; 
                        % end
                        %we extract z motor position to check if the movie
                        %is indeed a zCalibration (expect zStack)
                        
                        [zStep, ~] = tmp.getZPosMotor;
                        %TODO: Check other motor position (we do not want
                        %any other movement here.
                        
                        if zStep ~= 0
                            %if it is we store
                            Movie1 = Core.MPSRCalMovie(file , obj.MoviesCh1.cal2D, obj.MoviesCh1.info);
                            Movie1.calibrated = tmp.calibrated{1,1};
                            if count == 1
                                Movie1.giveInfo;
                            else
                                Movie1.info = obj.MoviesCh1.SRCalMovies.(['SRCal' num2str(1)]).getInfo; 
                            end
                            obj.MoviesCh1.SRCalMovies.(['SRCal' num2str(i-2)]) = Movie1;

                            Movie2 = Core.MPSRCalMovie(file , obj.MoviesCh2.cal2D, obj.MoviesCh2.info);
                            Movie2.calibrated = tmp.calibrated{1,2};
                            if count == 1
                                Movie2.giveInfo;
                            else
                                Movie2.info = obj.MoviesCh2.SRCalMovies.(['SRCal' num2str(1)]).getInfo; 
                            end
                            obj.MoviesCh2.SRCalMovies.(['SRCal' num2str(i-2)]) = Movie2;
                            
                        else
                            %if it is not we throw a warning message as it
                            %might be that many movie are in the main
                            %folder
                            error(['In ' folder2Mov(i).folder filesep folder2Mov(i).name ' no movement of the Z motor was found, the file is therefore ignored']);
                            
                        end
                        
                        
                    else
                        
                        warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);
                        
                    end
                end
            end
        end

        function SRAnalysis(obj, detectParam, trackParam)
            detectParam1 = detectParam{1};
            detectParam2 = detectParam{2};
            obj.MoviesCh1.retrieveSRCalData(detectParam1,trackParam,1);
            obj.MoviesCh2.retrieveSRCalData(detectParam2,trackParam,2);

            % calibratewwwnBB
            %% calc translation
            refPlane = 4;
            obj.MoviesCh1.corrTranslation(refPlane,1);
            obj.MoviesCh1.checkAccuracy(refPlane);      
            %% calc rotation
            obj.MoviesCh1.corrRotation(refPlane,1);
            obj.MoviesCh1.checkAccuracy(refPlane);

            obj.MoviesCh2.corrTranslation(refPlane,2);
            obj.MoviesCh2.checkAccuracy(refPlane);      
            %% calc rotation
            obj.MoviesCh2.corrRotation(refPlane,2);
            obj.MoviesCh2.checkAccuracy(refPlane);

            obj.CalcAccuracyChannels(refPlane);
        end

        function [rotMat,corrData] = CalcAccuracyChannels(obj,refPlane)
                nPlanes = size(obj.calib.SRCorrData{1, 1}, 1);
                figure()
                for p = 1:nPlanes
                    PartCh1 = [];
                    PartCh2 = [];
                    coords1s = [];
                    coords2s = [];
                    coords2t = [];
                    Im1 = [];
                    Im2 = [];
                    Im2Moved = [];

                    PartCh1 = obj.calib.SRCorrData{1, 1}{p, 1};
                    PartCh2 = obj.calib.SRCorrData{2, 1}{p, 1};

                    idx = find((diff(PartCh2.partNum)) < 0);
                    for z = 1:(size(idx, 1)+1)
                        if z == 1
                            Start = 1;
                        else
                            Start = idx(z-1)+1;
                        end

                        if z == size(idx, 1)+1
                            End = size(PartCh1, 1);
                        else
                            End = idx(z);
                        end
                        
                        Part1Selected = PartCh1(Start:End, :);
                        Part2Selected = PartCh2(Start:End, :);

                        coords1 = [Part1Selected.row, Part1Selected.col];
                        coords1 = sortrows(coords1);
                        coords2 = [Part2Selected.row, Part2Selected.col];
                        coords2 = sortrows(coords2);

                        D = pdist2(coords1, coords2);
                        threshold = 30;
                        [matches, costs] = matchpairs(D, threshold);
                        matched_coords1 = coords1(matches(:,1), :);
                        matched_coords2 = coords2(matches(:,2), :);
                                           
                        [~, ~, transform] = procrustes(matched_coords1, matched_coords2);
                        
                        coords2t = transform.b*matched_coords2*transform.T + transform.c;
    
                        subplot(2, nPlanes, p)
                        scatter(coords1(:,2), coords1(:,1),"MarkerEdgeColor", "#0072BD")
                        hold on 
                        scatter(coords2(:,2), coords2(:,1),"MarkerEdgeColor", "#D95319")
                        hold on
                        title(append('Plane ', num2str(p), ' raw'))
    
                        subplot(2, nPlanes, p+8)
                        scatter(matched_coords1(:,2), matched_coords1(:,1),"MarkerEdgeColor", "#0072BD")
                        hold on 
                        scatter(coords2t(:,2), coords2t(:,1),"MarkerEdgeColor", "#D95319")
                        hold on
                        title(append('Plane ', num2str(p), ' transf'))
       
                        Transformations{p,z} = transform;
                    end
                    
                end
                sgtitle('Transformations Channel1 - Channel2')

                for p = 1:nPlanes
                    for i = 1:z
                        T(:,:, i) = Transformations{p, i}.T;
                        b(i) = Transformations{p, i}.b;
                        c(:,:,i) = Transformations{p, i}.c(1,:);
                    end
                    Transf{p,1}.T = mean(T,3);
                    Transf{p,1}.b = mean(b);
                    Transf{p,1}.c = mean(c, 3);

                    Transf{p,2}.T = Transf{p,1}.T';
                    Transf{p,2}.b = 1/(Transf{p,1}.b);
                    Transf{p,2}.c = -(1 /Transf{p,1}.b) * Transf{p,1}.c * Transf{p,1}.T';
                end

                Transf = array2table(Transf, "VariableNames", {'Coords2toCoords1', 'Coords1toCoords2'});

                figure()
                for p = 1:nPlanes
                    PartCh1 = [];
                    PartCh2 = [];
                    coords1s = [];
                    coords2s = [];
                    coords2t = [];

                    PartCh1 = obj.calib.SRCorrData{1, 1}{p, 1};
                    PartCh2 = obj.calib.SRCorrData{2, 1}{p, 1};

                    coords1 = [PartCh1.row, PartCh1.col];
                    coords1 = sortrows(coords1);
                    coords2 = [PartCh2.row, PartCh2.col];
                    coords2 = sortrows(coords2);
                        
                    coords2t = Transf.Coords2toCoords1{1}.b*coords2*Transf.Coords2toCoords1{1}.T + Transf.Coords2toCoords1{1}.c;
    
                    subplot(2, nPlanes, p)
                    scatter(coords1(:,2), coords1(:,1),"MarkerEdgeColor", "#0072BD")
                    hold on 
                    scatter(coords2(:,2), coords2(:,1),"MarkerEdgeColor", "#D95319")
                    hold on
                    title(append('Plane ', num2str(p), ' raw'))
    
                    subplot(2, nPlanes, p+8)
                    scatter(coords1(:,2), coords1(:,1),"MarkerEdgeColor", "#0072BD")
                    hold on 
                    scatter(coords2t(:,2), coords2t(:,1),"MarkerEdgeColor", "#D95319")
                    hold on
                    title(append('Plane ', num2str(p), ' transf'));                    
                end


                obj.calib.corr{2,1}.Transformations = Transf;
                SRCal = obj.calib.corr{2,1};
                fileName = sprintf('%s%sSRCalibration2.mat',obj.path,'\');
                save(fileName,'SRCal');
        end
    end
end

