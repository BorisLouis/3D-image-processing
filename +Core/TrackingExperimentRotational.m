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
         function obj = TrackingExperimentRotational(folder2Data,cal2D,info,SRCalPath,zCalPath)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.cal2D = cal2D;
            obj.info = info;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath; 
            obj.SRCal2 = SRCalPath;
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
               

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration.mat');

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
                    obj.SRCal.path = SRCal;
                end
            end
            
        end

        function set.SRCal2(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal2.cal = SRCal;
                obj.SRCal2.path = SRCal;
            else
                assert(isfolder(SRCal), 'The given path is not a folder');

                if isfield(obj.info, 'multiModal')
                    if obj.info.multiModal == true
                        MultiModal = 1;
                    end
                elseif isfield(obj.info, 'multiTracking')
                    MultiModal = 1;
                end

                if MultiModal == true
                     [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration2.mat');

                    if isempty(file2Analyze)
                        error('No SR calibration file found in the given folder');
                    else
                        fileName = [file2Analyze.folder filesep file2Analyze.name];
                        cal = load(fileName);
                        field = fieldnames(cal);
                        cal = cal.(field{1});
                        assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                            'SR calibration is supposed to be a struct with 2 fields');
    
                        obj.SRCal2.cal = cal;
                        obj.SRCal2.path = SRCal;
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
                    
                    obj.info.multiModal = 1;
                    count = 0;
                    count = count+1;
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    tmp = Core.MPTrackingMovieRotational(file , obj.cal2D, obj.info, obj.SRCal.path, obj.ZCal.path);
                    
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
                if isempty(currentTrackMov.ROI)
                    currentTrackMov.findCandidatePos(detectParam);
                    f = 1;
                else
                    f = 2;
                end
                
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
        
       
        % Select particles only when visible at the same time in both
        % channels:
       

        function ConsolidateChannels3(obj)
            channel1 = obj.traces3D{1,1};
            channel2 = obj.traces3D{2,1};
            numTraces1 = size(channel1, 1);
            numTraces2 = size(channel2, 1);
        
            distances = zeros(numTraces1, numTraces2);
        
            for i = 1:numTraces1
                trace1 = channel1{i,1};
                for j = 1:numTraces2
                    trace2 = channel2{j,1}; 
                    
                    coords1 = table2array(trace1(:, 1:2)); 
                    time1 = table2array(trace1(:, 10));   
                    coords2 = table2array(trace2(:, 1:2));
                    time2 = table2array(trace2(:, 10)); 
                
                    common_time = intersect(time1, time2);
                    
                    coords1_common = coords1(ismember(time1, common_time), :);
                    coords2_common = coords2(ismember(time2, common_time), :);
                
                    if isempty(common_time)
                        avg_distance = Inf;
                        return;
                    end
                
                    distance = sqrt(sum((coords1_common - coords2_common).^2, 2));
                    distances(i, j) = mean(distance);
                end
            end
            
            [~, TracesIdx] = max(size(distances));
            [closest, closest_indices] = min(distances, [], TracesIdx);

            unique_vals = unique(closest_indices);
            n = 0;

            % while numel(unique_vals) < numel(closest_indices)
            %     n = n+1;
                duplicates = unique_vals(histc(closest_indices, unique_vals) > 1);
    
                for i = 1:numel(duplicates)
                    duplicate_value = duplicates(i);
                    duplicate_indices = find(closest_indices == duplicate_value);
                    [MaxDup, MaxDupIdx] = max(closest(duplicate_indices));
                    SearchPartner = duplicate_indices(MaxDupIdx);
    
                    PossPartners = distances(SearchPartner, :);
                    PossPartners(1, duplicate_value) = Inf;
    
                    [closest(SearchPartner, 1), closest_indices(SearchPartner, 1)] = min(PossPartners);
    
                    
                end
    
                unique_vals = unique(closest_indices);
            %     if n == 3
            %         break
            %     end
            % end

            closest_indices(:,2) = closest_indices(:,1);
            closest_indices(:,1) = [1:1:size(closest_indices, 1)].';
            closest_indices(:,3) = closest;
            
            traces3Dcommon = struct([]);
            for i = 1:size(closest_indices, 1)
                if closest_indices(i,3) < 350
                    Int1 = obj.traces3D{1,1}{closest_indices(i,1),1}.intensity;
                    Int2 = obj.traces3D{2,1}{closest_indices(i,2),1}.intensity;
    
                    traces3Dcommon{end+1,1} = Int1;
                    traces3Dcommon{end,2} = Int2;
    
                    Diff = Int1-Int2;
                    Ratio = Int1./Int2;
                    Time = [0:1:size(Diff,1)-1]*obj.info.expTime;
    
                    traces3Dcommon{end,3} = Diff;
                    traces3Dcommon{end,4} = Ratio;
                    traces3Dcommon{end,5} = Time.';
                end
            end

            obj.traces3Dcommon = traces3Dcommon;
        end
        
        function RotationalCalibration(obj)
            % AngFreq = 2*obj.info.RadTime*pi./180;
            % Model = 'a+b*sin(c*x+d)';
            % 
            
            Model = 'a+b*sin(1.72*x+c)';
            Fig = figure();
            AngFreq = 1.72;
            for i = 1:size(obj.traces3Dcommon,1)

                Diff = obj.traces3Dcommon{i,3};
                % Diff = Diff{1,1};
                Ratio = obj.traces3Dcommon{i,4};
                % Ratio = Ratio{1,1};
                Time = obj.traces3Dcommon{i,5};
                %obj.traces3Dcommon{i,5} = Time;

                %%%Fit on difference plot
                Mean = mean(Diff);
                Max = max(Diff);
                Min = min(Diff);
                StartPoints = [Mean, (Max-Min)./2, AngFreq, 180];
                Lower = [Min, (((Max-Mean) + (Mean - Min))./2)./2, 0, 0];
                Upper = [Max, Max-Min, AngFreq*2, 360];
    
                %[fitDiff, gofDiff] = fit(Time.', Diff, Model,'Lower', Lower, 'Upper', Upper, 'StartPoint', StartPoints);
                [fitDiff, gofDiff] = fit(Time, Diff, Model);
                coeff = coeffvalues(fitDiff);
                HeightDiff = coeff(1);
                AmpDiff = coeff(2);
                AngFreqDiff = 1.72;
                PhaseDiff = coeff(3);
                
                % if gofDiff.rsquare > 0.5
                obj.traces3Dcommon{i, 6} = HeightDiff;
                obj.traces3Dcommon{i, 7} = AmpDiff;
                obj.traces3Dcommon{i, 8} = AngFreqDiff;
                obj.traces3Dcommon{i, 9} = PhaseDiff;
                % else
                %     obj.traces3Dcommon{i, 6} = 'Failed fit';
                %     obj.traces3Dcommon{i, 7} = 'Failed fit';
                %     obj.traces3Dcommon{i, 8} = 'Failed fit';
                %     obj.traces3Dcommon{i, 9} = 'Failed fit';
                % end

                subplot(size(obj.traces3Dcommon,1), 2, i*2-1)
                plot(fitDiff, Time, Diff)

                %%%Fit on Ratio plot
                Mean = mean(Ratio);
                Max = max(Ratio);
                Min = min(Ratio);
                StartPoints = [Mean, (Max-Min)./2, AngFreq, 180];
                Lower = [Min, (((Max-Mean) + (Mean - Min))./2)./2, 0, 0];
                Upper = [Max, Max-Min, AngFreq*2, 360];

                % [fitRatio, gofRatio] = fit(Time.', Ratio, Model,'Lower', Lower, 'Upper', Upper, 'StartPoint', StartPoints);
                [fitRatio, gofRatio] = fit(Time, Ratio, Model);
                coeff = coeffvalues(fitRatio);
                HeightRatio = coeff(1);
                AmpRatio = coeff(2);
                AngFreqRatio = 1.72;
                PhaseRatio = coeff(3);
                
                % if gofRatio.rsquare > 0.5
                    obj.traces3Dcommon{i, 10} = HeightRatio;
                    obj.traces3Dcommon{i, 11} = AmpRatio;
                    obj.traces3Dcommon{i, 12} = AngFreqRatio;
                    obj.traces3Dcommon{i, 13} = PhaseRatio;
                % else
                %     obj.traces3Dcommon{i, 10} = 'Failed fit';
                %     obj.traces3Dcommon{i, 11} = 'Failed fit';
                %     obj.traces3Dcommon{i, 12} = 'Failed fit';
                %     obj.traces3Dcommon{i, 13} = 'Failed fit';
                % end

                subplot(size(obj.traces3Dcommon,1), 2, i*2)
                plot(fitRatio, Time, Ratio)

            end

            Model = 'a*(1+b*(cos(2*(x-c))))';
            Fig2 = figure();
            for i = 1:size(obj.traces3Dcommon,1)
                Time = obj.traces3Dcommon{i,5};
                Angle = Time*0.87266;

                traceCh1 = obj.traces3Dcommon{i,1};  
                f = fit(Angle, traceCh1, Model);
                coeff = coeffvalues(f);
                Nch(i,1) = coeff(1);
                Mch(i,1) = coeff(2);
                Phich(i,1) = coeff(3);
                subplot(size(obj.traces3Dcommon,1), 2, i*2-1)
                plot(f, Angle, traceCh1)
                obj.traces3Dcommon{i, 14} = Nch;
                obj.traces3Dcommon{i, 15} = Mch;

                traceCh2 = obj.traces3Dcommon{i,2};
                f = fit(Angle, traceCh2, Model);
                coeff = coeffvalues(f);
                Nch(i,2) = coeff(1);
                Mch(i,2) = coeff(2);
                Phich(i,2) = coeff(3);
                subplot(size(obj.traces3Dcommon,1), 2, i*2)
                plot(f, Angle, traceCh2)
                obj.traces3Dcommon{i, 16} = Nch;
                obj.traces3Dcommon{i, 17} = Mch;
            end

            obj.traces3Dcommon = cell2table(obj.traces3Dcommon, 'VariableNames', {'Int1', 'Int2',...
                'Diff', 'Ratio', 'Time', 'Height - Diff', 'Amp - Diff', 'Freq - Diff', 'Phase - Diff',...
                'Height - Ratio', 'Amp - Ratio', 'Freq - Ratio', 'Phase - Ratio', 'N - ch1', 'M - ch1', 'N - ch2', 'M - ch2'});
            CommonTraces = obj.traces3Dcommon;

            Filename = append(obj.path, filesep, 'CommonTraces.mat');
            save(Filename, 'CommonTraces');
            set(Fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            Figname = append(obj.path, filesep, 'DiffRatioFig.png');
            saveas(Fig, Figname)
            set(Fig2, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            Figname2 = append(obj.path, filesep, 'ChannelFig.png');
            saveas(Fig2, Figname2)
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
        
        function saveData(obj)
            
            trackRes = struct; 
            disp('Saving Data');
            
            if ~isempty(obj.traces3Dcommon)
                
                trackData = obj.traces3Dcommon;
                MSDs = obj.MSD;
                    
                if ~isempty(MSDs)
                    trackRes.MSD = MSDs;
                end
                
                trackRes.traces = trackData; 
                trackRes.info = obj.info;
                trackRes.path = obj.path;
                filename = [obj.path filesep 'trackResultsCommonCh.mat'];
                save(filename,'trackRes');
                disp('Common channeldata was succesfully saved');
            end

            if ~isempty(obj.traces3D)
                for i = 1:length(obj.traces3D)
                    trackData = obj.traces3D{i,1};
                    % MSDs = obj.MSD;
                    % 
                    % if ~isempty(MSDs)
                    %     trackRes.MSD = MSDs;
                    % end
                    
                    trackRes.traces = trackData; 
                    trackRes.info = obj.info;
                    trackRes.path = obj.path;
                    name = append('trackResults',num2str(i),'.mat');
                    filename = [obj.path filesep name];
                    save(filename,'trackRes');
                    disp('Data were succesfully saved');
                end    
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
