classdef SRCalibration < handle
    %SRCalibration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        ext
        info
        SRCalMovies
        cal2D
        calib
        accuracy
    end
    
    methods
        
        function obj = SRCalibration(path2SRCal,cal2D,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2SRCal.path;
            obj.ext  = path2SRCal.ext;
            obj.cal2D = cal2D;
            obj.info = info;
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
                obj.cal2D = cal2D;
                
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
                        
                        if count == 1
                            tmp.giveInfo;
                        else
                            tmp.info = obj.SRCalMovies.(['SRCal' num2str(1)]).getInfo; 
                        end
                        %we extract z motor position to check if the movie
                        %is indeed a zCalibration (expect zStack)
                        
                        [zStep, ~] = tmp.getZPosMotor;
                        %TODO: Check other motor position (we do not want
                        %any other movement here.
                        
                        if zStep ~= 0
                            %if it is we store
                            obj.SRCalMovies.(['SRCal' num2str(i-2)]) = tmp;
                            
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
        
        function retrieveSRCalData(obj,detectParam, trackParam)
            MultiModal = obj.info.multiModal + 1;
            fieldsN = fieldnames(obj.SRCalMovies);
            try
                nPlanes = obj.SRCalMovies.(fieldsN{1}).calibrated.nPlanes;
            catch
                nPlanes = obj.SRCalMovies.(fieldsN{1}).calibrated{1,1}.nPlanes;
            end
            %Checking user input
            assert(nargin==3, 'retrieveSRCalData expects 3 inputs, 1)detection Parameters, fit z parameter, tracking parameter');
            % assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            % assert(and(isstruct(trackParam),and(isfield(trackParam,'euDistPx'),isfield(trackParam,'commonPlanes'))),'trackParam is expected to be a struct with 2 fields "euDistPx", the tolerated euclidian distance between consecutive frames, "commonPlanes" nPlane to be considered as frame consistent ([1 9])');
            
            %Extraction of SRData
            nfields = numel(fieldnames(obj.SRCalMovies));
            for q = 1:MultiModal
                allData{q,1} = cell(nPlanes-1,2);
                calPerPlane{q,1} = cell(nPlanes,2);
            end
            for i = 1: nfields
                disp(['Retrieving data from SRCal file ' num2str(i) ' / ' num2str(nfields) ' ...']);
               
                currMov = obj.SRCalMovies.(fieldsN{i});
                %Molecule detection
                currMov.findCandidatePos(detectParam);

                 %SR fitting
                currMov.SRLocalizeCandidate(detectParam);
                
                %plane consolidation
                try
                    frames = 1:currMov.calibrated.nFrames;
                catch
                    frames = 1:currMov.calibrated{1,1}.nFrames;
                end
                currMov.consolidatePlanes(frames,detectParam);
                
                %getting Calibration data
                [SRCalibData,dataPerPlane] = currMov.getSRCalData(trackParam);
                
                %Get data for every particle every planes together stored
                %in allData.

                for q = 1:MultiModal      
                    for j = 1:length(allData{q,1})
                        allData{q,1}{j,1} = [allData{q,1}{j,1} ; SRCalibData{q,1}{j} ];
                        if ~isempty(SRCalibData{q,1}{j})
                            allData{q,1}{j,2} = [allData{q,1}{j,2} ; ones(height(SRCalibData{q,1}{j,1}),1)*i];
                        end
                    end
                    
                    for j = 1:length(dataPerPlane{q,1})
                        calPerPlane{q,1}{j,1}   = [calPerPlane{q,1}{j,1}; dataPerPlane{q,1}{j}];
                        if ~isempty(dataPerPlane{q,1}{j})
                            calPerPlane{q,1}{j,2}   = [calPerPlane{q,1}{j,2}; ones(height(dataPerPlane{q,1}{j,1}),1)*i];
                        end
                    end
                end
               
            end
            
            %reorder allData
            
            disp('=================> DONE <===================');
            %store the results of the data extraction
            obj.calib.trackParam = trackParam;
            obj.calib.data = allData;
            obj.calib.dataPerPlane = calPerPlane;
            
        end
        
        function [transMat,corrData] = corrTranslation(obj,refPlane)
            for q = 1:obj.info.multiModal+1
                assert(~isempty(obj.calib.data{q,1}),'You need to extract the SR calibration data before correction for translation');
                SRCalibData = obj.calib.data{q,1};
                data2Corr = obj.calib.dataPerPlane{q,1};
                %Calculate the translation
                [transMat] = Core.MPSRCalMovie.getTrans(SRCalibData);
                
                %Correct the data
                  %Correct the data for rotation
                  corrData = cell(size(data2Corr));
                for i = 1:length(data2Corr)
                    currentData = data2Corr{i};
                    [cData] = Core.MPSRCalMovie.applyTrans(currentData,transMat,refPlane);
                    corrData{i} = cData;
                end
                
                %storing
                obj.calib.SRCorrData{q,1} = corrData;
                obj.calib.corr{q,1}.trans = transMat;
                
                %Saving
                SRCal = obj.calib.corr{q,1};
                if q == 1
                   fileName = sprintf('%s%sSRCalibration1.mat',obj.path,'\');
                elseif q == 2
                   fileName = sprintf('%s%sSRCalibration2.mat',obj.path,'\');
                end
                save(fileName,'SRCal');
            end
        end
        
        function [rotMat,corrData] = corrRotation(obj,refPlane)
            for q = 1:obj.info.multiModal+1
            
                SRCalibData = obj.calib.data{q,1};
                data2Corr = obj.calib.dataPerPlane{q,1};
                %Calculate the rotation
                [rotMat] = Core.MPSRCalMovie.getRot(SRCalibData);
                
                corrData = cell(length(data2Corr),1);
                 %Correct the data for rotation
                for i = 1:length(data2Corr)
                    currentData = data2Corr{i};
                    [cData] = Core.MPSRCalMovie.applyRot(currentData,rotMat,refPlane);
                    corrData{i} = cData;
                end
    
                obj.calib.SRCorrData{q,1} = corrData;
                obj.calib.corr{q,1}.rot = rotMat;
    
                %Saving
                SRCal = obj.calib.corr{q,1};
                if q == 1
                   fileName = sprintf('%s%sSRCalibration1.mat',obj.path,'\');
                elseif q == 2
                   fileName = sprintf('%s%sSRCalibration2.mat',obj.path,'\');
                end
                save(fileName,'SRCal');
            end  
        end

        function [rotMat,corrData] = CalcAccuracyChannels(obj,refPlane)
            if obj.info.multiModal == 1
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
                        
                    coords2t = Transf.Coords2toCoords1{p}.b*coords2*Transf.Coords2toCoords1{p}.T + Transf.Coords2toCoords1{p}.c;
    
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
                
            else
                disp("No multimodal selected - no shift between channels to calculate")
            end
           
        end

        function applySRCal(obj, refPlane)
            for i = 1:numel(obj.SRCalMovies)
                currentMovie = obj.SRCalMovies.(append('SRCal', num2str(i)));
                for q = 1:(obj.info.multiModal + 1)
                    currentMovie.unCorrLocPos{q, 1}
                end
            end
        end

        
        function checkAccuracy(obj,refPlane)
            for q = 1:obj.info.multiModal+1
                data = obj.calib.dataPerPlane{q,1};
                corrData = obj.calib.SRCorrData{q,1};
                
                errRow = cell(length(data)-1,1);
                errCol = errRow;
                euclDist = errRow;
                errCRow = errRow;
                errCCol = errRow;
                euclCDist= errRow;
                
                allDist = [];
                allCorrDist = [];
                %Calculate the error on euclidian distance before and after
                %correction for each individual beads used for calibration
                for i = 1:length(data)-1
                    %Here we compare the data of plane x with ellipticity <1
                    %with the data of plane x+1 with ellipticity >1 before and
                    %after correction.
                    
                    mData = data{i};
                    nData = data{i+1};
    %                 mIdxMov = data{i,2};
    %                 nIdxMov = data{i+1,2};
    %                 
    %                 mTestIdx = mData.frame.^3 .* mData.partNum.^3 .*mIdxMov.^3;
    %                 nTestIdx = nData.frame.^3 .* nData.partNum.^3 .*nIdxMov.^3;
    %                 
    %                 idxM = ismember(mTestIdx,nTestIdx);
    %                 idxN = ismember(nTestIdx, mTestIdx);
                    if size(mData,1)<size(nData,1)
                        idxM = true(size(mData,1),1);
                        idxN = true(size(nData,1),1);
                        idxN(2:2:end) = false;
                    elseif size(mData,1)==size(nData,1)
                        
                        idxM = true(size(mData,1),1);
                        idxM(1:2:end) = false;
                        idxN = true(size(nData,1),1);
                        idxN(2:2:end) = false;
                    else
                        idxM = true(size(mData,1),1);
                        idxN = true(size(nData,1),1);
                       
                        idxM(1:2:end) = false;
                    end
    
    
                    errRow{i} = mData.row(idxM) - nData.row(idxN);
                    errCol{i} = mData.col(idxM) - nData.col(idxN);
                    euclDist{i} = sqrt(errRow{i}.^2+errCol{i}.^2);
                    
                    cData = corrData{i};
                    dData = corrData{i+1};
                    
                    errCRow{i}   =  cData.row(idxM) - dData.row(idxN);
                    errCCol{i}   =  cData.col(idxM) - dData.col(idxN);
                    euclCDist{i} = sqrt(errCRow{i}.^2+errCCol{i}.^2);
                    
                    allDist = [allDist; mData.row mData.col];
                    allCorrDist = [allCorrDist; cData.row cData.col];
                    
                end
                
                %Plot histogram of errors
                figure
                fullErrRow = [];
                fullCErrRow = [];
                errorPerPlane = zeros(length(data)-1,2);
                errorCorrPerPlane = zeros(length(data)-1,2);
                for i = 1:length(data)-1
                    [errN, errEdge] = histcounts(errRow{i},[-2:0.1:2]);
                    errBin = errEdge(1:end-1)+(errEdge(2)-errEdge(1));
                    [errCN, errCEdge] = histcounts(errCRow{i},[-2:0.1:2]);
                    errCBin = errCEdge(1:end-1)+(errCEdge(2)-errCEdge(1));
                    
                    subplot(1,3,1)
                    hold on
                    bar(errBin,errN)
                    xlim([-2 2])
                    title('Before SR Correction');
                    xlabel('X Error in Pixel');
                    ylabel('Occurennce');
                    legend;
                    hold off
                    subplot(1,3,2)
                    hold on
                    bar(errCBin,errCN)
                    xlim([-2 2])
                    title('After SR Correction');
                    xlabel('X Error in Pixel');
                    ylabel('Occurennce');
                    legend;
                    hold off
                    
                    fullErrRow = [fullErrRow; errRow{i}];
                    fullCErrRow = [fullCErrRow; errCRow{i}];
                    
                    errorPerPlane(i,:) = [mean(abs(errRow{i})) mean(errRow{i})];
                    errorCorrPerPlane(i,:) = [mean(abs(errCRow{i})) mean(errCRow{i})];
                    
                end
                
                [fErrN, fErrEdge] = histcounts(fullErrRow,[-2:0.1:2]);
                fErrBin = fErrEdge(1:end-1)+(fErrEdge(2)-fErrEdge(1));
                [fErrCN, fErrCEdge] = histcounts(fullCErrRow,[-2:0.2:2]);
                fErrCBin = fErrCEdge(1:end-1)+(fErrCEdge(2)-fErrCEdge(1));
                
                %Compare all errors before and after
                subplot(1,3,3)
                hold on
                plot(fErrBin,fErrN);
                plot(fErrCBin,fErrCN);
                legend({'Before Corr','After Corr'});
                title('Distribution before/after')
                xlim([-2 2])
                hold off
                
                figure
                plane = 1:7;
                height = min(errorPerPlane(:,2)):0.01: max(errorPerPlane(:,2));
                x = ones(1,length(height))*refPlane;
                subplot(1,2,1)
                hold on
                plot(plane,errorPerPlane(:,1));
                plot(plane,errorCorrPerPlane(:,1));
                plot(x,height,'k--')
                legend({'Before Correction','After Correction', 'Ref plane'});
                xlabel('Planes')
                ylabel('error compare to reference plane')
                hold off ;
                
                subplot(1,2,2)
                hold on
                plot(plane,errorPerPlane(:,2));
                plot(plane,errorCorrPerPlane(:,2));
                plot(x,height,'k--')
                legend({'Before Correction','After Correction', 'Ref plane'});
                xlabel('Planes')
                ylabel('error compare to reference plane')
                hold off ;
                
                %display some meaningful values to know what is happening
                avgAbsErr = mean(abs(nonzeros(fullErrRow)));
                avgErr = mean(nonzeros(fullErrRow));
                stdErr = std(abs(nonzeros(fullErrRow)));
                avgAbsCErr = mean(abs(nonzeros(fullCErrRow)));
                avgCErr = mean(nonzeros(fullCErrRow));
                stdCErr = std(abs(nonzeros(fullCErrRow)));
                
                fprintf('The absolute error before correction is: %0.3f pixel\n', mean(avgAbsErr));
                fprintf('The error before correction is: %d \n', mean(avgErr));
                fprintf('The accompanying standard deviation is: %0.3f pixel \n', mean(stdErr));
                fprintf('The absolute error after correction is: %0.3f pixel\n', mean(avgAbsCErr));
                fprintf('The error after correction is: %d \n', mean(avgCErr));
                fprintf('The accompanying standard deviation is: %0.3f pixel\n', mean(stdCErr));
                
                %Let us check that the correction worked or display warning
                if or(avgAbsErr<avgAbsCErr, stdErr<stdCErr)
                    warning(['Mean error or standard deviation was smaller before'...
                        ' correction, something might have gone wrong...'])
                end
                
                figure
                subplot(1,3,1)
                plot(allDist(:,2),allDist(:,1),'k+')
                subplot(1,3,2)
                plot(allCorrDist(:,2),allCorrDist(:,1),'r+');
                subplot(1,3,3)
                hold on
                plot(allCorrDist(:,2),allCorrDist(:,1),'r+');
                plot(allDist(:,2),allDist(:,1),'k+')
                hold off
            end
        end

        function ChannelAccuracy(obj)
        end
        
    end
    
    methods (Access = private)
        
        
    end
end

