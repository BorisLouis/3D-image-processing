classdef MPParticleMovie < Core.MPMovie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        
        candidatePos
        unCorrLocPos
        corrLocPos
        particles
        intData
        
    end
    
    methods
        function obj = MPParticleMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParamFull, frames)
            %Method to perform localization on each plane for each frame
            %Check if some candidate exists already in the folder (previously saved)
            for q = 1:(obj.info.multiModal+1)
                if iscell(detectParamFull)
                    if q == 1
                        detectParam = detectParamFull{1};
                    elseif q == 2
                        detectParam = detectParamFull{2};
                    end
                else
                    detectParam = detectParamFull;
                end

                switch nargin
                        case 2
                            
                            frames = 1: obj.calibrated{1,q}.nFrames;
                            disp('Running detection on every frame');
         
                            
                        case 3
                            
                            [frames] = obj.checkFrame(frames,obj.calibrated{1,q}.nFrames);
                            
                        otherwise
                            
                            error('too many inputs');
                            
                 end

                 folder = append('calibrated', num2str(q));

                 path = append(obj.raw.movInfo.Path, filesep, folder);
                 [run, candidate] = obj.existCandidate(obj.raw.movInfo.Path, '.mat');
                
                %if we only ask 1 frame we always run
                if length(frames) == 1
                    run = true;
                end
                if run
                   
                    
                    %Localization occurs here
                    assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first or load previous data');
                    assert(nargin>1,'not enough input argument or accept loading of previous data (if possible)');
                    [candidate] = obj.detectCandidate(detectParam,frames,q);
                    
                elseif ~isempty(candidate)
                else
                    %help message
                    disp('getCandidatePos is a function that detects features in a movie');
                    disp('To work, it needs to receive a structure containing 2 detection parameter:');
                    disp('delta which is the radius around which detection is performed usually 6 pixels');
                    disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                    disp('Typically between 20 and 200');
                    
                end
                %if we only ask 1 frame we do not save
                if length(frames) >1
                    %save the data
                    folder = append('calibrated', num2str(q));
                      
                    fileName = sprintf('%s%s%s%scandidatePos.mat',obj.raw.movInfo.Path,'\',folder,'\');
                    save(fileName,'candidate');
                else
                end

                obj.candidatePos{q,1} = candidate;
                obj.info.detectParam = detectParam;
            end
        end
        
        function [matched_particles] = getROIs(obj)
            for q = 1:obj.info.multiModal+1
                T = obj.candidatePos{1,1}{1,1};
                max_distance = 10;
                T.particle_count = zeros(height(T), 1);
                ROIsize = 10;
                particles = struct([]);
                          
                for i = 1:size(T,1)
                    Particle = [];
                    row_i = T.row(i);
                    col_i = T.col(i);
                    plane_i = T.plane(i);           
                    count_planes = 0;
                    Particle(1,:) = [round(row_i-round(ROIsize./2)), round(col_i-round(ROIsize./2)),...
                                     ROIsize, ROIsize,plane_i];
                    
                    for j = i+1:height(T)
                        row_j = T.row(j);
                        col_j = T.col(j);
                        plane_j = T.plane(j);
                        
                        if plane_i ~= plane_j & ~isnan(plane_i)  & ~isnan(plane_j)
                            distance = sqrt((row_i - row_j)^2 + (col_i - col_j)^2);
                            if distance <= max_distance
                                count_planes = count_planes + 1;
                                Particle(end+1, :) = [round(row_j-round(ROIsize./2)), round(col_j-round(ROIsize./2)),...
                                                        ROIsize, ROIsize,plane_j];
                                T.plane(j) = NaN;
                            end
                        end
                    end
                    T.plane(i) = NaN;
                    T.particle_count(i) = count_planes;
                    if size(Particle, 1) >= 3
                        particles{size(particles,1)+1,1} = Particle;
                    end
                end  
                ROI{q,1} = particles;
            end

            matched_particles = cell(max([size(ROI{1,1},1), size(ROI{2,1},1)]),2);
            for i = 1:size(ROI{1,1},1)
                particle1 = ROI{1,1}{i,1};
                xmin1 = particle1(:,1);
                ymin1 = particle1(:,2);
                plane1 = particle1(:,end);
                min_distance = inf;
                matched_particle = [];
                for j = 1:size(ROI{2,1},1)
                    particle2 = ROI{2,1}{j,1};
                    xmin2 = particle2(:,1);
                    ymin2 = particle2(:,2);
                    plane2 = particle2(:,end);
                    dist = [];
                    for k = 1:size(particle1,1)
                        x = xmin1(k);
                        y = ymin1(k);
                        plane = plane1(k);
                        Idx = find(plane2 == plane, 1, 'first');
                        if ~isnan(Idx)
                            x2 = xmin2(Idx);
                            y2 = xmin2(Idx);
                            dist(end+1) = sqrt((x - x2).^2 + (y-y2).^2);
                        end
                    end
                    distav = mean(dist);
                    if distav < min_distance
                        min_distance = distav;
                        matched_particle = particle2;
                    end
                end
                matched_particles{i,1} = particle1;
                matched_particles{i,2} = matched_particle;
            end      
            obj.ROI = matched_particles;
        end
            
            function [candidate] = getCandidatePos(obj, frames, q)
                %Extract the position of the candidate of a given frame
                [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                candidate = obj.candidatePos{q,1}{idx};
                
                if isempty(candidate)
                    
                    warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                    
                end
            end  
            
            function SRLocalizeCandidate(obj,roiSize,frames)
                for q = 1:(obj.info.multiModal+1)
                    assert(~isempty(obj.calibrated{1,q}),'Data should be calibrated to consolidate');
                    assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
                    assert(~isempty(obj.candidatePos{q,1}), 'No candidate found, please run findCandidatePos before consolidation');
                    folder = append('calibrated', num2str(q));

                    path = append(obj.raw.movInfo.Path, filesep, folder);
                    [run,locPos] = obj.existLocPos(path,'.mat');
                    
                    if run
                        switch nargin
        
                            case 1
                                roiSize = 6;
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running SRLocalization on every frame with ROI of 6 pixel radius');
        
                            case 2
        
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running SRLocalization on every frame');
        
                            case 3
        
                                [frames] = obj.checkFrame(frames);
        
                            otherwise
        
                                error('too many inputs');
        
                        end
                       
                        locPos = cell(size(obj.candidatePos{q,1}));
                        h = waitbar(0,'Fitting candidates ...');
                        nFrames = length(frames);
                        %Localization occurs here
                        for i = 1 : 1:nFrames
                            disp(['Fitting candidates: frame ' num2str(i) ' / ' num2str(nFrames)]);
                            idx = frames(i);
                            %#1 Extract Candidate Position for specific frame
                            [data] = obj.getFrame(idx, q);
                            if obj.info.rotational == 1
                                [frameCandidate] = obj.getCandidatePos(1,q);
                            else
                                [frameCandidate] = obj.getCandidatePos(idx, q);
                            end
                            
                            if isempty(frameCandidate)
                                
                                warning('Frame %d did not contain any candidate',idx);
                                locPos{i} = [];
                                
                            else
                                
                                locPos{i} = obj.superResLocFit(data,frameCandidate,roiSize);
                                
                            end
                            waitbar(i/nFrames,h,['Fitting candidates: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                        end
                        close(h)
                    else
                    end
                        %save the data

                    folder = append('calibrated',num2str(q));

                    fileName = sprintf('%s%s%s%sSRLocPos.mat',obj.raw.movInfo.Path,'\', folder, '\');
                    save(fileName,'locPos');
                    
                    %store in the object
                    obj.unCorrLocPos{q,1} = locPos;
                    obj.corrLocPos{q,1}   = locPos;
                end
            end
            
            function [locPos] = getLocPos(obj, frames, q)
                 %Extract the position of the candidate of a given frame
                [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
              
                locPos = obj.unCorrLocPos{q,1}{idx};
               
                if isempty(locPos)
                    
                    warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                    
                end
            end
            
            function consolidatePlanes(obj,frames,consThresh)
                for q = 1: obj.info.multiModal +1
                    %Consolidation refers to connect molecules that were localized
                    %at similar position in different plane on a single frame.
                    assert(~isempty(obj.calibrated{1,q}),'Data should be calibrated to consolidate');
                    assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
                    assert(~isempty(obj.candidatePos{q,1}), 'No candidate found, please run findCandidatePos before consolidation');
                    assert(~isempty(obj.unCorrLocPos{q,1}),'Localization needs to be performed before consolidation');
                   
                    %Check if some particles were saved already.
                    folder = append('calibrated', num2str(q));

                    path = append(obj.raw.movInfo.Path, filesep, folder);
                    [run, particle] = obj.existParticles(path, '.mat');
                    
                    if run
                        %Check the number of function input
                        switch nargin
                            case 1
                                
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running consolidation on every frame with roi of 6 pixel');
                                consThresh = 4;
                            case 2
                                [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                                consThresh = 4;                       
                            case 3
                                [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                                assert(isnumeric(consThresh),'Consolidation threshold should be numeric');
                            otherwise
                                
                                error('Something wrong with number of input');
                                
                        end
                        
                        nFrames = length(frames);
                        %allocate for storage
                        particleList = cell(1,obj.raw.maxFrame(1));
                        nParticles = zeros(1,obj.raw.maxFrame(1));
                        idx2TP = zeros(1,obj.raw.maxFrame(1));
                        h = waitbar(0,'Consolidating candidate ...');
                        
                        %Consolidation occurs here
                        for i = 1 :nFrames
                            disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                            idx = frames(i);
                            %#1 Extract localized Position for specific frame
                            [fCandMet] = obj.getLocPos(idx, q);
                            
                            if isempty(fCandMet)
                                
                                warning('Frame %d did not contain any localized positions',idx);
                                particleList{idx} = [];
                                %particleList{idx}{1} = nan(5);
                                nParticles(idx) = 0;
                                
                            else
                                  %#2 Consolidate the position of the given frame
                                  %across plane
                                  %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
                                  switch obj.info.zMethod
                                      case 'Intensity'
                                         %we do not do anythin at the moment.
                                         focusMetric = fCandMet.magX+fCandMet.magY;
                                         
                                      case '3DFit'
                                         %we do not do anythin at the moment.
                                         focusMetric = fCandMet.magX+fCandMet.magY;
                                      case 'PSFE' 
                                        [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet.ellip,fCandMet.fMetric);
                                  end
                                    %reformating to keep the same format as how the data is saved
                                    %later
                                    fCandMet.fMetric = focusMetric;
                                    
                                    %focusMetric((1-corrEllip)>0.3) = NaN;
                                    
                                    %Plane Consolidation occur here
                                    [part] = obj.planeConsolidation(fCandMet,focusMetric,consThresh, q);
        
                                    %we delete empty cells from the array
                                    idx2Empty = cellfun(@isempty,part);
                                    part(idx2Empty(:,1),:) = [];
                           
                                    particleList{idx} = part;
                                    nParticles(idx) = length(part);
                                    
                                    if ~isempty(part)
                                        idx2TP(idx) = idx;
                                    end
                                
                            end
                            waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
                        end
                        close(h);
                        
                        %#3 Storing List
                        particle.List       = particleList;
                        particle.nParticles = nParticles;
                        particle.tPoint     = nFrames;
                        particle.idx2TP     = nonzeros(idx2TP);
                        particle.Traces     = [];
                        particle.nTraces    = [];
                        
                        folder = append('calibrated', num2str(q));

                        fileName = sprintf('%s%s%s%sparticle.mat',obj.raw.movInfo.Path,'\', folder, '\');
                        save(fileName,'particle');
                        obj.particles{q,1} = particle;
                    elseif run == 0
                        obj.particles = particle;
                    end
                end
            end
            
            function [particle] = getParticles(obj,frames)
                for q = 1:obj.info.multiModal+1
                    %GetParticles
                    [idx] = obj.checkFrame(frames,obj.calibrated{1,q}.nFrames(1));
                    particle = obj.particles{q,1}.List{idx};
                    
                    if isempty(particle)
                        
                        warning('There was no particle found in this frames, check than you ran superResConsolidate method on this frame beforehand');
                        
                    end
                end
            end

            
            function showCandidate(obj,idx)
                %Display Candidate
                for q = 1:(obj.info.multiModal + 1)
                    if obj.info.rotational ~= 1
                        assert(length(idx)==1, 'Only one frame can be displayed at once');
                        [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                        assert(~isempty(obj.candidatePos{q,1}{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
                        
                        [frame] = getFrame(obj,idx,q);
                        
                        nImages = size(frame,3);
                       
                        nsFig = ceil(nImages/4);
                        
                        candidate = obj.getCandidatePos(idx,q);
                        
                    elseif obj.info.rotational == 1
                        [frame] = obj.calibrated{2,q};
                        nImages = size(frame,3);
                        nsFig = ceil(nImages/4);
                        candidate = obj.getCandidatePos(1,q);
                    end
                       
                    rowPos    = candidate.row;
                    colPos    = candidate.col;
                    planeIdx  = candidate.plane;
                    h = figure();
                    h.Name = sprintf('Frame %d',idx);
                    for i = 1:nImages
                        
                        subplot(2,nImages/nsFig,i)
                        hold on
                        imagesc(frame(:,:,i))
                        hold on
                        plot(colPos(planeIdx==i),rowPos(planeIdx==i),'g+','MarkerSize',10)
                        axis image;
                        % grid on;
                        % a = gca;
                        % a.XTickLabel = [];
                        % a.YTickLabel = [];
                        % a.GridColor = [1 1 1];
                        title({['Plane ' num2str(i)],sprintf(' Zpos = %0.3f',obj.calibrated{1,q}.oRelZPos(i))});
                        colormap('hot')
                        if obj.info.rotational == 1
                            hold on  
                            for h = 1:size(obj.ROI,1)
                                for z = 1:size(obj.ROI{h,q},1)
                                    ROI = obj.ROI{h,q}(z,:);
                                    plane = ROI(:,5);
                                    if plane == i;
                                        rectangle('Position',[ROI(idx, 2), ROI(idx,1), ROI(idx,3), ROI(idx,4)], 'EdgeColor', 'r');
                                        hold on
                                    end
                                end
                            end
                        else
                        end
                        hold on
                    end
                end
            end
            
            function showParticles(obj,idx)
                %display particles (after consolidation), On top of the
                %localization, consolidated particles are circled.
                assert(length(idx)==1, 'Only one frame can be displayed at once');
                [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                % Show Candidate
                obj.showCandidate(idx);
                
                if isempty(obj.particles)
                    
                    warning('You did not consolidate the candidate yet, please use the consolidate method before showing the particles');
                    
                else
                    
                    if isempty(obj.particles.List(idx))
                        
                        warning('The candidates of the requested frame were not consolidated yet, only showing the candidate');
                        
                    else
                        
                        roiSize = obj.particles.roiSize;
                        nParticles = obj.particles.nParticles(idx);
                        h = gcf;
                        nPlanes = obj.calibrated.nPlanes;
                        colors = rand(nParticles,3);
                        %Display circled
                        for i = 1 : nPlanes
                            subplot(2,nPlanes/2,i)
                            hold on
                            for j = 1 : nParticles
                                currPart = obj.particles.List{idx}{j};
                                if(~isempty(currPart(currPart.plane == i,:)))
                                    part2Plot = currPart(currPart.plane == i,:);
                                    plot(part2Plot.col,part2Plot.row,'o',...
                                        'LineWidth',2, 'MarkerSize',10, 'MarkerEdgeColor',colors(j,:));
                                end
                            end
                            hold off
                        end
                        %Here we display a zoom onto the particle visible on
                        %the specific frame onto the consolidated planes
                        [frame] = getFrame(obj,idx);
                        assert(isstruct(frame),'Error unknown data format, data should be a struct');
                        for i = 1:nParticles
                            
                            currPart = obj.particles.List{idx}{i};
                            %Remove rows containing NaNs
                            idx2NaN = isnan(currPart.row);
                            currPart(idx2NaN,:) = [];
                            planes = currPart.plane;
                            figure(20+i)
                            hold on
                            for j = 1 : length(planes)
                                jdx = planes(j);
                                currFrame = frame.(sprintf('plane%d',jdx));
                                ROI = EmitterSim.getROI(currPart.col(j), currPart.row(j),...
                                    roiSize, size(currFrame,2), size(currFrame,1));
                                subplot(1,length(planes),j)
                                imagesc(currFrame(ROI(3):ROI(4),ROI(1):ROI(2)));
                                title({['Particle ' num2str(i)],[ ' Plane ' num2str(jdx)]});
                                axis image
                                colormap('jet')
                                
                            end
                            hold off
                        end
                    end
                end
            end
            
            function [candidateList] = planeConsolidation(obj,candMet,focusMetric,consThresh,q)
                %Loop through all candidate of a given frame and match them
                %between frame until none can be match or all are matched.
                nPlanes = obj.calibrated{1,q}.nPlanes;
                counter = 1;
                nPart = 0;
                maxIt = size(candMet,1);
                zMethod = obj.info.zMethod;
                candidateList = cell(max(size(find(~isnan(focusMetric)))),1);
                %continue until the list is not empty
                while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                    
                    if counter> maxIt
                        
                        error('While loop ran for an unexpectedly long time, something might be wrong');
                        
                    end
                    
                    %Find candidate in best focus
                    [~,idx] = max(focusMetric);
                    currentPlane = candMet.plane(idx);
                    
                    switch nPlanes
                        case 1
                            planes2Check = [];
                            
                        otherwise
                            
                            %Check which planes are to be checked (currently 2 planes
                            %above and 2 planes below the given plane
                            planes2Check = currentPlane-nPlanes:currentPlane-1;
                            planes2Check = planes2Check(planes2Check>0);
                            planes2Check = [planes2Check currentPlane+1:currentPlane+nPlanes];
                            planes2Check = planes2Check(planes2Check<nPlanes+1);
                            
                    end
                    currentCand = candMet(idx,:);
                    direction = -1;%Start by checking above
                    
                    particle = array2table(nan(nPlanes,size(currentCand,2)));
                    particle.Properties.VariableNames = currentCand.Properties.VariableNames;
                    
                    particle(currentCand.plane,:) = currentCand;
                    nCheck = length(planes2Check);
                    camConfig = obj.calibrated{1,q}.camConfig;
                    for i = 1:nCheck
                        
                        cand = candMet(candMet.plane == planes2Check(i),:);
                        if(planes2Check(i) > currentPlane)
                            direction = +1;%check below (Plane 1 is the uppest plane 8 is lowest)
                        end
                        
                        [isPart] = Core.MPParticleMovie.isPartPlane(currentCand,cand,direction,consThresh,zMethod);
                        if ~all(isPart ==0)
                            id = cand.plane(isPart);
                            particle(id,:) = cand(isPart,:);
                        end
                        
                    end
                   
                    %We remove the particle(s) from the list
                    focusMetric(ismember(candMet(:,1), particle(:,1))) = [];
                    candMet(ismember(candMet(:,1), particle(:,1)),:) = [];
                    %format particle to be the same as before:
                    [particle] = obj.makeParticle(particle);
                    
                     %Check if the resulting configuration of the plane make
                    %sense e.g. no hole in the configuration
                     planeConfig = particle.plane;
                    [checkRes] = Core.MPParticleMovie.checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod);
                    
                    %Store
                    if checkRes
                        
                        nPart = nPart +1;
                        %store particle in a new list
                        candidateList{nPart,1} = particle;
                        
                    else
                        %Otherwise we remove it from the best focus search list
                        %by putting focus metric to NaN
                        %focusMetric(idx) = NaN;
                        
                    end
                    
                    counter = counter+1;
                    
                end
            end     
    end

     methods (Static)
       
        function [isPart]   = isPartPlane(current, next, direction,consThresh,zMethod)
            %This function aim at determining whether a candidate from one
            %plane and the another are actually the same candidate on
            %different plane or different candidate. The decision is based
            %on threshold on localization distance, ellipticity and focus
            %metric.
            
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');
            
            thresh = consThresh;
            [checkRes1] = Core.MPParticleMovie.checkEuDist([current.row, current.col],...
                [next.row, next.col],thresh);
            
            if strcmp(zMethod,'PSFE')
             % Test ellipticity
                [checkRes2] = Core.MPParticleMovie.checkEllipticity(current.ellip,...
                next.ellip,direction);
            
            elseif or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                %we do not test ellipticity here
                checkRes2 = checkRes1;
                
            else
                
                error('Unknown Z method for consolidation');
                
            end
            
            % Test focus Metric
            maxExpFM = current.fMetric+0.1*current.fMetric;
            checkRes3 = next.fMetric < maxExpFM;
            
            %isPart will only be true for particle that passes the 3 tests
            isPart = checkRes1.*checkRes2.*checkRes3;
            
            if(length(find(isPart))>1)
                
                warning('Could not choose which particle was the partner of the requested particle, killed them both');
                isPart(isPart==1) = 0;
            end
            
            if isempty(isPart)
                isPart = false;
            end
            
            isPart = logical(isPart);
            
        end
        
        function [checkRes] = checkEuDist(current,next,Thresh)
            %Use to check if the Euclidian distance is within reasonable
            %range
            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;
            
            if isempty(checkRes)
                checkRes = false;
            end
            
        end
        
        function [checkRes] = checkEllipticity(current, next, direction)
            %Use to check if the Ellipticity make sense with what we
            %expect from the behavior of the PSFEngineering plate
            switch direction
                
                case 1
                    
                    ellip = current < next+0.1*next;
                    
                case -1
                    
                    ellip = current +0.1 *current > next;
                    
            end
            
            checkRes = ellip;

            if isempty(checkRes)
                checkRes = false;
            end
            
          
            
        end
        
        function [checkRes] = checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            if or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                if nPlanes ==1
                    testPlanes = true;
                else
                    testPlanes = sum(~isnan(planeConfig))>=2;
                end
            elseif strcmp(zMethod,'PSFE')
                switch nPlanes
                    case 1
                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;
                    case 2

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;

                    case 4

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 2;

                    case 8

                        nPlanesEdgeFrange = 1;
                        nPlanesEdgeInterleaved = 2;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 3;
                    otherwise
                        error('Unknown number of planes, only expect 1,2,4,8')
                end
                %Let us test that we have consolidate the particle in at least
                %3 Planes
                isEdgePlane = or(~isempty(find(planeConfig==1,1)),~isempty(find(planeConfig==8,1)));


                switch camConfig
                    case 'fullRange'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

                    case 'interleaved'
                        if isEdgePlane
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeInterleaved;
                        else

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesInterleaved;
                        end
                    case 'equal'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

                    otherwise
                        error('unknown camera configuration');
                end
            
            end
            if testPlanes
                %We check that there is no "Gap" in the plane configuration
                %as it would not make sense.
%                 testConsec = diff(planeConfig(~isnan(planeConfig)));
%                 checkRes = length(testConsec==1)>=2;
                checkRes = true;
            else
                
                checkRes = false;
                
            end
            
      

        end
        
        function [int,SNR]  = getIntensity(ROI,sig)
            %extract central position
            center = [ceil(size(ROI,1)/2),ceil(size(ROI,2)/2)];
            rowPos = center(1);
            colPos = center(2);
            %we integrate at 3*sigma (take ROI
            roiSignal = ceil(sig);
            if roiSignal(1) > center-1
                roiSignal = [center-1 center-1];
            end
            %get the idx for the ROI to integrate
            rowIdx = rowPos-roiSignal(1):rowPos+roiSignal(1);
            colIdx = colPos-roiSignal(2):colPos+roiSignal(2);
            
            %Pixel to integrate for signal
            px2SumInt = ROI(rowIdx,colIdx);
            
             %get background pixels
            bkg = ROI;
            bkg(rowIdx,colIdx) = 0;
            px2SumBkg = bkg(bkg~=0);
            bkg = mean(px2SumBkg);
            bkgVar = std(px2SumBkg);
            %calculate signal
            int = px2SumInt - bkg;
            int = sum(sum(int));
            %SNR = max(max(px2SumInt))/bkgVar;
            SNR = sqrt(int);
            if or(SNR<0,int<0)
                int = sum(sum(px2SumInt));
                SNR = sqrt(int);
            end
           
        end
       
     end
     
     methods (Access = protected)
        %method linked to candidate
        function [run,candidate] = existCandidate(obj,Path,ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'candidatePos')==true)
                        candidate = load([file2Analyze(1).folder filesep 'candidatePos.mat']);
                        candidate = candidate.candidate;
                        
                        if size(candidate,1)== obj.calibrated.nFrames
                            run = false;
                        else
                           disp('Detection missing in some frames, rerunning detection');
                           candidate = [];
                           run = true;
                        end
                            
                    else
                
                        run = true;
                        candidate =[];
                
                    end
                case 'run'
                    run = true;
                    candidate =[];
            end
        end
        
        %method linked to fitting
        function [run,SRLocPos] = existLocPos(obj,Path,ext) 
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'SRLocPos')==true)
                        SRLocPos = load([file2Analyze(1).folder filesep 'SRLocPos.mat']);
                        name = fieldnames(SRLocPos);
                        SRLocPos = SRLocPos.(name{1});
                        run = false;
                    else
                
                         run = true;
                        SRLocPos =[];
                
                    end
                case 'run'
                     run = true;
                     SRLocPos =[];
            end    
        end
        
        %method Linked to particles/planeConsolidation
        function [run, particle] = existParticles(obj,Path, ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        particle = particle.particle;
                        run = false;
                    else
                
                        run = true;
                        particle = [];
                
                    end
                    
                case 'run'
                    
                     run = true;
                     particle = [];
                     
            end
        end
        
        %Methods linked to Candidate
        function [candidate] = detectCandidate(obj,detectParam,frames,q)
            %Do the actual localization
            assert(~isempty(obj.calibrated),'Data should be calibrated to detect candidate');
            assert(isstruct(detectParam),'Detection parameter should be a struct with two fields');
            nFrames = length(frames);
            currentCandidate = obj.candidatePos;
            detectionMethod = obj.info.detectionMethod;

            if(isempty(currentCandidate))
                
                candidate = cell(obj.calibrated{1,q}.nFrames,1);
                
            else
                
                candidate = currentCandidate;
                
            end
           
            
            %parameter for localization
            FWHM_pix = obj.info.FWHM_px;
            delta  = detectParam.delta;
            chi2   = detectParam.chi2;
            h = waitbar(0,'detection of candidates...');
            
            if obj.info.rotational == 1
                for i = 1:obj.raw.maxFrame(1)
                    waitbar(i./obj.raw.maxFrame(1), h, 'Rotational: Averaging frames')
                    AllFrames(:,:,:,i) = obj.getFrame(i,q);
                    if ~isempty(obj.ROI)
                        for k = 1:size(obj.ROI(:,q),1)
                            ROI = obj.ROI{k,q};
                            for z = 1:size(obj.ROI{k,q}, 1)
                                ParticleMovie(:,:,ROI(z,5),i,k) = AllFrames(ROI(z,1):ROI(z,1)+ROI(z,3),...
                                    ROI(z,2):ROI(z,2)+ROI(z,4),ROI(z,5),i);
                            end
                        end
                    end
                end
                MeanIm = mean(AllFrames, 4);
                obj.calibrated{2,q} = MeanIm;

                if exist('ParticleMovie','var')
                    for i = 1:size(ParticleMovie, 5)
                        obj.ParticlesROI{i,q} = ParticleMovie(:,:,:,:,i);
                    end
                end
            end


            for i = 1 : 1:nFrames
                
                position = table(zeros(500,1),zeros(500,1),zeros(500,1),...
                    zeros(500,1),'VariableNames',{'row', 'col', 'meanFAR','plane'});
                [volIm] = obj.getFrame(frames(i),q);
                nPlanes = size(volIm,3);
                


                for j = 1:nPlanes
                    currentIM = volIm(:,:,j);
                    if obj.info.rotational == 1
                        currentIM = MeanIm(:,:,j);
                    end
                    %localization occurs here
                     switch detectionMethod 
                        case 'MaxLR'
                            [ pos, meanFAR, ~ ] = Localization.smDetection(currentIM,...
                                delta, FWHM_pix, chi2 );
                            if ~isempty(pos)
                                startIdx = find(position.row==0,1,'First');
                                if isempty(startIdx)
                                    startIdx = length(position.row)+1;
                                end
                                pos(:,3) = meanFAR;
                                pos(:,4) = j;
                                position(startIdx:startIdx+size(pos,1)-1,:) = array2table(pos);
                            else
                            end
                         case 'Intensity'
                             
                             bwImage = imbinarize(currentIM./max(currentIM(:)));
                             
                             SE = strel('disk',5);
                             bwImage = imopen(bwImage,SE);
                             bwImage = bwareaopen(bwImage,300);
    
                             [ctr] = regionprops(bwImage,'Area','Centroid');

                             pos = cat(1,ctr.Centroid);
                             if ~isempty(pos)
                                startIdx = find(position.row==0,1,'First');
                                if isempty(startIdx)
                                    startIdx = length(position.row)+1;
                                end
                                 pos = flip(pos,2);
                                 pos(:,3) = NaN;
                                 pos(:,4) = j;
                                 position(startIdx:startIdx+size(pos,1)-1,:) = array2table(pos);
                             end
                     end
                           
                end
                
                idx = find(position.row==0,1,'First');
                if isempty(idx)
                    if obj.info.rotational == 1
                        candidate{1} = position;
                        break
                    else
                        candidate{frames(i)} = position;
                    end                   
                else
                    if obj.info.rotational == 1
                        candidate{1} = position(1:idx-1,:);
                        break
                    else
                        candidate{frames(i)} = position(1:idx-1,:);
                    end
                    
                end
                waitbar(i/nFrames,h,...
                    sprintf('detection of candidates in Frame %d/%d done',i,nFrames));
            end
            
            close(h);
        end
                
        function [candMet] = superResLocFit(obj,data,frameCandidate,roiSize)
            %Candidate metric are determined here (x,y,e,+focusmetric)
            delta = roiSize;
            
            %initialize table
            varNames = {'row','col','z','ellip','magX','magY','meanFAR','fMetric','gFitMet','plane'};
                candMet = table(zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    'VariableNames',varNames);
                sigSetup = [obj.info.sigma_px obj.info.sigma_px];
                
            for i = 1:size(frameCandidate,1)
                
                plane = frameCandidate.plane(i);
                planeData = double(data(:,:,plane));
                %Get the ROI
                [roi_lims] = EmitterSim.getROI(frameCandidate.col(i), frameCandidate.row(i),...
                    delta, size(planeData,2), size(planeData,1));
                ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                roiSize = size(ROI);
                if and(roiSize(1)==roiSize(2),mod(roiSize(1),2)==1)
                    if strcmpi(obj.info.fitMethod,'phasor')
                        %Phasor fitting to get x,y,e
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;

                    elseif strcmpi(obj.info.fitMethod,'Gauss')
                        [X,Y] = meshgrid(frameCandidate.col(i)-delta:frameCandidate.col(i)+...
                            delta,frameCandidate.row(i)-delta:frameCandidate.row(i)+delta);
                        domain(:,:,1) = X;
                        domain(:,:,2) = Y;

                        %Gauss (slower)
                        [gPar] = Localization.Gauss.MultipleFitting(ROI,frameCandidate.col(i),...
                            frameCandidate.row(i),domain,1);%data,x0,y0,domain,nbOfFit
                        colPos = gPar(5); %Should be directly the position of the particle as we
                            %gave above the domain of the ROI in the space of the image
                        rowPos = gPar(6);

                        row = rowPos - round(frameCandidate.row(i));
                        col = colPos - round(frameCandidate.col(i));

                        e = gPar(3)/gPar(2);

                        magX = 0;
                        magY = 0;
                    else
                        %Phasor fitting to get x,y,e
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;
                    end

                     %LRT focus metric
                    [fMetric,~] = Localization.likelihoodRatioTest(ROI,sigSetup,[row col]);

                    if magX>=magY
                        sig(1) = sigSetup(1) * magX/magY;
                        sig(2) = sigSetup(2);
                    else
                        sig(1) = sigSetup(1);
                        sig(2) = sigSetup(2) * magY/magX;
                    end

                    [int,SNR] = obj.getIntensity(ROI,sigSetup);
                    %LRT focus metric
                    [gFitMet,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);

                    %storing info
                    candMet.row(i) = rowPos;
                    candMet.col(i) = colPos;
                    candMet.z(i) = 0;
                    candMet.ellip(i) = e;
                    candMet.magX(i) = magX;
                    candMet.magY(i) = magY;
                    candMet.intensity(i) = int;
                    candMet.SNR(i) = SNR;
                    candMet.meanFAR(i) = frameCandidate.meanFAR(i);
                    candMet.fMetric(i) = fMetric;
                    candMet.gFitMet(i) = gFitMet;
                    candMet.plane(i) = plane;
                else
                    
                    %storing info
                    candMet.row(i)          = NaN;
                    candMet.col(i)          = NaN;
                    candMet.z(i)            = NaN;
                    candMet.ellip(i)        = NaN;
                    candMet.magX(i)         = NaN;
                    candMet.magY(i)         = NaN;
                    candMet.intensity(i)    = NaN;
                    candMet.SNR(i)          = NaN;
                    candMet.meanFAR(i)      = NaN;
                    candMet.fMetric(i)      = NaN;
                    candMet.gFitMet(i)      = NaN;
                    candMet.plane(i)        = NaN;
                end
            end
            %remove NaN's  
            idx = isnan(candMet.row);
            candMet(idx,:) = [];

        end
        
        function [doAvg]  = checkDoAverage(obj,ellip)
            camConfig = obj.calibrated{1,1}.camConfig;
            switch camConfig
                case 'fullRange'

                    if and(ellip>0.8,ellip<1.25)

                        doAvg = false;
                    else
                        doAvg = true;
                    end

                case 'interleaved'

                        doAvg = true;

                case 'equal'

                    if and(ellip>0.8,ellip<1.25)

                        doAvg = false;
                    
                    else
                        
                        doAvg = true;
                    
                    end

                otherwise
                    error('unknown camera config');
            end
        end
        
        function [newPart] = makeParticle(~,particleData)
            newPart = array2table(nan(5,size(particleData,2)));
            newPart.Properties.VariableNames = particleData.Properties.VariableNames;
            %store best focus in center
            [~,idx] = nanmax(particleData.fMetric);
            if idx-2 > 0
                newPart(1,:) = particleData(idx-2,:);
            end
            
            if idx-1 > 0
                newPart(2,:) = particleData(idx-1,:);
            end
            
            newPart(3,:) = particleData(idx,:);
            
            if idx+1 <= 8 && idx+1<= height(particleData)
                newPart(4,:) = particleData(idx+1,:);
            end
            
            if idx+2 <= 8 && idx+2 <= height(particleData)
                newPart(5,:) = particleData(idx+2,:);
            end
        
        end
       
     end
end

