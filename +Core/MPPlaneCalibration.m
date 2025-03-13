classdef MPPlaneCalibration < handle
    %PlaneCalibration is a class that will hold the information of several
    %MPCalMovie and be able to process these information to create a
    %calibration file from the data of all the movies
    
    properties (SetAccess = 'private')
        path
        ext
        MPCalibrations
        info
        allCal
        cal        
    end
    
    methods
        
        function obj = MPPlaneCalibration(path2MPCal,info)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
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
            assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            
            obj.path = path;
            
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
                    tmp = Core.MPCalibration(file, obj.info);
                    obj.MPCalibrations.(['MPCal' num2str(i-2)]) = tmp;

                else

                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);

                end

            end
            disp('=======> DONE ! <========')
        end
        
        function calcIndivCal(obj)
            
            fieldN = fieldnames(obj.MPCalibrations);
            nfields = numel(fieldN);
      
            for i = 1: nfields
                
                if isfield(obj.info,'nChan')
                    nChan = obj.info.nChan;
                else
                    nChan = 4;
                end
                
                disp(['Retrieving data from MPCalibration file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                
                obj.MPCalibrations.(fieldN{i}).calc(nChan);
                currentCal = obj.MPCalibrations.(fieldN{i}).getCal;
                
                obj.allCal(i).file = currentCal.file;
                
            end
        end
        
        function calcCombinedCal(obj)
            disp('Combining data from different calibration...');
            allData = obj.allCal;
            nFiles = length(allData);
            nPlanes1 = length(allData(1).file.neworder1);
            allROI1 = zeros([size(allData(1).file.ROI1) nFiles]);
            allFocusMet1 =  zeros([size(allData(1).file.focusMet1) nFiles]);
            allFit1 = zeros([size(allData(1).file.fit1) nFiles]);
            allNewOrder1 =  zeros([size(allData(1).file.neworder1) nFiles]);
            allICorrF1   =  zeros([size(allData(1).file.neworder1) nFiles]);
            inFocus1 = allData(1).file.inFocus1;
            inFocus1 = rmfield(inFocus1,{'frame', 'zpos'});
            allRelZpos1 = zeros(nPlanes1,nFiles);
            allZpos1 = zeros(nPlanes1,nFiles);
  

            if allData(1).file.multiModal == true
                nPlanes2 = length(allData(1).file.neworder2);
                allROI2 = zeros([size(allData(1).file.ROI2) nFiles]);
                allROI2FullCam = zeros([size(allData(1).file.ROI2FullCam) nFiles]);
                allFocusMet2 =  zeros([size(allData(1).file.focusMet2) nFiles]);
                allFit2 = zeros([size(allData(1).file.fit2) nFiles]);
                allNewOrder2 =  zeros([size(allData(1).file.neworder2) nFiles]);
                allICorrF2   =  zeros([size(allData(1).file.neworder2) nFiles]);
                inFocus2 = allData(1).file.inFocus2;
                inFocus2 = rmfield(inFocus2,{'frame', 'zpos'});
                allRelZpos2 = zeros(nPlanes2,nFiles);
                allZpos2 = zeros(nPlanes2,nFiles);
            else 
            end

            for i = 1:nFiles
                allROI1(:,:,i) = allData(i).file.ROI1;
                %allFocusMet(:,:,i) = allData(i).file.focusMet;
                %allFit(:,:,i) = allData(i).file.fit;
                %allNewOrder(:,:,i) = allData(i).file.neworder; 
                allICorrF1(:,:,i)  = allData(i).file.Icorrf1;   
                tmp1 = cell2mat({allData(i).file.inFocus1.zpos});
                allRelZpos1(:,i) = tmp1-mean(tmp1);
                allZpos1(:,i) = tmp1;

                if allData(1).file.multiModal == true
                    allROI2(:,:,i) = allData(i).file.ROI2;
                    allROI2FullCam(:,:,i) = allData(i).file.ROI2FullCam;
                    %allFocusMet(:,:,i) = allData(i).file.focusMet;
                    %allFit(:,:,i) = allData(i).file.fit;
                    %allNewOrder(:,:,i) = allData(i).file.neworder; 
                    allICorrF2(:,:,i)  = allData(i).file.Icorrf2;   
                    tmp2 = cell2mat({allData(i).file.inFocus2.zpos});
                    allRelZpos2(:,i) = tmp2-mean(tmp2);
                    allZpos2(:,i) = tmp2;
                end
            end

            ROI1 = floor(nanmean(allROI1,3));
            RelZPos1 = nanmean(allRelZpos1,2);
            tmp1 = num2cell(RelZPos1');
            [inFocus1(1,:).relZPos] = tmp1{:};
            tmp1 = num2cell(nanmean(allZpos1,2)');
            [inFocus1(1,:).zPos] = tmp1{:};
            obj.cal.nFiles = nFiles;
            obj.cal.fullPath = obj.path;
            obj.cal.file = obj.allCal(1).file;
            obj.cal.file = rmfield(obj.cal.file,{'focusMet1','fit1','Zpos'});
            obj.cal.file.ROI1 = ROI1;
            obj.cal.file.Icorrf1 = squeeze(nanmean(allICorrF1,3));
            obj.cal.file.inFocus1 = inFocus1;

            if allData(1).file.multiModal == true
                ROI2 = floor(nanmean(allROI2,3));
                ROI2FullCam = floor(nanmean(allROI2FullCam,3));
                RelZPos2 = nanmean(allRelZpos2,2);
                tmp2 = num2cell(RelZPos2');
                [inFocus2(1,:).relZPos] = tmp2{:};
                tmp2 = num2cell(nanmean(allZpos2,2)');
                [inFocus2(1,:).zPos] = tmp2{:};
                obj.cal.nFiles = nFiles;
                obj.cal.fullPath = obj.path;
                %obj.cal.file = obj.allCal(1).file;
                %obj.cal.file = rmfield(obj.cal.file,{'focusMet2','fit2','Zpos'});
                ROI2(:,3) = ones(size(ROI2, 1), 1)*min(ROI2(1,3), ROI1(1,3));
                ROI2(:,4) = ones(size(ROI2, 1), 1)*min(ROI2(1,4), ROI1(1,4));
                ROI1(:,3) = ones(size(ROI2, 1), 1)*min(ROI2(1,3), ROI1(1,3));
                ROI1(:,4) = ones(size(ROI2, 1), 1)*min(ROI2(1,4), ROI1(1,4));
                obj.cal.file.ROI2 = ROI2;
                obj.cal.file.ROI2FullCam = ROI2FullCam;
                obj.cal.file.Icorrf2 = squeeze(mean(allICorrF1,3));
                obj.cal.file.inFocus2 = inFocus2; 


                % for j = 1:length(allTransformations)
                %    obj.cal.file.Transformation{j,1}.Dimensionality = nanmean(allTransformations{j,1}.Dimensionality, 3);
                %    obj.cal.file.Transformation{j,1}.Scale = nanmean(allTransformations{j,1}.Scale, 3);
                %    obj.cal.file.Transformation{j,1}.RotationAngle = nanmean(allTransformations{j,1}.RotationAngle, 3);
                %    obj.cal.file.Transformation{j,1}.Translation = nanmean(allTransformations{j,1}.Translation, 3);
                %    obj.cal.file.Transformation{j,1}.R = nanmean(allTransformations{j,1}.R, 3);
                %    obj.cal.file.Transformation{j,1}.A = nanmean(allTransformations{j,1}.A, 3);
                %    %obj.cal.file.Transformation{j,1} = simtform2d(mean(allTransformations{j,1}.Scale, 3), mean(allTransformations{j,1}.RotationAngle, 3), mean(allTransformations{j,1}.Translation, 3))
                % end
            else
            end
            
            % if allData(1).file.multiModal == true
            %     [Transformations, ZDiff] = mpSetup.cali.getTransformationMultiModal(obj);
            %     obj.cal.file.Transformation = Transformations;
            %     avZDiff = nanmean(ZDiff,1);
            %     tmp2 = cell2mat({obj.cal.file.inFocus2.zPos}) - avZDiff;
            %     RelPos = tmp2 - mean(tmp2);
            %     for m = 1:size(obj.cal.file.inFocus2, 2)
            %         obj.cal.file.inFocus2(m).relZPos = RelPos(m);
            %         obj.cal.file.inFocus2(m).zPos = tmp2(m);
            %     end
            % end

            [~] = obj.determineCAMConfig;
            disp('================>DONE<====================');
            
            
            
        end
        
        function save(obj)
            cal2D = obj.getCal;
            filePath = obj.path;
            fileName = [filePath filesep '2DCal.mat'];
            
            save(fileName,'cal2D')
            disp('The calibration was succesfully saved');
        end
        
        function showCal(obj,idx)
            fields = fieldnames(obj.MPCalibrations);
            mov2Use = obj.MPCalibrations.(fields{idx});
            focusMet1 = mov2Use.cal.file.focusMet1;
            fit1 = mov2Use.cal.file.fit1(:,2:2:end);
            fitZ1 = mov2Use.cal.file.fit1(:,1:2:end);
            ZPos = mov2Use.cal.file.Zpos;
            color = hsv(8);
            FocusZ1 = {mov2Use.cal.file.inFocus1.zpos};

            figure()
            subplot(double(mov2Use.cal.file.multiModal)+1,1,1)
            hold on
            leg1 = cell(1,size(focusMet1,2));
            height1 =  max(max(focusMet1));
            y = 1:height1;
            for i = 1 : size(focusMet1,2)
                [~,idx] = max(fit1(:,i));
                scatter(ZPos-mean(ZPos),focusMet1(:,i),[],color(i,:),'filled')
                plot(fitZ1(:,i)-mean(ZPos),fit1(:,i),'Color', color(i,:),'LineWidth',2.5,'HandleVisibility','off')
               
                x = ones(1,length(y))*(FocusZ1{i}-mean(ZPos));
                plot(x(:),y(:),'k--','HandleVisibility','off');
                
                leg1{i} = ['Cam' num2str(obj.cal.file.inFocus1(i).cam) ' - Plane' num2str(obj.cal.file.inFocus1(i).ch)];
                
            end
            ylim([min(min(focusMet1)), max(max(focusMet1))]);
            xlim([round(min(ZPos-mean(ZPos))), round(max(ZPos-mean(ZPos)))]);
            title('Setup Plane Calibration');
            ylabel('Intensity (a.u.)')
            xlabel('z position (um)')
            legend(leg1)
            
            if mov2Use.cal.file.multiModal == true
                hold on
                focusMet2 = mov2Use.cal.file.focusMet2;
                fit2 = mov2Use.cal.file.fit2(:,2:2:end);
                fitZ2 = mov2Use.cal.file.fit2(:,1:2:end);
                FocusZ2 = {mov2Use.cal.file.inFocus2.zpos};

                subplot(double(mov2Use.cal.file.multiModal)+1,1,2)
                hold on
                leg2 = cell(1,size(focusMet2,2));
                height2 =  max(max(focusMet2));
                y = 1:height2;
                for i = 1 : size(focusMet2,2)
                    [~,idx] = max(fit2(:,i));
                    scatter(ZPos-mean(ZPos),focusMet2(:,i),[],color(i,:),'filled')
                    plot(fitZ2(:,i)-mean(ZPos),fit2(:,i),'Color', color(i,:),'LineWidth',2.5,'HandleVisibility','off')
                    x = ones(1,length(y))*(FocusZ2{i}-mean(ZPos));
                    plot(x(:),y(:),'k--','HandleVisibility','off');
                    
                    leg2{i} = ['Cam' num2str(obj.cal.file.inFocus2(i).cam) ' - Plane' num2str(obj.cal.file.inFocus2(i).ch)];                  
               end
               ylim([min(min(focusMet2)), max(max(focusMet2))]);
               xlim([round(min(ZPos-mean(ZPos))), round(max(ZPos-mean(ZPos)))]);
               title('Setup Plane Calibration - planes 9-16');
               ylabel('Intensity (a.u.)')
               xlabel('z position (um)')
               legend(leg2)
            else
            end
        end
        
        function [camConfig] = determineCAMConfig(obj)
        multiModal = obj.allCal(1).file.multiModal;
        if multiModal == true
            j = 2
        else 
            j = 1
        end

        for i = 1:j
            if i == 1
                inFocus = obj.cal.file.inFocus1;
            elseif i == 2
                inFocus = obj.cal.file.inFocus2;
            end
            relZPos = cell2mat({inFocus.relZPos});
            cams    = cell2mat({inFocus.cam});
            
            if i == 2
                meanCam1 = mean(relZPos(cams==3));
                meanCam2 = mean(relZPos(cams==4));
            else
                meanCam1 = mean(relZPos(cams==1));
                meanCam2 = mean(relZPos(cams==2));
            end
            
            dCam     = abs(meanCam1-meanCam2)*1000;
            
            relZPos = relZPos(obj.cal.file.neworder1);
            minDist = min(abs(diff(relZPos(2:end))))*1000;
            %If distance between camera center is above 500 nm, most likely is
            %full cam
            if dCam > 500
                camConfig = 'fullRange';
            %if minDist is relatively large and dCam is small it is probably
            %interleaved
            elseif and(minDist>80, dCam<300)
                camConfig = 'interleaved';
            %if minDist is small then probably equal
            elseif minDist<80
                camConfig = 'equal';
            else
                error('Something is wrong with your distance between planes')
            end
            
            obj.cal.camConfig{i} = camConfig;
        end

        if multiModal == true 
            assert(strcmpi(obj.cal.camConfig(1), obj.cal.camConfig(2)), 'Different camera configuration for lower and upper planes. Distance between planes 1-8 is not corresponding to distance between planes 9-16')
        else
        end
        end
        
        function offTarget(obj)
            newOrder1 = obj.cal.file.neworder1;
            zPos1 = abs([obj.cal.file.inFocus1.zPos]);
            zPos1 = zPos1(newOrder1);
            switch obj.cal.camConfig{1,1}
                case 'fullRange'
                    distBetweenCamPlanes1 = abs(mean(diff(zPos1(1:end/2))) + mean(diff(zPos1(end/2+1:end))))/2;
                    target1    = distBetweenCamPlanes1;
                    distBetweenCam1 = abs(zPos1(5)-zPos1(4));
                    offTarget1 = distBetweenCam1 - target1;
                   
                case 'interleaved'
                    distBetweenCamPlanes1 = abs(mean(diff(zPos1(1:2:end))) + mean(diff(zPos1(2:2:end))))/2;
                    target1    = distBetweenCamPlanes1/2;
                    distBetweenPlane1 = abs(diff(zPos1));
                    offTarget1 = distBetweenPlane1 - target1;
                    offTarget1 = mean(abs(offTarget1));
  
                case 'equal'
                    
                    distBetweenPlanes1 = abs(diff(zPos1));
                    distBetweenPlanes1 = distBetweenPlanes1(1:2:end);
                    offTarget1 = mean(distBetweenPlanes1);
            end
           fprintf('The difference between the target and the current plane conformation (planes 1-8) \nis %d nm',round(offTarget1*1000));
            if obj.cal.file.multiModal == true
                newOrder2 = obj.cal.file.neworder2;
                zPos2 = abs([obj.cal.file.inFocus2.zPos]);
                zPos2 = zPos1(newOrder2);
                switch obj.cal.camConfig{1,2}
                    case 'fullRange'
                        distBetweenCamPlanes2 = abs(mean(diff(zPos2(1:end/2))) + mean(diff(zPos2(end/2+1:end))))/2;
                        target2    = distBetweenCamPlanes2;
                        distBetweenCam2 = abs(zPos2(5)-zPos2(4));
                        offTarget2 = distBetweenCam1 - target1;
                       
                    case 'interleaved'
                        distBetweenCamPlanes2 = abs(mean(diff(zPos2(1:2:end))) + mean(diff(zPos2(2:2:end))))/2;
                        target2    = distBetweenCamPlanes2/2;
                        distBetweenPlane2 = abs(diff(zPos2));
                        offTarget2 = distBetweenPlane1 - target2;
                        offTarget2 = mean(abs(offTarget2));
      
                    case 'equal'
                        
                        distBetweenPlanes2 = abs(diff(zPos2));
                        distBetweenPlanes2 = distBetweenPlanes2(1:2:end);
                        offTarget2 = mean(distBetweenPlanes2);
                end
                fprintf('The difference between the target and the current plane conformation (planes 9-16) \nis %d nm',round(offTarget2*1000));

            else
            end
        end
        
        
        
     end
end