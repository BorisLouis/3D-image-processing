classdef MPMovie < Core.Movie
    %From this class we assume that we are dealing with multiplane data.
    %The movie can thus be calibrated if calibration data is provided.
    %Extension will be centered around localization but SOFI movie could
    %also inherit from this object
    
    properties (SetAccess = 'protected')
        cal2D %calibration data
        calibrated %Path to calibrated data and Info
    end
    
    methods
        
        function obj = MPMovie(raw,cal,info)
            if ~isfield(info,'multiModal')
                info.multiModal = 0;
            end

            
            obj = obj@Core.Movie(raw,info);
            
            obj.cal2D = cal;
            
        end
        
        function set.cal2D(obj,cal2D)
            if isempty(cal2D)
                obj.cal2D;
                disp('No 2D calibration received, assuming it is not multiplane Data');
            else
                if ischar(cal2D)
                    assert(isfolder(cal2D), 'The path given is not a folder, 2DCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
                    
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
                    
                else
                    assert(and(isstruct(cal2D), and(isfield(cal2D,'camConfig'),isfield(cal2D,'file'))),...
                        '2D calibration is supposed to be a struct with 4 fields');
                    obj.cal2D = cal2D;
                    
                end
            end
            
        end
        
        function set.calibrated(obj,calibrated)
            
            if isstruct(calibrated)
                test = calibrated;
                calibrated = [];
                calibrated{1,1} = test;
            end
            assert(iscell(calibrated),'Calibrated is expected to be a struct');
            obj.calibrated = calibrated;
            
        end
        
        function calibrate(obj)
            %Method that the user will call to calibrate the data
            folderContent = dir(obj.raw.movInfo.Path);
            idx2Calibrated = contains({folderContent.name}, 'calibrated');
            
            %if there is only 1 diff value, this value must be 0 and thus
            %calibrated folder does not exist thus, we calibrate
            if length(unique(idx2Calibrated))<2
                
                disp('Calibrating the dataset');
                [calibrate] = obj.applyCalib;
                disp('Data is now calibrated');
                
                %if there is 2 value, then calibrated folder exist and then we
                %check if a calibration file is in there.
            elseif length(unique(idx2Calibrated))==2
                for q = 1:obj.info.multiModal+1
                    name = append('calibrated', num2str(q))

                    fullPath = [obj.raw.movInfo.Path filesep name];
                    [file2Analyze] = obj.getFileInPath(fullPath, '.mat');
                    
                    if (~isempty(file2Analyze))
                        if obj.info.calibrate == true
                            %need to delete existing tif files
                            folder = file2Analyze.folder;
                            filePattern = fullfile(folder,'*.tif');
                            theFiles = dir(filePattern);
                            fprintf(1,'Deleting calibration in folder: %s \n',folder);
                            
                            for i = 1: length(theFiles)
                                name = theFiles(i).name;
                                fullFName = fullfile(folder,name);
                                fprintf(1,'Deleting %s \n',name);
                                delete(fullFName);
                            end
                            
                            disp('Calibrating the dataset');
                            [calibrate] = obj.applyCalib;
                            disp('Data is now calibrated');
                        else
                            
                            
                            [file] = obj.getFileInPath (fullPath,'.tif');
                            
                            if ~isempty(file) %only work with 8 Planes now
                                %If data is already calibrated we load it
                                disp('The dataset is already calibrated, Loading from existing file');
                                
                                if length(file) ~= 8
                                    
                                    warning('Did not find 8 planes, if you are not using the prism it is okay, otherwise you might want to recalibrate');
                                    
                                end
                            end
                            
                            % find the calibrated.mat file in the folder
                            for i = 1:length(file2Analyze)
                                k = strfind(convertCharsToStrings(file2Analyze(i).name), 'calibrated'); 
                                if k == 1
                                    idx = i;
                                end
                            end
                            fullpath = [file2Analyze(idx).folder filesep file2Analyze(idx).name];
                            tmp = load(fullpath);
                            calibrate{q} = tmp.calib;
                            fieldN = fieldnames(calibrate{q}.filePath);
                            for i = 1: length(fieldN)
                                
                                currentPath = calibrate{q}.filePath.(fieldN{i});
                                idx1 = strfind(fullPath,'\calibrated');
                                idx2 = strfind(currentPath,'\calibrated');
                                newPath = [fullPath(1:idx1-1) currentPath(idx2(1):end)];
                                calibrate{q}.filePath.(fieldN{i}) = newPath;
                                
                            end
                            
                            if isfield(calibrate{q},'transPath')
                                fieldN = fieldnames(calibrate{q}.transPath);
                                
                                for i = 1: length(fieldN)
                                    
                                    currentPath = calibrate{q}.transPath.(fieldN{i});
                                    idx1 = strfind(fullPath,'\calibrated');
                                    idx2 = strfind(currentPath,'\calibrated');
                                    newPath = [fullPath(1:idx1-1) currentPath(idx2(1):end)];
                                    calibrate{q}.transPath.(fieldN{i}) = newPath;
                                    
                                end
                            end
                            
                            disp('Done');
                        end
                        
                    else
                        
                        %If no calibration was found we calibrate again
                        disp('Calibrating the dataset');
                        [calibrate] = obj.applyCalib;
                        disp('Data is now calibrated');
                        
                    end
                end
            else
                
                error('Something is wrong with your calibrated directory');
                    
            end
            
            obj.calibrated = calibrate;
            
        end
        
        function h = showFrame(obj,idx,scaleBar,idx2Plane)
            switch nargin
                case 3
                    idx2Plane = [];
                case 4
                otherwise
                    error('Not enough input argument');
            end
            IScale = 1;
            %Adapted method from the Core.Movie one, its behavior changed
            %depending on whether the data has been calibrated or not.
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            %check the frame requested
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %Get the data of the requested frame
            [frame] = getFrame(obj,idx);
            pxSize = obj.info.pxSize/1000;%in µm
            scaleBarPx = scaleBar/pxSize;
            planes = size(frame,3);
            
            if ~isempty(idx2Plane)
                for i = 1 :length(planes)
                    if i ~=idx2Plane
                        frame = rmfield(frame,planes{i});
                    end
                    
                end
                nsFig = 1;
                
            else
                nsFig = ceil(planes/4);
            end
            
            if ~isempty(obj.calibrated)
                
                zPos = obj.calibrated.oRelZPos;
                
            else
                
                zPos = zeros(planes,1);
                
            end
            
            %Displaying occur hear
            h = figure(1);
            h.Name = sprintf('Frame %d',idx);
            
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            
            for i = 1:planes
                currentIM = frame(:,:,i);
                subplot(nsFig,planes/nsFig,i)
                
                if strcmp(obj.info.type,'transmission')
                    colormap('gray');
                    currentIM = imcomplement(currentIM);
                else
                    colormap('jet')
                end
                
                imagesc(currentIM)
                ax = gca;
                caxis(ax,[min(min(min(frame))), IScale*max(max(max(frame)))]);
                
                hold on
                x = size(currentIM,2)-scaleBarPx-(0.05*size(currentIM,2)):size(currentIM,2)-0.05*size(currentIM,2);
                y = ones(1,length(x))*size(currentIM,1)-0.05*size(currentIM,2);
                text(mean(x),mean(y)-0.04*size(currentIM,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',8);
                plot(x,y,'-w','LineWidth',2.5);
                axis image;
                
                %grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({['plane ',num2str(i)],sprintf(' Zpos = %0.3f',zPos(i))});
                hold off
                
                
            end
        end
        
        function [data, CorrFactor] = getFrame(obj,idx, q)
            %Allow the user to extract data from a specific frame, behavior
            %depends on the calibration
            assert(length(idx)==1,'Only one frame at a time');
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %Behavior depend on status
            CorrFactor = [];
            if isempty(obj.calibrated)
                
                [data] = getFrame@Core.Movie(obj,idx);
                
            elseif isstruct(obj.calibrated{1,q})
                fieldsN = fieldnames(obj.calibrated{1,q}.filePath);                
                %data = zeros(obj.calibrated{1,q}.Height,obj.calibrated{1,q}.Width,numel(fieldsN));
                for i = 1:numel(fieldsN)
                    %Load plane
                    [mov] = Load.Movie.tif.getframes(obj.calibrated{1,q}.filePath.(fieldsN{i}),idx);
                    bg = double(imgaussfilt(mov, 15));
                    mov = double(mov) - bg;
                    data(:,:,i) = mov;
                    if q == 2
                        [mov1] = Load.Movie.tif.getframes(obj.calibrated{1,1}.filePath.(fieldsN{i}),idx);
                        bg1 = double(imgaussfilt(mov1, 15));
                        Corr = mean(bg1./bg, 'all');
                        mov = mov.*Corr;
                    end
                end
            end
        end
        
        function [data] = getPlane(obj,idx,frame)
            %Allow the user to extract data from a specific plane, behavior
            %depends on the calibration
            switch nargin
                case 1
                    frame = 1:obj.raw.maxFrame(1);
                case 2
            end
            assert(and(idx<9,idx>=1),'plane should be between 1 and 8');
            
            [data] = Load.Movie.tif.getframes(obj.calibrated.filePath.(sprintf('plane%d',idx)),frame);
            
            
        end
        
        function [calibrated] = getCalibrated(obj)
            
            calibrated = obj.calibrated;
            
        end
        
        function playMovie(obj)
            %TODO: code this function
        end
        
        function [xStep, xPosMotor] = getXPosMotor(obj)
            if strcmpi(obj.raw.ext,'.ome.tif')
                %small function that extract the zStep from the info in the raw
                nFrames = obj.raw.maxFrame(1);
                xPosMotor = zeros(nFrames,1);
                
                for i = 1 : obj.raw.maxFrame(1)
                    
                    xPosMotor(i) = obj.raw.frameInfo(2*i).Pos(1);
                    
                end
                
                xStep = diff(xPosMotor);
            else
                xStep = zeros(obj.raw.maxFrame(1),1);
                xPosMotor = zeros(obj.raw.maxFrame(1),1);
            end
        end
        
        function [yStep, yPosMotor] = getYPosMotor(obj)
            if strcmpi(obj.raw.ext,'.ome.tif')
                %small function that extract the zStep from the info in the raw
                nFrames = obj.raw.maxFrame(1);
                yPosMotor = zeros(nFrames,1);
                
                for i = 1 : obj.raw.maxFrame(1)
                    
                    yPosMotor(i) = obj.raw.frameInfo(2*i).Pos(2);
                    
                end
                
                yStep = diff(yPosMotor);
            else
                yStep = zeros(obj.raw.maxFrame(1),1);
                yPosMotor = zeros(obj.raw.maxFrame(1),1);
            end
            
        end
        
        function [zStep, zPosMotor] = getZPosMotor(obj)
            if strcmpi(obj.raw.ext,'.ome.tif')
                %small function that extract the zStep from the info in the raw
                nFrames = obj.raw.maxFrame(1);
                zPosMotor = zeros(nFrames,1);
                
                for i = 1 : obj.raw.maxFrame(1)
                    
                    zPosMotor(i) = obj.raw.frameInfo(2*i).Pos(3);
                    
                end
                
                zStep = diff(zPosMotor);
            else
                zStep = zeros(obj.raw.maxFrame(1),1);
                zPosMotor = zeros(obj.raw.maxFrame(1),1);
            end
            
        end
        
    end
    
    methods (Static)
        
    end
    
    methods (Access = private)
        
        function [calib] = applyCalib(obj)
            
            frameInfo = obj.raw.frameInfo;
            movInfo   = obj.raw.movInfo;
            maxFrame = obj.raw.movInfo.maxFrame(1);
            step = 800;
            frame2Load = obj.info.frame2Load;
            if ischar(frame2Load)
                frame2Load = 1:maxFrame;
            else
                frame2Load = Core.Movie.checkFrame(frame2Load,maxFrame);
            end
            
            endFrame = max(frame2Load);
            startFrame = frame2Load(1);
            
            nStep = ceil((endFrame-startFrame)/step);
            frame = startFrame:step:endFrame;
            disp('Calibrating the dataset');
            for i = 1:nStep
                disp(['Calibrating step ' num2str(i)]);
                if i < nStep
                    
                    cFrame = frame(i):frame(i+1)-1;
                else
                    
                    cFrame = frame(i):endFrame;
                    
                end
                %change behavior depending on extension:
                switch obj.raw.movInfo.ext
                    case '.ome.tif'
                        assert(~isempty(obj.cal2D.file),'no calibration file provided, cannot calibrate');
                        MP = true;
                        % load the raw data
                        [ movC1, movC2] = Load.Movie.ome.load( frameInfo, movInfo, cFrame );

                        %applying the calibration
                        if exist('ROI','var') == 1
                            [data,isTransmission,ROI, Bg] = mpSetup.cali.apply( movC1, movC2, obj.cal2D.file, obj.info, ROI);
                        else 
                            [data,isTransmission,~, Bg] = mpSetup.cali.apply( movC1, movC2, obj.cal2D.file, obj.info);
                        end
                        %saving data per plane and info to cal
                        [calib] = obj.saveCalibrated(data,endFrame,isTransmission,MP, Bg);
                        
                    otherwise
                        extName = strrep(obj.raw.movInfo.ext,'.','');
                        
                        MP = false;%not a multiplane data
                        [movC1] = Load.Movie.(extName).getFrame(obj.raw.fullPath,cFrame);
                        %reformat in X,Y,P,T format for the saveCalibrated
                        %part
                        if extName == 'mpg' %#ok<BDSCA>
                            data(:,:,1,1:size(movC1,3)) = uint32(movC1);
                        else
                            data(:,:,1,1:size(movC1,3)) = uint16(movC1);
                        end
                        isTransmission = false;
                        [calib] = obj.saveCalibrated(data,endFrame,isTransmission,MP);
                       
                end
            end
            
        end
        
        function [calib] = saveCalibrated(obj,data,maxFrame,isTransmission,MP, Bg)
            
            calDir1 = [obj.raw.movInfo.Path filesep 'calibrated1'];
            calTransDir1 = [calDir1 filesep 'Transmission'];
            mkdir(calDir1);
            mkdir(calTransDir1);
            %Save the resulting planes in separated TIF and save a txt info
            %file
            if MP
                cal = obj.cal2D.file;
            end
            fid = fopen([calDir1 filesep 'CalibratedInfo.txt'],'w');
            fprintf(fid,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
            calib1.mainPath = calDir1;
            if islogical(isTransmission)
                Test = isTransmission;
                isTransmission = [];
                isTransmission{1,1} = Test;
            end

            calib1.nPlanes = sum(~isTransmission{1,1});

            idx2Plane = 1;

            if isstruct(data)
            elseif iscell(data)
            else
                test = data;
                data = [];
                data{1,1} = test;
            end

            
            for i = 1:size(data{1,1},3)
                data2Store1 = squeeze(data{1,1}(:,:,i,:));

                isTrans1 = isTransmission{1,1}(i);
                if strcmpi(obj.info.type,'transmission')
                    data2Store1 = imcomplement(data2Store1);
                end
                fieldN = sprintf('plane%d',i);
                fName = sprintf('calibratedPlane%d.tif',i);
                if isTrans1
                    fPathTiff = [calTransDir1 filesep fName];
                    calib1.transPath.(fieldN) = fPathTiff;
                else
                    fPathTiff = [calDir1 filesep fName];
                    calib1.filePath.(fieldN) = fPathTiff;
                end
                
                calib1.nFrames = maxFrame;
                t1 = Tiff(fPathTiff, 'a');
                
                if isa(data2Store1,'uint32')
                    t1 = dataStorage.writeTiff(t1,data2Store1,32);
                else
                    t1 = dataStorage.writeTiff(t1,data2Store1,16);
                end
                
                t1.close;
                
                %We also write a few info about the calibrated data in a
                %text file
                if MP
                    fprintf(fid,...
                        'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
                        i,cal.inFocus1(cal.neworder1(i)).cam,...
                        cal.inFocus1(cal.neworder1(i)).ch,...
                        cal.ROI1(cal.neworder1(i),1),...
                        cal.ROI1(cal.neworder1(i),1)+...
                        cal.ROI1(cal.neworder1(i),3),...
                        cal.inFocus1(cal.neworder1(i)).relZPos);
                    calib1.camConfig = obj.cal2D.camConfig{1,1};
                else
                    fprintf(fid,...
                        'No Calibration was performed on this data as only a single plane was provided. It is likely to be coming from the widefield setup');
                    calib1.camConfig{1,1} = 'None';
                end
                if isTrans1
                else
                    if MP
                        calib1.oRelZPos(idx2Plane) =  cal.inFocus1(cal.neworder1(i)).relZPos;
                    else
                        calib1.oRelZPos(idx2Plane) = 0;
                    end
                    idx2Plane = idx2Plane+1;
                end
                
                
            end
            calib1.Width  = size(data{1,1},2);
            calib1.Height = size(data{1,1},1);
            calib1.Backgrounds = Bg.Ch1;
            calib = calib1;
            fclose(fid);
            fName = [calDir1 filesep 'calibrated1.mat'];
            save(fName,'calib');
            
            if MP == 1
                if cal.multiModal == true
                   calDir2 = [obj.raw.movInfo.Path filesep 'calibrated2'];
                   calTransDir2 = [calDir2 filesep 'Transmission'];
                   mkdir(calDir2);
                   mkdir(calTransDir2);
                   fid2 = fopen([calDir2 filesep 'CalibratedInfo.txt'],'w');
                   fprintf(fid2,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
                   calib2.mainPath = calDir2;
                   calib2.nPlanes = sum(~isTransmission{2,1});
                   idx2Plane = 1;
    
                   for i = 1:size(data{2,1},3)
                        data2Store = squeeze(data{2,1}(:,:,i,:));
                        for j = 1:size(data2Store, 3)
                            Min = min(data2Store(:,:,j), [], 'all');
                            Min2 = mean(data2Store(:,:,j), 'all');
                            %tform = simtform2d(cal.Transformation{i,1}.Scale, cal.Transformation{i,1}.RotationAngle, cal.Transformation{i,1}.Translation);
                            %data2 = imwarp(double(data2Store(:,:,j)),tform,"OutputView",imref2d(size(double(data2Store1(:,:,j)))));
                            data2 = data2Store(:,:,j);
                            data2(data2 < Min) = Min2;
                            data2Store2(:,:,j) = uint16(data2);
    
                        end
                        isTrans2 = isTransmission{2,1}(i);
                        if strcmpi(obj.info.type,'transmission')
                            data2Store2 = imcomplement(data2Store2);
                        end
                        fieldN = sprintf('plane%d',i);
                        fName = sprintf('calibratedPlane%d.tif',i);
                        if isTrans2
                            fPathTiff = [calTransDir2 filesep fName];
                            calib2.transPath.(fieldN) = fPathTiff;
                        else
                            fPathTiff = [calDir2 filesep fName];
                            calib2.filePath.(fieldN) = fPathTiff;
                        end
                        
                        calib2.nFrames = maxFrame;
                        t2 = Tiff(fPathTiff, 'a');
                        
                        if isa(data2Store2,'uint32')
                            t2 = dataStorage.writeTiff(t2,data2Store2,32);
                        else
                            t2 = dataStorage.writeTiff(t2,data2Store2,16);
                        end
                        
                        t2.close;
                        
                        %We also write a few info about the calibrated data in a
                        %text file
                        if MP
                            fprintf(fid,...
                                'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
                                i+8,cal.inFocus2(cal.neworder2(i)).cam,...
                                cal.inFocus2(cal.neworder2(i)).ch,...
                                cal.ROI2(cal.neworder2(i),1),...
                                cal.ROI2(cal.neworder2(i),1)+...
                                cal.ROI2(cal.neworder2(i),3),...
                                cal.inFocus2(cal.neworder2(i)).relZPos);
                            calib2.camConfig = obj.cal2D.camConfig{1,2};
                        else
                            fprintf(fid,...
                                'No Calibration was performed on this data as only a single plane was provided. It is likely to be coming from the widefield setup');
                            calib2.camConfig{2,1} = 'None';
                        end
                        if isTrans2
                        else
                            if MP
                                calib2.oRelZPos(idx2Plane) =  cal.inFocus2(cal.neworder2(i)).relZPos;
                            else
                                calib2.oRelZPos(idx2Plane) = 0;
                            end
                            idx2Plane = idx2Plane+1;
                        end
                        
                        
                    end
                    calib2.Width  = size(data2,2);
                    calib2.Height = size(data2,1);
                    calib2.Backgrounds = Bg.Ch2;
                    calib2.RatioCh1Ch2 = Bg.Ratio;
                    fclose(fid2);
                    calib = struct([]);
                    calib = calib2;
                    fName = [calDir2 filesep 'calibrated2.mat'];
                    save(fName,'calib');
    
                else
                end
            
            calib = struct([]);

            calib = {calib1, calib2};
            end
        end    
    end
end

