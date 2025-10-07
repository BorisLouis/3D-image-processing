classdef MPDDMMovie < Core.MPMovie
    %MPDDMMOVIE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        AllFrames
        DDMOutput
        FitOutput
        MSD
    end
    
    methods
        function obj = MPDDMMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
        end
        
        function getFullFrames(obj)
            for c = 1:obj.calibrated{1, 1}.nPlanes
                h = waitbar(0, 'initializing');
                for frame = 1:obj.raw.movInfo.maxFrame
                    Frame = double(Load.Movie.tif.getframes(obj.calibrated{1, 1}.filePath.(append('plane', num2str(c))), frame));
                    if strcmp(obj.info.ddmParam.CorrectBleaching, 'on')
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame + bleaching correction -- frame ',...
                            num2str(frame), '/', num2str(obj.raw.movInfo.maxFrame), ' plane ', num2str(c)));
                        if frame == 1
                            RefInt = mean(Frame(Frame ~= 0));
                        else
                            Frame = Frame * (RefInt ./ mean(Frame(Frame ~= 0)));
                        end
                    else
                        waitbar(frame./obj.raw.movInfo.maxFrame, h, append('Load frame - no bleaching correction -- frame ',...
                            num2str(frame), '/', num2str(obj.raw.movInfo.maxFrame)));
                    end
                    obj.AllFrames(:,:,frame,c) = Frame;
                end
                close(h)
            end
        end

        function mainDDM(obj,varargin)
            for c = obj.calibrated{1,1}.nPlanes
                DDMOutputfile = append(obj.calibrated{1, 1}.mainPath, filesep, 'DDMOutput' , num2str(c), '.mat');
                if ~exist(DDMOutputfile)
                    run = 1;
                else
                    if strcmp(obj.info.runMethod, 'run')
                        run = 1;
                        delete(DDMOutputfile);
                    else
                        run = 0;
                    end
                end

                if run == 1
                    p = inputParser;  
                    addOptional(p, 'ROI',  [1, size(obj.AllFrames,1), size(obj.AllFrames,1);  %Default ROI as image size
                                            1, size(obj.AllFrames,2), size(obj.AllFrames,2);
                                            1 ,size(obj.AllFrames,3) ,size(obj.AllFrames,3)]);
                    addOptional(p, 'Padsize',zeros(1, 3));
                    addOptional(p, 'NumBins', 200);
                    addOptional(p, 'CriticalAngle',0);
                    FrameSize = size(obj.AllFrames);
                    FrameSize(3) = 1;
                    ExpTime = obj.info.ddmParam.ExpTime; 
                    DDMOutput = [];
                    AnisotropyOutput = [];

                    for dt=1:obj.raw.movInfo.maxFrame-1
                        AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3));
                        [AvgFFT] = obj.CalculateDelta(AvgFFT,FrameSize,dt, p.Results.ROI, c); 
                        
                        [RadiallyAveragedDDMSignal, AnisotropyValues]=  obj.AverageRadialy3D(AvgFFT,FrameSize, p.Results.CriticalAngle);
                        % QImages{dt, :} = reconstructions;
                        DDMOutput(:,1) =[NaN ; RadiallyAveragedDDMSignal(:,1)];
                        DDMOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)];
    
                        AnisotropyOutput(:,1) = [NaN ; RadiallyAveragedDDMSignal(:,1)];
                        AnisotropyOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)]; 
                        waitbar(dt./obj.DDMInfo.FramesToAnalyze,f,append('Calculating AvgFFT frame by frame - Timelag ', num2str(dt),...
                            ' out of ', num2str(obj.DDMInfo.nFrames-1)));
                    end
                else
                    disp('Found DDMOuput file - loading it');
                    DDMOutput = load(DDMOutputfile);
                    DDMOutput = DDMOutput.DDMOutput;
                    obj.DDMOutput{c,1} = DDMOutput;
                    disp('Found DDMOuput file - Done');
                end
            end
        end

        function fitDDM(obj)
        end

        function MSDoptimisation(obj)
        end

         function [AverageDDMValueAtR, anisotropy_values] = AverageRadialy3D(obj, AvgFFT,FrameSize,critangle)
            % Averages the scattering function in 3d radially. Returns the
            % average values in function of time and q-vector.

            [RadialValueInQSpace, ValidRange]=obj.Get3DGrid(FrameSize,critangle);
            % BinSize = max(RadialValueInQSpace,[],'all')/NumBins; 
            NumBins = round((obj.DDMInfo.QmaxSelect - obj.DDMInfo.QminSelect)./(4*pi./(min(size(RadialValueInQSpace))*obj.DDMInfo.PixelSize))-1);
            MinQTracked = min(RadialValueInQSpace(:));
            Min = max(MinQTracked, obj.DDMInfo.QminSelect);
            BinSize = (obj.DDMInfo.QmaxSelect - Min)./NumBins;
            FoundRadii = 0;
            R = Min;
            next = 1;
            AverageDDMValueAtR = nan(NumBins,2);
            AvgFFT = gather(AvgFFT);

            [N, M] = size(AvgFFT);
            cx = floor(N/2) + 1;
            cy = floor(M/2) + 1;
            num_theta = 360;
            theta_vals = linspace(0, 2*pi, num_theta);

            %%% Calculate Average
            while ~isempty(FoundRadii)
                R = R+BinSize;
                if R <= obj.DDMInfo.QmaxSelect
                    if size(AvgFFT) == size(RadialValueInQSpace)
                        for t = 1:num_theta
                            x = round(cx + R * cos(theta_vals(t)));
                            y = round(cy + R * sin(theta_vals(t)));
                            if x >= 1 && x <= N && y >= 1 && y <= M
                                intensity_profile(t) = AvgFFT(x,y);
                            end
                        end

                        I_max = max(intensity_profile);
                        I_min = min(intensity_profile);
                        anisotropy_values(next,:) = [R, (I_max - I_min) / (I_max + I_min)];

                        FoundRadii = AvgFFT(RadialValueInQSpace>=R-BinSize & RadialValueInQSpace<R &  ValidRange );
                        AverageDDMValueAtR(next,:) = [R ,  nanmean(FoundRadii, 'all')];
                        next = next+1;

                    elseif size(AvgFFT)
                        error('Differential image and reference frame q vector do not have same dimensions')
                    end
                else
                    break
                end
            end
         end

         function [AvgFFT] =  CalculateDelta(obj,AvgFFT,FrameSize,dt,ROI,c)
            % DiffFrames = {};
            counts = 0;
            for t=1:obj.DDMInfo.nFrames-dt
                if obj.IsCudaDevice==1
                    FrameDelta =gpuArray(obj.AllFrames{c, t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{c,t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))); 
                else
                    FrameDelta =(obj.AllFrames{c, t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{c,t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2)));  
                end
                % DiffFrames{1, end+1} = fftshift(fft2(FrameDelta-mean(FrameDelta,'all')));
                % DiffFrames{2, end} = FrameDelta-mean(FrameDelta,'all');
                FrameDelta = abs( fftshift(fftn(FrameDelta-mean(FrameDelta,'all'),[FrameSize(1),FrameSize(2),FrameSize(3)]))).^2;                  
                AvgFFT = AvgFFT + FrameDelta;  
                counts = counts + 1;
            end
            AvgFFT = AvgFFT./counts;
        end 
    end
end

