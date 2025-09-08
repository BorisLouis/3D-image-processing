clc ;
clear ;
close all;

%% USER INPUT
[FilePath, Experiment, FilenameRaw, Dimension, expTime, Temp, Radius1, Radius2, DiffFit, MinSize, Ext, ParticleType, path2RotCal] = UserInput.CalcMSDinfoGUI;

%% Loading
f = waitbar(0, 'Initializing...');
MainFolder = dir(FilePath);
name = append(Experiment, Ext);
MainFolder([MainFolder.isdir] ~= 1) = [];

AllMovieResults = [];
for j = 3 : size(MainFolder,1)
    try
        folder = dir(append(MainFolder(j).folder, filesep, append(MainFolder(j).name)));
        idx = contains({folder.name},FilenameRaw);
        folder(~idx) = [];

        if strcmp(Experiment, 'Dual color tracking')
            Loop = 2;
        else
            Loop = 1;
        end
        for l = 1:Loop
            if Loop == 2
                if l == 1
                    Filename = append(FilenameRaw, '1');
                    Radius = Radius1;
                elseif l == 2
                    Filename = append(FilenameRaw, '2');
                    Radius = Radius2;
                end
            else
                Filename = FilenameRaw;
                Radius = Radius1;
            end
    
            f2Load = [folder(1).folder filesep Filename '.mat'];
            Path = folder(1).folder;
            
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});

            if strcmp(Path(end-4:end), 'n0_1')
                Temp = 303.15;
            elseif strcmp(Path(end-4:end), 'n2_1')
                Temp = 304.15;
            elseif strcmp(Path(end-4:end), 'n4_1')
                Temp = 305.15;
            elseif strcmp(Path(end-4:end), 'n6_1')
                Temp = 306.15;
            elseif strcmp(Path(end-4:end), 'n8_1')
                Temp = 307.15;
            elseif strcmp(Path(end-4:end), '10_1')
                Temp = 308.15;
            elseif strcmp(Path(end-4:end), '13_1')
                Temp = 296.15;
            end
        
        %% Processing
            if ~strcmp(Experiment, 'Rotational Tracking')
                try
                    allHeight = cellfun(@height,data.traces(:,1));
                    idx = allHeight>MinSize;
                    currMov = data(idx, 1);
                catch
                    allHeight = cellfun(@height,data(:,1));
                    idx = allHeight>MinSize;
                    currMov = data(idx, 1);
                end
                if isempty(currMov)
                    error(append('No traces found that are longer than MinSize (', num2str(MinSize), ' datapoints)'))
                end
                allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
                    'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
                    'vZ',0,'vR',0);
                if strcmp(Experiment, 'Tracking-Segmentation')
                    allRes.Mask = 0;
                elseif strcmp(Experiment, 'Tracking-Phase')
                    allRes.Phase = 0;
                end
                allRes(length(currMov)).msdX = [];
                maxLength = max(allHeight);
                allMSDX = zeros(length(currMov),maxLength-1);
                allMSDY = allMSDX;
                allMSDZ = allMSDY;
                allMSDR = allMSDY;
                for i = 1:length(currMov)
                    waitbar(i./length(currMov), f, append('Calculating diffusion & doing microrheology - Movie ', num2str(j-2), ' out of ', num2str(size(MainFolder,1) - 2)));
                    currPart = currMov{i};
                
                    coord = [currPart.col, currPart.row, currPart.z];
                    CM = mean(coord,1);
                    coord = coord-CM;
                
                    %in X
                    AvStep = MSD.getAvStepSize(coord(:,1)/10^3); 
                    msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
                    tau = (1:length(msdx))'*expTime;
                    allMSDX(i,1:length(msdx)) = msdx;
                    DX   = MSD.getDiffCoeff(msdx,tau,DiffFit,'1D');
                    nX   = MSD.getViscosity(DX,Radius,Temp);
                    %aX   = MSD.getDiffTypeAlpha(msdx,expTime);
                    aX   = MSD.getDiffTypeAlpha2(msdx,expTime, AvStep);
                    vX   = abs(coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime)); %um/s
                
                    %inY
                    AvStep = MSD.getAvStepSize(coord(:,2)/10^3); 
                    msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
                    allMSDY(i,1:length(msdy)) = msdy;
                    DY   = MSD.getDiffCoeff(msdy,tau,DiffFit,'1D');
                    nY   = MSD.getViscosity(DY,Radius,Temp);
                    %aY   = MSD.getDiffTypeAlpha(msdy,expTime);
                    aY   = MSD.getDiffTypeAlpha2(msdy,expTime, AvStep);
                    vY   = abs(coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime)); %um/s
                
                    %inZ
                    if strcmp(Dimension, '3D')
                        AvStep = MSD.getAvStepSize(coord(:,3)/10^3); 
                        msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
                        allMSDZ(i,1:length(msdz)) = msdz;
                        DZ   = MSD.getDiffCoeff(msdz,tau,DiffFit,'1D');
                        nZ   = MSD.getViscosity(DZ,Radius,Temp);
                        %aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
                        aZ   = MSD.getDiffTypeAlpha2(msdz,expTime, AvStep);
                        vZ   = abs(coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime)); %um/s
                    else
                        msdz = [];
                        allMSDZ(i, 1:1) = NaN;
                        DZ = NaN;
                        nZ = NaN;
                        aZ = NaN;
                        vZ = NaN;
                    end
                        
                
                
                    %inR
                    if strcmp(Dimension, '3D')
                        AvStep = MSD.getAvStepSize(coord(:, 1:3)/10^3); 
                        msdr = MSD.calc(coord(:, 1:3)/10^3);%convert to um;
                        allMSDR(i,1:length(msdr)) = msdr;
                        DR   = MSD.getDiffCoeff(msdr,tau,DiffFit,'3D');
                    elseif strcmp(Dimension, '2D')
                        AvStep = MSD.getAvStepSize(coord(:, 1:2)/10^3); 
                        msdr = MSD.calc(coord(:, 1:2)/10^3);%convert to um;
                        allMSDR(i,1:length(msdr)) = msdr;
                        DR   = MSD.getDiffCoeff(msdr,tau,DiffFit,'2D');
                    end          
                    nR   = MSD.getViscosity(DR,Radius,Temp);
                    %aR   = MSD.getDiffTypeAlpha(msdr,expTime);
                    aR   = MSD.getDiffTypeAlpha2(msdr,expTime, AvStep);
                    dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2 +...
                        (coord(1,3)-coord(end,3))^2);
                    vR = dR/10^3/(length(coord)*expTime); %um/s
                
                    allRes(i).msdX = msdx;% in um^2
                    allRes(i).msdY = msdy;
                    allRes(i).msdZ = msdz;
                    allRes(i).msdR = msdr;
                    allRes(i).tau = tau; % in sec
            
                    allRes(i).DX   = DX;% in um^2 /sec
                    allRes(i).DY   = DY;% in um^2 /sec
                    allRes(i).DZ   = DZ;% in um^2 /sec
                    allRes(i).DR   = DR;% in um^2 /sec
                
                    allRes(i).nX   = nX;
                    allRes(i).nY   = nY;
                    allRes(i).nZ   = nZ;
                    allRes(i).nR   = nR;
                    
                    allRes(i).aX   = aX;
                    allRes(i).aY   = aY;
                    allRes(i).aZ   = aZ;
                    allRes(i).aR   = aR;
                    
                    allRes(i).vX   = vX;
                    allRes(i).vY   = vY;
                    allRes(i).vZ   = vZ;
                    allRes(i).vR   = vR;
                    
                    
                    allRes(i).num  = length(msdx);
            
                    if strcmp(Experiment, 'Tracking-Segmentation')
                        allRes(i).Mask = round(mean(currPart.InSegment));
                    elseif strcmp(Experiment, 'Tracking-Phase')
                        allRes(i).Phase = mean(currPart.Phase);
                        allRes(i).IntPhaseCh = mean(currPart.IntPhaseCh);
                        allRes(i).GradientMagnitude = mean(currPart.GradientMagnitude);
                        allRes(i).LocalVariance = mean(currPart.LocalVariance);
                        allRes(i).SharpnessLaplacian = mean(currPart.SharpnessLaplacian);
                    end
                end
            
                %%
                allMSDR(allMSDR == 0) = NaN;
                meanMSDR = nanmean(allMSDR,1);
                tau = (1:length(meanMSDR))'*expTime;
                DR   = MSD.getDiffCoeff(meanMSDR,tau,DiffFit,Dimension);
                nR   = MSD.getViscosity(DR,Radius,Temp);
                
                disp(['The diffusion coefficient is ', num2str(DR), ' \mum^2/s and the viscosity is ' num2str(nR) ' cp']);
                %%
            
                if strcmp(Experiment, 'Tracking-Segmentation')
                    name = append('msdResSegmentation', Ext);
                elseif strcmp(Experiment, 'Tracking-Phase')
                    name = append('msdResPhase', Ext);
                elseif strcmp(Experiment, 'Tracking')
                    name = append('msdRes', Filename(end), Ext);
                elseif strcmp(Experiment, 'Dual color tracking')
                    name = append('msdRes', Filename(end), Ext);
                else
                    error('Please specify type of experiment: Tracking, Tracking-Segmentation, Tracking-Phase');
                end
            
                FileNameToSave = append(Path, filesep, name);
                save(FileNameToSave,'allRes');
                disp(append('== Data succesfully saved - movie ', num2str(j-2), ' out of ', num2str(size(MainFolder, 1)-2), ' =='));
            
                AllMovieResults = [AllMovieResults, allRes];
            else
                %% This is for rotational tracking
                allHeight = cellfun(@height,data.Coord1);
                idx = allHeight>MinSize;
                currMov = data(idx,:);
                allRes = struct('msadTheta',0,'msadPhi',0,'msadr',0,'tau',0,'DTheta',0,'DPhi',0,'Dr',0,...
                    'nTheta',0,'nPhi',0,'nr',0,'vTheta',0,'vPhi',0,'vr',0);
                allRes(size(currMov, 1)).masdTheta = [];
                maxLength = max(allHeight);
                allmsadTheta = zeros(size(currMov, 1),maxLength-1);
                allmsadPhi = allmsadTheta;
                allmsad = allmsadPhi;
    
                % calibration = load(append(Path2Rotcal, filesep, 'RotCalib.mat'));
                % name = fieldnames(calibration);
                % calibration = calibration.(name{1,1});
                
                for i = 1:size(currMov,1)
                    currPart = currMov(i,:);
                
                    %%% calculate angels
                    TotInt = currPart.Int1{1,1} + currPart.Int2{1,1};
                    I1 = currPart.Int1{1,1}./TotInt;
                    I2 = currPart.Int2{1,1}./TotInt;
                    Diff = I1 - I2;
                
                    Theta = 0.5*real(acos(Diff./calibration.I0_mean));
                    Phi = real(acos(sqrt(TotInt./calibration.TotI0_mean)));
                
                    coord = [Theta, Phi];
                
                    tau = currPart.Time{1,1};
    
                    %For Theta
                    msadTheta = MSD.Rotational.calc(coord(:,1), tau, expTime);
                    allmsadTheta(i,1:size(msadTheta, 2)) = msadTheta(1,:);
                    DTheta   = MSD.Rotational.getDiffCoeff(msadTheta,tau,DiffFit,'2D');
                    nTheta   = MSD.Rotational.getViscosity(DTheta,Radius,ParticleType, Temp);
                    vTheta   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s
                
                    %For Phi
                    msadPhi = MSD.Rotational.calc(coord(:,2), tau, expTime);
                    allmsadPhi(i,1:size(msadPhi, 2)) = msadPhi(1,:);
                    DPhi   = MSD.Rotational.getDiffCoeff(msadPhi,tau,DiffFit,'2D');
                    nPhi   = MSD.Rotational.getViscosity(DPhi,Radius,ParticleType, Temp);
                    vPhi   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s
                
                    %For both
                    msadr = MSD.Rotational.calc(coord, tau, expTime);%convert to um;
                    allMSDR(i,1:size(msadr,2)) = msadr(1,:);
                    DR   = MSD.Rotational.getDiffCoeff(msadr,tau,DiffFit,'3D');
                    nR   = MSD.Rotational.getViscosity(DR,Radius,ParticleType,Temp);
                    dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2);
                    vR = dR/10^3/(length(coord)*expTime); %rad/s
                
                    allRes(i).msadTheta = msadTheta;% in rad^2 
                    allRes(i).msadPhi = msadPhi;
                    allRes(i).msadr = msadr;
                    allRes(i).tau = tau; % in sec
                
                    allRes(i).DTheta   = DTheta;% in rad^2 /sec
                    allRes(i).DPhi   = DPhi;% in rad^2 /sec
                    allRes(i).Dr  = DR;% in rad^2 /sec
                
                    allRes(i).nTheta   = nTheta;
                    allRes(i).nPhi   = nPhi;
                    allRes(i).nr   = nR;
                    
                    allRes(i).vTheta   = vTheta;
                    allRes(i).vPhi   = vPhi;
                    allRes(i).vr   = vR;
                    
                    allRes(i).num  = length(msadTheta);
                end
                
                %%
                
                meanMSDR = nanmean(allMSDR,1);
                tau = (1:length(meanMSDR))'*expTime;
                DR   = nanmean([allRes.DR]);
                nR   = nanmean([allRes.nR]);
                disp(['The diffusion coefficient is ' num2str(DR) 'µm^2s^-^1' '\n'...
                    'the viscosity is ' num2str(nR) ' cp' '\n' ...
                    'for movie '  MainFolder(j).name]);
                %%
                if strcmp(Experiment, 'Dual color tracking')
                    if l == 1
                        filenameAllRes = [folder(1).folder filesep 'msadRes1.mat'];
                    elseif l == 2
                        filenameAllRes = [folder(1).folder filesep 'msadRes2.mat'];
                    end
                else
                    filenameAllRes = [folder(1).folder filesep 'msadRes.mat'];
                end
                save(filenameAllRes,'allRes');
    
                AllMovieResults = [AllMovieResults, allRes];
            end
        end
    catch
        disp(append('Failed to calculate movie ', MainFolder(j).name));
    end
end
close(f)

save(FilePath, "AllMovieResults");

if strcmp(Experiment, 'Tracking-Segmentation')
    name = 'msdResSegmentation';
    AllMovieResultsMask = AllMovieResults([AllMovieResults.Mask] == 1);
    AllMask = rmfield(AllMovieResultsMask, {'msdX', 'msdY', 'msdZ', 'msdR', 'tau'});
    AllMask = struct2table(AllMask);
    writetable(AllMask,append(FilePath, filesep, name, '.xlsx'),'Sheet','data - mask');
    AllMovieResultsNoMask = AllMovieResults([AllMovieResults.Mask] == 0);
    AllNoMask = rmfield(AllMovieResultsNoMask, {'msdX', 'msdY', 'msdZ', 'msdR', 'tau'});
    AllNoMask = struct2table(AllNoMask);
    writetable(AllNoMask,append(FilePath, filesep, name, '.xlsx'),'Sheet','data - no mask');
    writecell({AllMovieResultsMask.msdR}',append(FilePath, filesep, name, '.xlsx'),'Sheet','msdR - mask');
    writecell({AllMovieResultsNoMask.msdR}',append(FilePath, filesep, name, '.xlsx'),'Sheet','msdR - no mask');
elseif strcmp(Experiment, 'Tracking-Phase')
    name = 'msdResPhase';
elseif strcmp(Experiment, 'Tracking')
    name = 'msdRes';
elseif strcmp(Experiment, 'Rotational Tracking')
    name = 'msadResRot';
end
AllMovieResultsTable = struct2table(AllMovieResults);
writetable(AllMovieResultsTable,append(FilePath, filesep, name, '.xlsx'),'Sheet',1,'Range','D1');



