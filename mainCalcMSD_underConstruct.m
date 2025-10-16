clc ;
clear ;
close all;

%% USER INPUT
[FilePath, Experiment, FilenameRaw, Dimension, expTime, Temp, Radius1, Radius2, DiffFit, MinSize, Ext, ParticleType, path2RotCal, CutTraces, ExpModel] = UserInput.CalcMSDinfoGUI;

%% Loading
f = waitbar(0, 'Initializing...');
MainFolder = dir(FilePath);
name = append(Experiment, Ext);
MainFolder([MainFolder.isdir] ~= 1) = [];

AllMovieResults = [];
CorrectionsTheta = [];
CorrectionsPhi = [];
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
            if strcmp(Experiment, 'Dual color tracking')
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
            % 
            % if strcmp(Path(end-4:end), 'n0_1')
            %     Temp = 303.15;
            % elseif strcmp(Path(end-4:end), 'n2_1')
            %     Temp = 304.15;
            % elseif strcmp(Path(end-4:end), 'n4_1')
            %     Temp = 305.15;
            % elseif strcmp(Path(end-4:end), 'n6_1')
            %     Temp = 306.15;
            % elseif strcmp(Path(end-4:end), 'n8_1')
            %     Temp = 307.15;
            % elseif strcmp(Path(end-4:end), '10_1')
            %     Temp = 308.15;
            % elseif strcmp(Path(end-4:end), '13_1')
            %     Temp = 296.15;
            % end

        %% Processing
            if ~strcmp(Experiment, 'Rotational Tracking')
                %%% cut up traces
                if ~isnan(CutTraces)
                    Length = 100;
                    for step = 1:8
                        CurrTraceCell = {};
                        for i = 1:size(data.traces)
                            CurrTrace = data.traces{i, 1};
                            CurrTraceCutted = CurrTrace(ismember(CurrTrace.t, (step*Length-(Length-1):step*Length)), :);    
                            if ~isempty(CurrTraceCutted)
                                CurrTraceCell{end+1,1} = CurrTraceCutted;
                            end
                        end
                        dataMatrix{step, 1} = CurrTraceCell;
                    end
                else
                    try
                        dataMatrix{1,1} = data.traces;
                    catch
                        dataMatrix{1,1} = data;
                    end
                end

                for k = 1:size(dataMatrix, 1)
                    DataCurrent = dataMatrix{k, 1};
                    try
                        allHeight = cellfun(@height,dataMatrix{k,1}(:,1));
                        idx = allHeight>MinSize;
                        currMov = dataMatrix{k,1}(idx, 1);
                    catch
                        allHeight = cellfun(@height,dataMatrix(:,1));
                        idx = allHeight>MinSize;
                        currMov = data.traces(idx, 1);
                    end
                    if ~isempty(currMov)
                    
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
                            waitbar(i./length(currMov), f, append('Diffusion & microrheology - Movie ', num2str(j-2), ' out of ', num2str(size(MainFolder,1) - 2), ' step ', num2str(l)));
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
                                try
                                    allRes(i).Phase = nanmean(currPart.Phase);
                                    allRes(i).IntPhaseCh = nanmean(currPart.IntPhaseCh);
                                    allRes(i).GradientMagnitude = nanmean(currPart.GradientMagnitude);
                                    allRes(i).LocalVariance = nanmean(currPart.LocalVariance);
                                    allRes(i).SharpnessLaplacian = nanmean(currPart.SharpnessLaplacian);
                                catch
                                    allRes(i).Phase = NaN;
                                    allRes(i).IntPhaseCh = NaN;
                                    allRes(i).GradientMagnitude = NaN;
                                    allRes(i).LocalVariance = NaN;
                                    allRes(i).SharpnessLaplacian = NaN;
                                end
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
                            if ~isnan(CutTraces)
                                name = append('msdResSegmentation_Step', num2str(k), Ext);
                            else
                                name = append('msdResSegmentation', Ext);
                            end
                        elseif strcmp(Experiment, 'Tracking-Phase')
                            if ~isnan(CutTraces)
                                name = append('msdResPhase_Step', num2str(k), Ext);
                            else
                                name = append('msdResPhase', Ext);
                            end
                        elseif strcmp(Experiment, 'Tracking')
                            if ~isnan(CutTraces)
                                name = append('msdRes_Step', num2str(k), '_', Filename(end), Ext);
                            else
                                name = append('msdRes', Filename(end), Ext);
                            end
                        elseif strcmp(Experiment, 'Dual color tracking')
                            if ~isnan(CutTraces)
                                name = append('msdRes_Step', num2str(k), '_', Filename(end), Ext);
                            else
                                name = append('msdRes', Filename(end), Ext);
                            end
                        else
                            error('Please specify type of experiment: Tracking, Tracking-Segmentation, Tracking-Phase');
                        end
                    
                        FileNameToSave = append(Path, filesep, name);
                        save(FileNameToSave,'allRes');
                        disp(append('== Data succesfully saved - movie ', num2str(j-2), ' out of ', num2str(size(MainFolder, 1)-2), ' =='));
                    
                        AllMovieResults = [AllMovieResults, allRes];
                    else
                        disp(append('No traces found that are longer than MinSize (', num2str(MinSize), ' datapoints)'))
                    end
                end
            else
                %% This is for rotational tracking
                allHeight = cellfun(@height,data.Coord1);
                idx = allHeight>MinSize;
                currMov = data(idx,:);
                allRes = struct('GTheta',0,'GPhi',0,'tau',0,'DTheta',0,'DPhi',0,'Dr',0,...
                    'nTheta',0,'nPhi',0,'nr',0,'vTheta',0,'vPhi',0,'vr',0, 'TauCTheta', 0, 'TauCPhi', 0, 'num', 0);
                allRes(size(currMov, 1)).masdTheta = [];
                maxLength = max(allHeight);
                allGTheta = zeros(size(currMov, 1),maxLength-1);
                allGPhi = allGTheta;
                
                if strcmp(MainFolder(j).name(1:4), '3min')
                    CorrFactor = 31;
                elseif strcmp(MainFolder(j).name(1:4), '5min')
                    CorrFactor = 37;
                elseif strcmp(MainFolder(j).name(1:4), '7min')
                    CorrFactor = 29;
                elseif strcmp(MainFolder(j).name(1:4), '9min')
                    CorrFactor = 1;
                elseif strcmp(MainFolder(j).name(1:4), '11mi')
                    CorrFactor = 1.7;
                elseif strcmp(MainFolder(j).name(1:4), '13mi')
                    CorrFactor = 2.5;
                elseif strcmp(MainFolder(j).name(1:4), '15mi')
                    CorrFactor = 12;
                elseif strcmp(MainFolder(j).name(1:4), '17mi')
                    CorrFactor = 19;
                elseif strcmp(MainFolder(j).name(1:4), '19mi')
                    CorrFactor = 22.5;
                end
                
                for i = 1:size(currMov,1)
                    waitbar(i./size(currMov,1), f, append('Diffusion & microrheology - Movie ', num2str(j-2), ' out of ', num2str(size(MainFolder,1) - 2)));
                    Diff = [];
                    TotInt = [];
                    GTheta = [];
                    GPhi = [];
                    tau = [];
                    currPart = currMov(i,:);
                
                    %%% calculate angels
                    TotInt = currPart.Int1{1,1} + currPart.Int2{1,1};
                    I1 = currPart.Int1{1,1} ./ TotInt;
                    I2 = currPart.Int2{1,1} ./ TotInt;
                    Diff = I1 - I2; 

                    eta_truth = 484;
                    %%% for Theta
                    [GTheta,tau] = MSD.Rotational.GetAutoCorrelation(Diff,100,expTime);                   
                    allGTheta(i,1:size(GTheta, 2)) = GTheta';
                    tau(isnan(GTheta)) = [];
                    GTheta(isnan(GTheta)) = [];
                    if size(GTheta, 2) > MinSize
                        try
                            if strcmp(ExpModel, 'Test')
                                [Model] = MSD.Rotational.TestModels(GTheta, tau);
                                %[DTheta] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, Model,0);
                                [DTheta, corrParamsTheta,TauCTheta] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, Model, 0, eta_truth, 'Bipyramid', 'Theta', MinSize);
                            else
                                %[DTheta] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, ExpModel,0);
                                [DTheta, corrParamsTheta, TauCTheta] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, ExpModel, 0, eta_truth, 'Bipyramid', 'Theta', MinSize);
                            end
                            nTheta  = MSD.Rotational.getViscosity(DTheta.*25195./CorrFactor,Radius1,ParticleType, Temp, 'Theta');
                            Dcalctheta = (3*1.380649*10^(-23)*Temp*log(Radius1(1)/Radius1(2)))./(pi*nTheta*Radius1(1)^3)*1000;
                            TauCTheta = (((((1./(2*pi*Dcalctheta))*2*pi))./8));
                            vTheta = MSD.Rotational.GetRotationalSpeed(GTheta, tau);
                        catch
                            DTheta = NaN;
                            nTheta = NaN;
                            vTheta = NaN;
                            TauCTheta = NaN;
                            corrParamsTheta = [];
                        end
                    else
                        DTheta = NaN;
                        nTheta = NaN;
                        vTheta = NaN;
                        TauCTheta = NaN;
                        corrParamsTheta = [];
                    end

                    %%% for Phi
                    [GPhi,tau] = MSD.Rotational.GetAutoCorrelation(TotInt,100,expTime);
                    allGPhi(i,1:size(GPhi, 2)) = GPhi';
                    tau(isnan(GPhi)) = [];
                    GPhi(isnan(GPhi)) = [];
                    if size(GPhi,2) > MinSize
                        try
                            if strcmp(ExpModel, 'Test')
                                [Model] = MSD.Rotational.TestModels(GPhi, tau);
                                %[DPhi] = MSD.Rotational.GetDiffusion(GPhi, tau, Radius1, Temp, Dimension, Model,0);
                                [DPhi, corrParamsPhi, TauCPhi] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, Model, 0, eta_truth, 'Bipyramid', 'Theta', MinSize);
                            else
                                %[DPhi] = MSD.Rotational.GetDiffusion(GPhi, tau, Radius1, Temp, Dimension, ExpModel,0);
                                [DPhi, corrParamsPhi, TauCPhi] = MSD.Rotational.GetDiffusion(GTheta, tau, Radius1, Temp, Dimension, ExpModel, 0, eta_truth, 'Bipyramid', 'Theta', MinSize);
                            end
                            nPhi  = MSD.Rotational.getViscosity(DPhi.*25195./CorrFactor,Radius1,ParticleType, Temp, 'Phi');
                            DcalcPhi = (3*1.380649*10^(-23)*Temp*log(Radius1(1)/Radius1(2)))./(pi*nPhi*Radius1(1)^3)*1000;
                            TauCPhi = (((((1./(2*pi*DcalcPhi))*2*pi))./8));
                            vPhi = MSD.Rotational.GetRotationalSpeed(GPhi, tau);
                        catch
                            DPhi = NaN;
                            nPhi = NaN;
                            vPhi = NaN;
                            TauCPhi = NaN;
                            corrParamsPhi = [];
                        end
                    else
                        DPhi = NaN;
                        nPhi = NaN;
                        vPhi = NaN;
                        TauCPhi = NaN;
                        corrParamsPhi = [];
                    end
                
                
                    %For both
                    Dr = mean([DTheta, DPhi]);
                    nr = mean([nTheta, nPhi]);
                    vr = sqrt(vTheta^2 + vPhi^2);
                
                    allRes(i).GTheta = GTheta;% in rad^2 
                    allRes(i).GPhi = GPhi;
                    allRes(i).tau = tau; % in sec
                
                    allRes(i).DTheta   = DTheta;% in rad^2 /sec
                    allRes(i).DPhi   = DPhi;% in rad^2 /sec
                    allRes(i).Dr  = Dr;% in rad^2 /sec
                
                    allRes(i).nTheta   = nTheta;
                    allRes(i).nPhi   = nPhi;
                    allRes(i).nr   = nr;
                    
                    allRes(i).vTheta   = vTheta;
                    allRes(i).vPhi   = vPhi;
                    allRes(i).vr   = vr;

                    allRes(i).TauCTheta = TauCTheta;
                    allRes(i).TauCPhi = TauCPhi;
                    
                    allRes(i).num  = length(GTheta);

                    CorrectionsTheta = [CorrectionsTheta; corrParamsTheta];
                    CorrectionsPhi = [CorrectionsPhi; corrParamsPhi];
                end
                
                %%
                
                DR   = nanmean([allRes.Dr]);
                nR   = nanmean([allRes.nr]);
                disp(['The diffusion coefficient is ' num2str(DR) 'µm^2s^-^1' '\n'...
                    'the viscosity is ' num2str(nR) ' cp' '\n' ...
                    'for movie '  MainFolder(j).name]);
                %%
                filenameAllRes = [folder(1).folder filesep 'msadRes.mat'];
                save(filenameAllRes,'allRes');
    
                AllMovieResults = [AllMovieResults, allRes];
            end
        end
    catch
        disp(append('Failed to calculate movie ', MainFolder(j).name));
    end
end
% mean(CorrectionsTheta,1)
% mean(CorrectionsPhi,1)
close(f)
disp(append('TauTheta = ', num2str(mean([AllMovieResults.TauCTheta], 'omitnan'))));
disp(append('TauPhi = ', num2str(mean([AllMovieResults.TauCPhi], 'omitnan'))));
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



