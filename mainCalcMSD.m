clc ;
clear ;
close all;

%% USER INPUT
[FilePath, Experiment, Filename, Dimension, expTime, Temp, Radius, DiffFit, MinSize, Ext] = UserInput.CalcMSDinfoGUI;

%% Loading
f = waitbar(0, 'Initializing...');
MainFolder = dir(FilePath);
name = append(Experiment, Ext);
MainFolder([MainFolder.isdir] ~= 1) = [];

AllMovieResults = [];
for j = 3 : size(MainFolder,1)

    folder = dir(append(MainFolder(j).folder, filesep, append(MainFolder(j).name)));
    idx = contains({folder.name},Filename);
    folder(~idx) = [];
    
    f2Load = [folder(1).folder filesep folder(1).name];
    Path = folder(1).folder;
    
    tmpData = load(f2Load);
    name = fieldnames(tmpData);
    data = tmpData.(name{1});
    
    %% Processing
    allHeight = cellfun(@height,data(:,1));
    idx = allHeight>MinSize;
    currMov = data(idx, 1);
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
        msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
        tau = (1:length(msdx))'*expTime;
        allMSDX(i,1:length(msdx)) = msdx;
        DX   = MSD.getDiffCoeff(msdx,tau,DiffFit,'1D');
        nX   = MSD.getViscosity(DX,Radius,Temp);
        aX   = MSD.getDiffTypeAlpha(msdx,expTime);
        vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s
    
        %inY
        msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
        allMSDY(i,1:length(msdy)) = msdy;
        DY   = MSD.getDiffCoeff(msdy,tau,DiffFit,'1D');
        nY   = MSD.getViscosity(DY,Radius,Temp);
        aY   = MSD.getDiffTypeAlpha(msdy,expTime);
        vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s
    
        %inZ
        if strcmp(Dimension, '3D')
            msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
            allMSDZ(i,1:length(msdz)) = msdz;
            DZ   = MSD.getDiffCoeff(msdz,tau,fitRDiff,'1D');
            nZ   = MSD.getViscosity(DZ,R,T);
            aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
            vZ   = coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime); %um/s
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
            msdr = MSD.calc(coord(:, 1:3)/10^3);%convert to um;
            allMSDR(i,1:length(msdr)) = msdr;
            DR   = MSD.getDiffCoeff(msdr,tau,DiffFit,'3D');
        elseif strcmp(Dimension, '2D')
            msdr = MSD.calc(coord(:, 1:2)/10^3);%convert to um;
            allMSDR(i,1:length(msdr)) = msdr;
            DR   = MSD.getDiffCoeff(msdr,tau,DiffFit,'2D');
        end          
        nR   = MSD.getViscosity(DR,Radius,Temp);
        aR   = MSD.getDiffTypeAlpha(msdr,expTime);
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
            allRes(i).Phase = mean(currPart.InSegment);
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
        name = append('msdRes', Ext);
    else
        error('Please specify type of experiment: Tracking, Tracking-Segmentation, Tracking-Phase');
    end

    FileName = append(Path, filesep, name);
    save(FileName,'allRes');
    disp(append('== Data succesfully saved - movie ', num2str(j-2), ' out of ', num2str(size(MainFolder, 1)-2), ' =='));

    AllMovieResults = [AllMovieResults, allRes];
end
close(f)

FileName = append(FilePath, filesep, name);
save(FileName, "AllMovieResults");
if strcmp(Experiment, 'Tracking-Segmentation')
    name = 'msdResSegmentation';
elseif strcmp(Experiment, 'Tracking-Phase')
    name = 'msdResPhase';
elseif strcmp(Experiment, 'Tracking')
    name = 'msdRes';
end
AllMovieResultsTable = struct2table(AllMovieResults);
writetable(AllMovieResultsTable,append(FilePath, filesep, name, '.xlsx'),'Sheet',1,'Range','D1');



