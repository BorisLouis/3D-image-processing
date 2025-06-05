clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.010; %in sec
Temp = 296.15; %temperature in Kelvin
ParticleType = 'Bipyramid'; %Bipyramid, ellipsoid, rod, cilinder,...
R = [184, 92]; %Long axis, short axis in nm
fitRDiff = 3; %in number of data
minSize = 10; %frames
ext = '.mat';
path2RotCal = 'E:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';

%% Path info
MainFolder = 'E:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\Glycerol';
SubFolder = {'glycerol 80', 'glycerol 85', 'glycerol 90', 'glycerol 95', 'glycerol 100'}; % 'glycerol_80', 'glycerol_85', 'glycerol_90','glycerol_95', 
SubsubFolder = {'sample1', 'sample2', 'sample3','sample4', 'sample5'}; %


f = waitbar(0,'Initializing');
DiffIMatrix = [];
DiffTotMatrix = [];
for r = 1:numel(SubFolder)
    Visc = [];
    for o = 1:numel(SubsubFolder)
        try
            path = append(MainFolder, filesep, SubFolder{r}, filesep, SubsubFolder{o});
    
            %% Loading
            folder = dir(path);
            idx = contains({folder.name},'trackResultsCommonCh.mat');
            folder(~idx) = [];
            f2Load = [folder(1).folder filesep folder(1).name];
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});
            
            calibration = load(append(path2RotCal, filesep, 'RotCalib.mat'));
            name = fieldnames(calibration);
            calibration = calibration.(name{1,1});
            %% Processing
            currMov =  data(1).traces;
            currMov = table2cell(currMov);
            allHeight = cellfun(@height,currMov(:,1));
            % allHeight = size(currMov,1);
            idx = allHeight>minSize;
            currMov = currMov(idx,:);
            allRes = struct('msadTheta',0,'msadPhi',0,'msadr',0,'tau',0,'DTheta',0,'DPhi',0,'Dr',0,...
                'nTheta',0,'nPhi',0,'nr',0,'vTheta',0,'vPhi',0,'vr',0, 'Totamp', 0);
            allRes(length(currMov)).masdTheta = [];
            maxLength = max(allHeight);
            allmsadTheta = zeros(length(currMov),maxLength-1);
            allmsadPhi = allmsadTheta;
            allmsad = allmsadPhi;
            
            for i = 1:size(currMov,1)
                try
                    waitbar(i./size(currMov,1),f, append('doing microrheology: sample ', num2str(o + (r - 1)*(numel(SubsubFolder))), ' out of ', num2str(numel(SubFolder)*numel(SubsubFolder))));
                    currPart = currMov(i,:);
                
                    %%% calculate angels
                    TotInt = currPart{1,3} + currPart{1,4};
                    I1 = currPart{1,3}./TotInt;
                    I2 = currPart{1,4}./TotInt;
                    Diff = I1 - I2;
                    Time = currPart{1,5};

                    Phi = 0.5*real(acos(sqrt(TotInt/(calibration.TotI_mean))));
                    Theta = 0.25*real(acos(Diff./calibration.I_mean));          
                    coord = [Theta, Phi];

                    % For Theta
                    tau = Time;
                    [msadTheta, tau] = MSD.Rotational.calc2(coord(:,1), tau, expTime);
                    DTheta   = MSD.Rotational.getDiffCoeff(msadTheta,tau,fitRDiff,'2D');
                    nTheta   = MSD.Rotational.getViscosity(DTheta,R,ParticleType, Temp);
                    vTheta   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi;

                    DiffImean = msadTheta(1,1);


                    allRes(i).msadTheta = msadTheta;% in rad^2 
                    allRes(i).tau = tau; % in sec

                    allRes(i).DTheta   = DTheta;% in rad^2 /sec
                    allRes(i).nTheta   = nTheta;
                    allRes(i).vTheta   = vTheta;
                    allRes(i).num  = length(msadTheta);
                    allRes(i).diffI = DiffImean;
                    Visc(end+1, 1) = abs(nTheta);


                catch
                end
                    
            end

            disp(append('Diff = ', num2str(nanmean([allRes.DTheta]))))
            disp(append('Visc = ', num2str(nanmean([allRes.nTheta]))))

            DResultsTheta(o, r) = nanmedian([allRes.DTheta]);
            nResultsTheta(o, r) = nanmedian([allRes.nTheta]);
        catch
            DResultsTheta(o, r) = nan;
            nResultsTheta(o, r) = nan;
        end
    end
  
    disp(append('The viscosity is ', num2str(nanmedian(Visc)), ' cP'))
    disp(append('number of traces is ', num2str(numel(Visc))))
end
close(f)