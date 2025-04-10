clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.010; %in sec
Temp = 296.15; %temperature in Kelvin
ParticleType = 'Bipyramid'; %bipyramid, ellipsoid, rod, cilinder,...
R = [184, 92]; %Long axis, short axis in nm
fitRDiff = 4; %in number of data
minSize = 20; %frames
ext = '.mat';
path2RotCal = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';

%% Path info
MainFolder = 'S:\Rotational Tracking\20250407_AuBPs_184s92_glycerol\Glycerol';
SubFolder = {'glycerol_80', 'glycerol_85', 'glycerol_90','glycerol_95', 'glycerol_100'}; % 'glycerol_80', 'glycerol_85', 'glycerol_90','glycerol_95', 
SubsubFolder = {'sample1', 'sample2', 'sample3','sample4', 'sample5'}; %


f = waitbar(0,'Initializing');
DiffIMatrix = [];
DiffTotMatrix = [];
for r = 1:numel(SubFolder)
    Visc = [];
    DiffTotCol = [];
    DiffICol = [];
    stepsizes = [];
    for o = 1:numel(SubsubFolder)
        % try
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
            allHeight = cellfun(@height,currMov(:,1));
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
                    TotInt = currPart{1,1} + currPart{1,2};
                    I1 = currPart{1,1}./TotInt;
                    I2 = currPart{1,2}./TotInt;
                    Diff = I1 - I2;
                    Time = currPart{1,5};

                    Phi = real(acos(sqrt(TotInt./calibration.TotI0_mean)));
                    ampI = calibration.I0_mean*(cos(Phi)).^2;
                
                    Theta = 0.5*real(acos(Diff./calibration.I0_mean));          
                    coord = [Theta, Phi];

                    %  DiffTot = diff(Phi);
                    % DiffI = diff(Theta);
                    % ik = abs(diff(Time)-0.01) < 0.0002; 
                    % DiffTotmean = abs(mean(DiffTot(ik)));
                    % DiffImean = abs(mean(DiffI(ik)));

                    % For Theta
                    tau = Time;
                    [msadTheta] = MSD.Rotational.calc(coord(:,1), tau, expTime);
                    allmsadTheta(i,1:length(msadTheta)) = msadTheta(1,:);
                    DTheta   = MSD.Rotational.getDiffCoeff(msadTheta,tau,fitRDiff,'2D');
                    nTheta   = MSD.Rotational.getViscosity(DTheta,R,ParticleType, Temp);
                    vTheta   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s
    
                    % %For Phi
                    tau = Time;
                    msadPhi = MSD.Rotational.calc(coord(:,2), tau, expTime);
                    tau = (1:length(msadPhi))'*expTime;
                    allmsadPhi(i,1:length(msadPhi)) = msadPhi(1,:);
                    DPhi   = MSD.Rotational.getDiffCoeff(msadPhi,tau,fitRDiff,'2D');
                    nPhi   = MSD.Rotational.getViscosity(DPhi,R,ParticleType, Temp);
                    vPhi   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s
                    
                    %For both
                    tau = Time;
                    msadr = MSD.Rotational.calc(coord, tau, expTime);%convert to um;
                    allMSDR(i,1:length(msadr)) = msadr(1,:);
                    DR   = MSD.Rotational.getDiffCoeff(msadr,tau,fitRDiff,'3D');
                    nR   = MSD.Rotational.getViscosity(DR,R,ParticleType,Temp);
                    dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2);
                    vR = dR/10^3/(length(coord)*expTime); %rad/s

                   
                    DiffTotmean = msadPhi(1,1);
                    DiffImean = msadTheta(1,1);
                    
                
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
                    
                    allRes(i).diffTot = DiffTotmean;
                    allRes(i).diffI = DiffImean;
                    Visc(end+1, 1) = abs(nR);
                catch
                end
              
            end
            
            filename = [path filesep 'msadRes.mat'];
            save(filename,'allRes');
            disp(append('Phi visc = ', num2str(median([allRes.nPhi]))))
            disp(append('Theta visc = ', num2str(median([allRes.nTheta]))))
            disp(append('Total visc = ', num2str(median([allRes.nr]))))
            disp(append('stepsize= ', num2str(median([stepsizes]))))
            % disp(append('tot max = ', num2str(median([allRes.Dr]))))
            % disp(append('tot min = ', num2str(median([allRes.totMin]))))

            DiffTotCol = [DiffTotCol; mean([allRes.DPhi])];
            DiffICol = [DiffICol; mean([allRes.DTheta])];
    end
   
    DiffTotMatrix = [DiffTotMatrix, DiffTotCol];
    DiffIMatrix = [DiffIMatrix, DiffICol];
    disp(append('The viscosity is ', num2str(nanmean(Visc)), ' cP'))
    disp(append('number of traces is ', num2str(numel(Visc))))
end
 close(f)