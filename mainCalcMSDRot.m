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
path = 'S:\Rotational Tracking\20250303_AuBPS_184x92_glycerol\AuBPs_184x92_in_glycerol\3_cP\sample5';
path2RotCal = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';

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
    'nTheta',0,'nPhi',0,'nr',0,'vTheta',0,'vPhi',0,'vr',0);
allRes(length(currMov)).masdTheta = [];
maxLength = max(allHeight);
allmsadTheta = zeros(length(currMov),maxLength-1);
allmsadPhi = allmsadTheta;
allmsad = allmsadPhi;

for i = 1:size(currMov,1)
    currPart = currMov(i,:);

    %%% calculate angels
    TotInt = currPart{1,1} + currPart{1,2};
    I1 = currPart{1,1}./TotInt;
    I2 = currPart{1,2}./TotInt;
    Diff = I1 - I2;

    Theta = 0.5*real(acos(Diff./calibration.I0_mean));
    Phi = real(acos(sqrt(TotInt./calibration.TotI0_mean)));

    coord = [Theta, Phi];

    %For Theta
    msadTheta = MSD.Rotational.calc(coord(:,1));
    tau = (1:length(msadTheta))'*expTime;
    allmsadTheta(i,1:length(msadTheta)) = msadTheta;
    DTheta   = MSD.Rotational.getDiffCoeff(msadTheta,tau,fitRDiff,'2D');
    nTheta   = MSD.Rotational.getViscosity(DTheta,R,ParticleType, Temp);
    vTheta   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s

    %For Phi
    msadPhi = MSD.Rotational.calc(coord(:,2));
    tau = (1:length(msadPhi))'*expTime;
    allmsadPhi(i,1:length(msadPhi)) = msadPhi;
    DPhi   = MSD.Rotational.getDiffCoeff(msadPhi,tau,fitRDiff,'2D');
    nPhi   = MSD.Rotational.getViscosity(DPhi,R,ParticleType, Temp);
    vPhi   = coord(1,1) - coord(end,1)/(length(coord)*expTime)*180/pi; %degrees/s

    %For both
    msadr = MSD.calc(coord);%convert to um;
    allMSDR(i,1:length(msadr)) = msadr;
    DR   = MSD.Rotational.getDiffCoeff(msadr,tau,fitRDiff,'3D');
    nR   = MSD.Rotational.getViscosity(DR,R,ParticleType,Temp);
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
DR   = MSD.Rotational.getDiffCoeff(meanMSDR,tau,fitRDiff,'3D');
nR   = MSD.Rotational.getViscosity(DR,R,ParticleType,Temp);
disp(['The diffusion coefficient is ' num2str(DR) 'µm^2s^-^1' '\n'...
    'the viscosity is ' num2str(nR) ' cp' '\n']);
%%
filename = [path filesep 'msadRes.mat'];
save(filename,'allRes');
h = msgbox('Data succesfully saved');