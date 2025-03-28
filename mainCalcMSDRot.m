clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.010; %in sec
Temp = 296.15; %temperature in Kelvin
AuBP = [184, 92]; %Long axis, short axis in nm
fitRDiff = 4; %in number of data
minSize = 15; %frames
ext = '.mat';
path = 'S:\Rotational Tracking\20250303_AuBPS_184x92_glycerol\AuBPs_184x92_in_glycerol\3_cP\sample1';
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
allRes = struct('masdTheta',0,'masdPhi',0,'masd',0,'tau',0,'DrTheta',0,'DrPhi',0,'Dr',0,...
    'nrTheta',0,'nrPhi',0,'nr',0,'vrTheta',0,'vrPhi',0,'vr',0);
allRes(length(currMov)).masdTheta = [];
maxLength = max(allHeight);
allmasdTheta = zeros(length(currMov),maxLength-1);
allmasdPhi = allmasdTheta;
allmasd = allmasdPhi;

for i = 1:length(currMov)
    currPart = currMov(i,:);

    %%% calculate angels
    TotInt = currPart{1,1} + currPart{1,2};
    I1 = currPart{1,1}./TotInt;
    I2 = currPart{1,2}./TotInt;
    Diff = currPart{1,1} + currPart{1,2};

    Theta = 0.5*real(acos(Diff./calibration.I0_mean));
    Phi = real(acos(sqrt(TotInt./calibration.TotI0_mean)));

    coord = [Theta, Phi];

    %For Theta
    masdTheta = MSD.Rotational.calc(coord(:,1));
    tau = (1:length(msdx))'*expTime;
    allMSDX(i,1:length(msdx)) = msdx;
    DX   = MSD.getDiffCoeff(msdx,tau,fitRDiff,'1D');
    nX   = MSD.getViscosity(DX,R,Temp);
    aX   = MSD.getDiffTypeAlpha(msdx,expTime);
    vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s

    %For Phi
    msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
    allMSDY(i,1:length(msdy)) = msdy;
    DY   = MSD.getDiffCoeff(msdy,tau,fitRDiff,'1D');
    nY   = MSD.getViscosity(DY,R,Temp);
    aY   = MSD.getDiffTypeAlpha(msdy,expTime);
    vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s

    %For both
    msdr = MSD.calc(coord/10^3);%convert to um;
    allMSDR(i,1:length(msdr)) = msdr;
    DR   = MSD.getDiffCoeff(msdr,tau,fitRDiff,'3D');
    nR   = MSD.getViscosity(DR,R,Temp);
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
end

%%

meanMSDR = nanmean(allMSDR,1);
tau = (1:length(meanMSDR))'*expTime;
DR   = MSD.getDiffCoeff(meanMSDR,tau,fitRDiff,'3D');
nR   = MSD.getViscosity(DR,R,Temp);
nRg = nanmean(allRg);
disp(['The diffusion coefficient is ' num2str(DR) 'µm^2s^-^1' '\n'...
    'the viscosity is ' num2str(nR) ' cp' '\n'...
    'the radius of gyration is ' num2str(nRg) ' µm']);
%%
filename = [path filesep 'msdRes1.mat'];
save(filename,'allRes');
h = msgbox('Data succesfully saved');