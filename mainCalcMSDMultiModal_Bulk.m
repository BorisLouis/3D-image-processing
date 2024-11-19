clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.010; %in sec
Temp = 296.15; %temperature in Kelvin
R1 = 0.100; %Radius of particle in um for channel 1 (planes 1-8);
R2 = 0.050; %Radius of particle in um for channel 1 (planes 9-16);
MultiModal = 'MultiColor'; %MultiColor or Rotational tracking
fitRDiff = 4; %in number of data
minSize = 20; %frames
ext = '.mat';
MainFolder = 'G:\multicolor_polarization\multicolor\20241105_polymerisation_dual_color_PS_air_obj';
SubFolders = {'sample1', 'sample2', 'sample3'};
SubsubFolders = {'3 min', '5 min', '7 min', '9 min', '11 min', '13 min', '15 min'};

gn3 = [];
gD3 = [];
ga3 = [];
gRg3 = [];
gn5 = [];
gD5 = [];
ga5 = [];
gRg5 = [];
gn7 = [];
gD7 = [];
ga7 = [];
gRg7 = [];
gn9 = [];
gD9 = [];
ga9 = [];
gRg9 = [];
gn11 = [];
gD11  = [];
ga11 = [];
gRg11 = [];
gn13 = [];
gD13 = [];
ga13 = [];
gRg13 = [];
gn15 = [];
gD15 = [];
ga15 = [];
gRg15 = [];
rn3 = [];
rD3 = [];
ra3 = [];
rRg3 = [];
rn5 = [];
rD5 = [];
ra5 = [];
rRg5 = [];
rn7 = [];
rD7 = [];
ra7 = [];
rRg7 = [];
rn9 = [];
rD9 = [];
ra9 = [];
rRg9 = [];
rn11 = [];
rD11  = [];
ra11 = [];
rRg11 = [];
rn13 = [];
rD13 = [];
ra13 = [];
rRg13 = [];
rn15 = [];
rD15 = [];
ra15 = [];
rRg15 = [];

for m = 1:numel(SubFolders)
    SubFolder = SubFolders{m};
    for u = 1:numel(SubsubFolders)
        SubsubFolder = SubsubFolders{u};
        path = append(MainFolder, filesep, SubFolder, filesep, SubsubFolder);
        try
            %% Loading
            R = R1;
            folder = dir(path);
            idx = contains({folder.name},'trackResults1.mat');
            folder(~idx) = [];
            
            f2Load = [folder(1).folder filesep folder(1).name];
            
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});
            
            %% Processing
            currMov =  data(1).traces;
            allHeight = cellfun(@height,currMov(:,1));
            idx = allHeight>minSize;
            currMov = currMov(idx,1);
            allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
                'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
                'vZ',0,'vR',0);
            allRes(length(currMov)).msdX = [];
            maxLength = max(allHeight);
            allMSDX = zeros(length(currMov),maxLength-1);
            allMSDY = allMSDX;
            allMSDZ = allMSDY;
            allMSDR = allMSDY;
            allRg = zeros(length(currMov),1);
            for i = 1:length(currMov)
                currPart = currMov{i};
            
                coord = [currPart.col, currPart.row, currPart.z];
                CM = mean(coord,1);
                coord = coord-CM;
            
                %in X
                msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
                tau = (1:length(msdx))'*expTime;
                allMSDX(i,1:length(msdx)) = msdx;
                DX   = MSD.getDiffCoeff(msdx,tau,fitRDiff,'1D');
                nX   = MSD.getViscosity(DX,R,Temp);
                aX   = MSD.getDiffTypeAlpha(msdx,expTime);
                vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s
            
                %inY
                msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
                allMSDY(i,1:length(msdy)) = msdy;
                DY   = MSD.getDiffCoeff(msdy,tau,fitRDiff,'1D');
                nY   = MSD.getViscosity(DY,R,Temp);
                aY   = MSD.getDiffTypeAlpha(msdy,expTime);
                vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s
            
                %inZ
                msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
                allMSDZ(i,1:length(msdz)) = msdz;
                DZ   = MSD.getDiffCoeff(msdz,tau,fitRDiff,'1D');
                nZ   = MSD.getViscosity(DZ,R,Temp);
                aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
                vZ   = coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime); %um/s
            
            
                %inR
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
            
                %%% Calculate Rg
                Coord = coord./(10^3);
                N = size(Coord, 1);
            
                avX = sum(Coord(:,1))./N;
                avY = sum(Coord(:,2))./N;
                avZ = sum(Coord(:,3))./N;
            
                %substract averages from coordinates
                Coord(:,1) = Coord(:,1) - avX;
                Coord(:,2) = Coord(:,2) - avY;
                Coord(:,3) = Coord(:,3) - avZ;
            
                %build gyration tensor for each trajectory
                T(1,1) = (1/N)*sum((Coord(:,1)).^2);
                T(1,2) = (1/N)*sum((Coord(:,1)).*(Coord(:,2)));
                T(1,3) = (1/N)*sum((Coord(:,1)).*(Coord(:,3)));
                T(2,1) = (1/N)*sum((Coord(:,2)).*(Coord(:,1)));
                T(2,2) = (1/N)*sum((Coord(:,2)).^2);
                T(2,3) = (1/N)*sum((Coord(:,2)).*(Coord(:,3)));
                T(3,1) = (1/N)*sum((Coord(:,3)).*(Coord(:,1)));
                T(3,2) = (1/N)*sum((Coord(:,3)).*(Coord(:,2)));
                T(3,3) = (1/N)*sum((Coord(:,3)).^2);
            
                e = eig(T);
                Rg = sqrt((e(1).^2) + e(2).^2 + e(3).^2);
            
                allRes(i).Rg = Rg;
                allRg(i,1) = Rg;
                
                
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

            for h = 1:size(allRes,2)
                if strcmp(SubsubFolder, '3 min')
                    rn3(end+1, 1) = allRes(h).nR;
                    rD3(end+1, 1) = allRes(h).DR;
                    ra3(end+1, 1) = allRes(h).aR;
                    rRg3(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '5 min')
                    rn5(end+1, 1) = allRes(h).nR;
                    rD5(end+1, 1) = allRes(h).DR;
                    ra5(end+1, 1) = allRes(h).aR;
                    rRg5(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '7 min')
                    rn7(end+1, 1) = allRes(h).nR;
                    rD7(end+1, 1) = allRes(h).DR;
                    ra7(end+1, 1) = allRes(h).aR;
                    rRg7(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '9 min')
                    rn9(end+1, 1) = allRes(h).nR;
                    rD9(end+1, 1) = allRes(h).DR;
                    ra9(end+1, 1) = allRes(h).aR;
                    rRg9(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '11 min')
                    rn11(end+1, 1) = allRes(h).nR;
                    rD11(end+1, 1) = allRes(h).DR;
                    ra11(end+1, 1) = allRes(h).aR;
                    rRg11(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '13 min')
                    rn13(end+1, 1) = allRes(h).nR;
                    rD13(end+1, 1) = allRes(h).DR;
                    ra13(end+1, 1) = allRes(h).aR;
                    rRg13(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '15 min')
                    rn15(end+1, 1) = allRes(h).nR;
                    rD15(end+1, 1) = allRes(h).DR;
                    ra15(end+1, 1) = allRes(h).aR;
                    rRg15(end+1, 1) = allRes(h).Rg;
                end
            end
            
            
            %% Loading
            R = R2;
            folder = dir(path);
            idx = contains({folder.name},'trackResults2.mat');
            folder(~idx) = [];
            
            f2Load = [folder(1).folder filesep folder(1).name];
            
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});
            
            %% Processing
            currMov =  data(1).traces;
            allHeight = cellfun(@height,currMov(:,1));
            idx = allHeight>minSize;
            currMov = currMov(idx,1);
            allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
                'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
                'vZ',0,'vR',0);
            allRes(length(currMov)).msdX = [];
            maxLength = max(allHeight);
            allMSDX = zeros(length(currMov),maxLength-1);
            allMSDY = allMSDX;
            allMSDZ = allMSDY;
            allMSDR = allMSDY;
            allRg = zeros(length(currMov),1);
            for i = 1:length(currMov)
                currPart = currMov{i};
            
                coord = [currPart.col, currPart.row, currPart.z];
                CM = mean(coord,1);
                coord = coord-CM;
            
                %in X
                msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
                tau = (1:length(msdx))'*expTime;
                allMSDX(i,1:length(msdx)) = msdx;
                DX   = MSD.getDiffCoeff(msdx,tau,fitRDiff,'1D');
                nX   = MSD.getViscosity(DX,R,Temp);
                aX   = MSD.getDiffTypeAlpha(msdx,expTime);
                vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s
            
                %inY
                msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
                allMSDY(i,1:length(msdy)) = msdy;
                DY   = MSD.getDiffCoeff(msdy,tau,fitRDiff,'1D');
                nY   = MSD.getViscosity(DY,R,Temp);
                aY   = MSD.getDiffTypeAlpha(msdy,expTime);
                vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s
            
                %inZ
                msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
                allMSDZ(i,1:length(msdz)) = msdz;
                DZ   = MSD.getDiffCoeff(msdz,tau,fitRDiff,'1D');
                nZ   = MSD.getViscosity(DZ,R,Temp);
                aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
                vZ   = coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime); %um/s
            
            
                %inR
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
            
                %%% Calculate Rg
                Coord = coord./(10^3);
                N = size(Coord, 1);
            
                avX = sum(Coord(:,1))./N;
                avY = sum(Coord(:,2))./N;
                avZ = sum(Coord(:,3))./N;
            
                %substract averages from coordinates
                Coord(:,1) = Coord(:,1) - avX;
                Coord(:,2) = Coord(:,2) - avY;
                Coord(:,3) = Coord(:,3) - avZ;
            
                %build gyration tensor for each trajectory
                T(1,1) = (1/N)*sum((Coord(:,1)).^2);
                T(1,2) = (1/N)*sum((Coord(:,1)).*(Coord(:,2)));
                T(1,3) = (1/N)*sum((Coord(:,1)).*(Coord(:,3)));
                T(2,1) = (1/N)*sum((Coord(:,2)).*(Coord(:,1)));
                T(2,2) = (1/N)*sum((Coord(:,2)).^2);
                T(2,3) = (1/N)*sum((Coord(:,2)).*(Coord(:,3)));
                T(3,1) = (1/N)*sum((Coord(:,3)).*(Coord(:,1)));
                T(3,2) = (1/N)*sum((Coord(:,3)).*(Coord(:,2)));
                T(3,3) = (1/N)*sum((Coord(:,3)).^2);
            
                e = eig(T);
                Rg = sqrt((e(1).^2) + e(2).^2 + e(3).^2);
            
                allRes(i).Rg = Rg;
                allRg(i,1) = Rg;
                
                
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
            filename = [path filesep 'msdRes2.mat'];
            save(filename,'allRes');
            h = msgbox('Data succesfully saved');

            for h = 1:size(allRes,2)
                if strcmp(SubsubFolder, '3 min')
                    gn3(end+1, 1) = allRes(h).nR;
                    gD3(end+1, 1) = allRes(h).DR;
                    ga3(end+1, 1) = allRes(h).aR;
                    gRg3(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '5 min')
                    gn5(end+1, 1) = allRes(h).nR;
                    gD5(end+1, 1) = allRes(h).DR;
                    ga5(end+1, 1) = allRes(h).aR;
                    gRg5(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '7 min')
                    gn7(end+1, 1) = allRes(h).nR;
                    gD7(end+1, 1) = allRes(h).DR;
                    ga7(end+1, 1) = allRes(h).aR;
                    gRg7(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '9 min')
                    gn9(end+1, 1) = allRes(h).nR;
                    gD9(end+1, 1) = allRes(h).DR;
                    ga9(end+1, 1) = allRes(h).aR;
                    gRg9(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '11 min')
                    gn11(end+1, 1) = allRes(h).nR;
                    gD11(end+1, 1) = allRes(h).DR;
                    ga11(end+1, 1) = allRes(h).aR;
                    gRg11(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '13 min')
                    gn13(end+1, 1) = allRes(h).nR;
                    gD13(end+1, 1) = allRes(h).DR;
                    ga13(end+1, 1) = allRes(h).aR;
                    gRg13(end+1, 1) = allRes(h).Rg;
                elseif strcmp(SubsubFolder, '15 min')
                    gn15(end+1, 1) = allRes(h).nR;
                    gD15(end+1, 1) = allRes(h).DR;
                    ga15(end+1, 1) = allRes(h).aR;
                    gRg15(end+1, 1) = allRes(h).Rg;
                end
            end
        catch
            continue
        end
    end
end