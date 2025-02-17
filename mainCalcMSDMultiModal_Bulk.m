clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.0305; %in sec
Temp = 296.15; %temperature in Kelvin
<<<<<<< HEAD
R1 = 0.050; %Radius of particle in um for channel 1 (planes 1-8);
=======
<<<<<<< HEAD
R1 = 0.11345; %Radius of particle in um for channel 1 (planes 1-8);
R2 = 0.11415; %Radius of particle in um for channel 1 (planes 9-16);
=======
R1 = 0.100; %Radius of particle in um for channel 1 (planes 1-8);
>>>>>>> 08411330cc35a0f416504520e8c1f4351498a589
R2 = 0.100; %Radius of particle in um for channel 1 (planes 9-16);
>>>>>>> ce14b8dba61579ab08ba8ab6e1098281c0fbcf9a
MultiModal = 'MultiColor'; %MultiColor or Rotational tracking
fitRDiff = 4; %in number of data
minSize = 30; %frames
ext = '.mat';
<<<<<<< HEAD
MainFolder = 'S:\DDM_TestData';
SubFolders = {'PS_100_low_conc', 'PS_200_low_conc'};
SubsubFolders = {'sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'};
=======
<<<<<<< HEAD
MainFolder = 'D:\Dual Color';
SubFolders = {'20250121', '20250122'};
SubsubFolders = {'PS_200_green_PS_100_red'}; %'Multicolor_particles','PS_200_green_PS_100_red' 
SubsubsubFolders = {'0_min_measurements', 'sample1', 'sample2', 'sample3'};
SubsubsubsubFolders = {'0_min1', '0_min2', '0_min3', '0_min4', '0_min5', ...
'3_min', '5_min', '7_min', '9_min', '11_min', '13_min', '15_min', '17_min'};
=======
MainFolder = 'S:\Dual Color\20250122\Multicolor_particles';
SubFolders = {'In_water'};
SubsubFolders = {'0_min1', '0_min2', '0_min3', '0_min4', '0_min6'};
>>>>>>> ce14b8dba61579ab08ba8ab6e1098281c0fbcf9a
>>>>>>> 08411330cc35a0f416504520e8c1f4351498a589


gn0 = [];
gD0 = [];
ga0 = [];
gRg0 = [];
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
gn17 = [];
gD17 = [];
ga17 = [];
gRg17 = [];
rn0 = [];
rD0 = [];
ra0 = [];
rRg0 = [];
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
rn17 = [];
rD17 = [];
ra17 = [];
rRg17 = [];

<<<<<<< HEAD
for m = 1:numel(SubFolders)
    SubFolder = SubFolders{m};
    for u = 1:numel(SubsubFolders)
        SubsubFolder = SubsubFolders{u};
        path = append(MainFolder, filesep, SubFolder, filesep, SubsubFolder);
        % try
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
            %currMov =  data(1).traces;
            currMov = data.traces{1,1};
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
            
            
            % %% Loading
            % R = R2;
            % folder = dir(path);
            % idx = contains({folder.name},'trackResults2.mat');
            % folder(~idx) = [];
            % 
            % f2Load = [folder(1).folder filesep folder(1).name];
            % 
            % tmpData = load(f2Load);
            % name = fieldnames(tmpData);
            % data = tmpData.(name{1});
            % 
            % %% Processing
            % currMov =  data(1).traces;
            % allHeight = cellfun(@height,currMov(:,1));
            % idx = allHeight>minSize;
            % currMov = currMov(idx,1);
            % allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
            %     'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
            %     'vZ',0,'vR',0);
            % allRes(length(currMov)).msdX = [];
            % maxLength = max(allHeight);
            % allMSDX = zeros(length(currMov),maxLength-1);
            % allMSDY = allMSDX;
            % allMSDZ = allMSDY;
            % allMSDR = allMSDY;
            % allRg = zeros(length(currMov),1);
            % for i = 1:length(currMov)
            %     currPart = currMov{i};
            % 
            %     coord = [currPart.col, currPart.row, currPart.z];
            %     CM = mean(coord,1);
            %     coord = coord-CM;
            % 
            %     %in X
            %     msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
            %     tau = (1:length(msdx))'*expTime;
            %     allMSDX(i,1:length(msdx)) = msdx;
            %     DX   = MSD.getDiffCoeff(msdx,tau,fitRDiff,'1D');
            %     nX   = MSD.getViscosity(DX,R,Temp);
            %     aX   = MSD.getDiffTypeAlpha(msdx,expTime);
            %     vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s
            % 
            %     %inY
            %     msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
            %     allMSDY(i,1:length(msdy)) = msdy;
            %     DY   = MSD.getDiffCoeff(msdy,tau,fitRDiff,'1D');
            %     nY   = MSD.getViscosity(DY,R,Temp);
            %     aY   = MSD.getDiffTypeAlpha(msdy,expTime);
            %     vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s
            % 
            %     %inZ
            %     msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
            %     allMSDZ(i,1:length(msdz)) = msdz;
            %     DZ   = MSD.getDiffCoeff(msdz,tau,fitRDiff,'1D');
            %     nZ   = MSD.getViscosity(DZ,R,Temp);
            %     aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
            %     vZ   = coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime); %um/s
            % 
            % 
            %     %inR
            %     msdr = MSD.calc(coord/10^3);%convert to um;
            %     allMSDR(i,1:length(msdr)) = msdr;
            %     DR   = MSD.getDiffCoeff(msdr,tau,fitRDiff,'3D');
            %     nR   = MSD.getViscosity(DR,R,Temp);
            %     aR   = MSD.getDiffTypeAlpha(msdr,expTime);
            %     dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2 +...
            %         (coord(1,3)-coord(end,3))^2);
            %     vR = dR/10^3/(length(coord)*expTime); %um/s
            % 
            % 
            %     allRes(i).msdX = msdx;% in um^2
            %     allRes(i).msdY = msdy;
            %     allRes(i).msdZ = msdz;
            %     allRes(i).msdR = msdr;
            %     allRes(i).tau = tau; % in sec
            % 
            % 
            %     allRes(i).DX   = DX;% in um^2 /sec
            %     allRes(i).DY   = DY;% in um^2 /sec
            %     allRes(i).DZ   = DZ;% in um^2 /sec
            %     allRes(i).DR   = DR;% in um^2 /sec
            % 
            %     allRes(i).nX   = nX;
            %     allRes(i).nY   = nY;
            %     allRes(i).nZ   = nZ;
            %     allRes(i).nR   = nR;
            % 
            %     allRes(i).aX   = aX;
            %     allRes(i).aY   = aY;
            %     allRes(i).aZ   = aZ;
            %     allRes(i).aR   = aR;
            % 
            %     allRes(i).vX   = vX;
            %     allRes(i).vY   = vY;
            %     allRes(i).vZ   = vZ;
            %     allRes(i).vR   = vR;
            % 
            %     %%% Calculate Rg
            %     Coord = coord./(10^3);
            %     N = size(Coord, 1);
            % 
            %     avX = sum(Coord(:,1))./N;
            %     avY = sum(Coord(:,2))./N;
            %     avZ = sum(Coord(:,3))./N;
            % 
            %     %substract averages from coordinates
            %     Coord(:,1) = Coord(:,1) - avX;
            %     Coord(:,2) = Coord(:,2) - avY;
            %     Coord(:,3) = Coord(:,3) - avZ;
            % 
            %     %build gyration tensor for each trajectory
            %     T(1,1) = (1/N)*sum((Coord(:,1)).^2);
            %     T(1,2) = (1/N)*sum((Coord(:,1)).*(Coord(:,2)));
            %     T(1,3) = (1/N)*sum((Coord(:,1)).*(Coord(:,3)));
            %     T(2,1) = (1/N)*sum((Coord(:,2)).*(Coord(:,1)));
            %     T(2,2) = (1/N)*sum((Coord(:,2)).^2);
            %     T(2,3) = (1/N)*sum((Coord(:,2)).*(Coord(:,3)));
            %     T(3,1) = (1/N)*sum((Coord(:,3)).*(Coord(:,1)));
            %     T(3,2) = (1/N)*sum((Coord(:,3)).*(Coord(:,2)));
            %     T(3,3) = (1/N)*sum((Coord(:,3)).^2);
            % 
            %     e = eig(T);
            %     Rg = sqrt((e(1).^2) + e(2).^2 + e(3).^2);
            % 
            %     allRes(i).Rg = Rg;
            %     allRg(i,1) = Rg;
            % 
            % 
            %     allRes(i).num  = length(msdx);
            % end
            % 
            % %%
            % 
            % meanMSDR = nanmean(allMSDR,1);
            % tau = (1:length(meanMSDR))'*expTime;
            % DR   = MSD.getDiffCoeff(meanMSDR,tau,fitRDiff,'3D');
            % nR   = MSD.getViscosity(DR,R,Temp);
            % nRg = nanmean(allRg);
            % disp(['The diffusion coefficient is ' num2str(DR) 'µm^2s^-^1' '\n'...
            %     'the viscosity is ' num2str(nR) ' cp' '\n'...
            %     'the radius of gyration is ' num2str(nRg) ' µm']);
            % %%
            % filename = [path filesep 'msdRes2.mat'];
            % save(filename,'allRes');
            % h = msgbox('Data succesfully saved');
            % 
            % for h = 1:size(allRes,2)
            %     if strcmp(SubsubFolder, '3 min')
            %         gn3(end+1, 1) = allRes(h).nR;
            %         gD3(end+1, 1) = allRes(h).DR;
            %         ga3(end+1, 1) = allRes(h).aR;
            %         gRg3(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '5 min')
            %         gn5(end+1, 1) = allRes(h).nR;
            %         gD5(end+1, 1) = allRes(h).DR;
            %         ga5(end+1, 1) = allRes(h).aR;
            %         gRg5(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '7 min')
            %         gn7(end+1, 1) = allRes(h).nR;
            %         gD7(end+1, 1) = allRes(h).DR;
            %         ga7(end+1, 1) = allRes(h).aR;
            %         gRg7(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '9 min')
            %         gn9(end+1, 1) = allRes(h).nR;
            %         gD9(end+1, 1) = allRes(h).DR;
            %         ga9(end+1, 1) = allRes(h).aR;
            %         gRg9(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '11 min')
            %         gn11(end+1, 1) = allRes(h).nR;
            %         gD11(end+1, 1) = allRes(h).DR;
            %         ga11(end+1, 1) = allRes(h).aR;
            %         gRg11(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '13 min')
            %         gn13(end+1, 1) = allRes(h).nR;
            %         gD13(end+1, 1) = allRes(h).DR;
            %         ga13(end+1, 1) = allRes(h).aR;
            %         gRg13(end+1, 1) = allRes(h).Rg;
            %     elseif strcmp(SubsubFolder, '15 min')
            %         gn15(end+1, 1) = allRes(h).nR;
            %         gD15(end+1, 1) = allRes(h).DR;
            %         ga15(end+1, 1) = allRes(h).aR;
            %         gRg15(end+1, 1) = allRes(h).Rg;
            %     end
            % end
        % catch
        %     continue
        % end
=======
for t = 1:numel(SubFolders)
    for r = 1:numel(SubsubFolders)
        for a = 1:numel(SubsubsubFolders)
            for c = 1:numel(SubsubsubsubFolders)
                path = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r},...
                    filesep, SubsubsubFolders{a}, filesep, SubsubsubsubFolders{c});
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
                        if strcmp(SubsubsubsubFolders{c}, '0_min1')
                            gn0(end+1, 1) = allRes(h).nR;
                            gD0(end+1, 1) = allRes(h).DR;
                            ga0(end+1, 1) = allRes(h).aR;
                            gRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0_min2')
                            gn0(end+1, 1) = allRes(h).nR;
                            gD0(end+1, 1) = allRes(h).DR;
                            ga0(end+1, 1) = allRes(h).aR;
                            gRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0_min3')
                            gn0(end+1, 1) = allRes(h).nR;
                            gD0(end+1, 1) = allRes(h).DR;
                            ga0(end+1, 1) = allRes(h).aR;
                            gRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0_min4')
                            gn0(end+1, 1) = allRes(h).nR;
                            gD0(end+1, 1) = allRes(h).DR;
                            ga0(end+1, 1) = allRes(h).aR;
                            gRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0_min5')
                            gn0(end+1, 1) = allRes(h).nR;
                            gD0(end+1, 1) = allRes(h).DR;
                            ga0(end+1, 1) = allRes(h).aR;
                            gRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '3_min')
                            gn3(end+1, 1) = allRes(h).nR;
                            gD3(end+1, 1) = allRes(h).DR;
                            ga3(end+1, 1) = allRes(h).aR;
                            gRg3(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '5_min')
                            gn5(end+1, 1) = allRes(h).nR;
                            gD5(end+1, 1) = allRes(h).DR;
                            ga5(end+1, 1) = allRes(h).aR;
                            gRg5(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '7_min')
                            gn7(end+1, 1) = allRes(h).nR;
                            gD7(end+1, 1) = allRes(h).DR;
                            ga7(end+1, 1) = allRes(h).aR;
                            gRg7(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '9_min')
                            gn9(end+1, 1) = allRes(h).nR;
                            gD9(end+1, 1) = allRes(h).DR;
                            ga9(end+1, 1) = allRes(h).aR;
                            gRg9(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '11_min')
                            gn11(end+1, 1) = allRes(h).nR;
                            gD11(end+1, 1) = allRes(h).DR;
                            ga11(end+1, 1) = allRes(h).aR;
                            gRg11(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '13_min')
                            gn13(end+1, 1) = allRes(h).nR;
                            gD13(end+1, 1) = allRes(h).DR;
                            ga13(end+1, 1) = allRes(h).aR;
                            gRg13(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '15_min')
                            gn15(end+1, 1) = allRes(h).nR;
                            gD15(end+1, 1) = allRes(h).DR;
                            ga15(end+1, 1) = allRes(h).aR;
                            gRg15(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '17_min')
                            gn17(end+1, 1) = allRes(h).nR;
                            gD17(end+1, 1) = allRes(h).DR;
                            ga17(end+1, 1) = allRes(h).aR;
                            gRg17(end+1, 1) = allRes(h).Rg;
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
                        if strcmp(SubsubsubsubFolders{c}, '0min_1')
                            rn0(end+1, 1) = allRes(h).nR;
                            rD0(end+1, 1) = allRes(h).DR;
                            ra0(end+1, 1) = allRes(h).aR;
                            rRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0min_2')
                            rn0(end+1, 1) = allRes(h).nR;
                            rD0(end+1, 1) = allRes(h).DR;
                            ra0(end+1, 1) = allRes(h).aR;
                            rRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0min_3')
                            rn0(end+1, 1) = allRes(h).nR;
                            rD0(end+1, 1) = allRes(h).DR;
                            ra0(end+1, 1) = allRes(h).aR;
                            rRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0min_4')
                            rn0(end+1, 1) = allRes(h).nR;
                            rD0(end+1, 1) = allRes(h).DR;
                            ra0(end+1, 1) = allRes(h).aR;
                            rRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '0min_5')
                            rn0(end+1, 1) = allRes(h).nR;
                            rD0(end+1, 1) = allRes(h).DR;
                            ra0(end+1, 1) = allRes(h).aR;
                            rRg0(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '3_min')
                            rn3(end+1, 1) = allRes(h).nR;
                            rD3(end+1, 1) = allRes(h).DR;
                            ra3(end+1, 1) = allRes(h).aR;
                            rRg3(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '5_min')
                            rn5(end+1, 1) = allRes(h).nR;
                            rD5(end+1, 1) = allRes(h).DR;
                            ra5(end+1, 1) = allRes(h).aR;
                            rRg5(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '7_min')
                            rn7(end+1, 1) = allRes(h).nR;
                            rD7(end+1, 1) = allRes(h).DR;
                            ra7(end+1, 1) = allRes(h).aR;
                            rRg7(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '9_min')
                            rn9(end+1, 1) = allRes(h).nR;
                            rD9(end+1, 1) = allRes(h).DR;
                            ra9(end+1, 1) = allRes(h).aR;
                            rRg9(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '11_min')
                            rn11(end+1, 1) = allRes(h).nR;
                            rD11(end+1, 1) = allRes(h).DR;
                            ra11(end+1, 1) = allRes(h).aR;
                            rRg11(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '13_min')
                            rn13(end+1, 1) = allRes(h).nR;
                            rD13(end+1, 1) = allRes(h).DR;
                            ra13(end+1, 1) = allRes(h).aR;
                            rRg13(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '15_min')
                            rn15(end+1, 1) = allRes(h).nR;
                            rD15(end+1, 1) = allRes(h).DR;
                            ra15(end+1, 1) = allRes(h).aR;
                            rRg15(end+1, 1) = allRes(h).Rg;
                        elseif strcmp(SubsubsubsubFolders{c}, '17_min')
                            rn17(end+1, 1) = allRes(h).nR;
                            rD17(end+1, 1) = allRes(h).DR;
                            ra17(end+1, 1) = allRes(h).aR;
                            rRg17(end+1, 1) = allRes(h).Rg;
                        end
                    end
                catch
                    continue
                end
            end
        end

        ExcelNameG = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r}, filesep, 'ResultsGreen.xlsx');
        writematrix(gD0, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'A2:A5000');
        writematrix(gD3, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'B2:B5000');
        writematrix(gD5, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'C2:C5000');
        writematrix(gD7, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'D2:D5000');
        writematrix(gD9, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'E2:E5000');
        writematrix(gD11, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'F2:F5000');
        writematrix(gD13, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'G2:GA5000');
        writematrix(gD15, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'H2:H5000');
        writematrix(gD17, ExcelNameG, 'Sheet', 'Diffusion', 'Range', 'I2:I5000');
        
        writematrix(gn0, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'A2:A5000');
        writematrix(gn3, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'B2:B5000');
        writematrix(gn5, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'C2:C5000');
        writematrix(gn7, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'D2:D5000');
        writematrix(gn9, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'E2:E5000');
        writematrix(gn11, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'F2:F5000');
        writematrix(gn13, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'G2:GA5000');
        writematrix(gn15, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'H2:H5000');
        writematrix(gn17, ExcelNameG, 'Sheet', 'Viscosity', 'Range', 'I2:I5000');
        
        writematrix(ga0, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'A2:A5000');
        writematrix(ga3, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'B2:B5000');
        writematrix(ga5, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'C2:C5000');
        writematrix(ga7, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'D2:D5000');
        writematrix(ga9, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'E2:E5000');
        writematrix(ga11, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'F2:F5000');
        writematrix(ga13, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'G2:GA5000');
        writematrix(ga15, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'H2:H5000');
        writematrix(ga17, ExcelNameG, 'Sheet', 'AnExp', 'Range', 'I2:I5000');
        
        writematrix(gRg0, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'A2:A5000');
        writematrix(gRg3, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'B2:B5000');
        writematrix(gRg5, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'C2:C5000');
        writematrix(gRg7, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'D2:D5000');
        writematrix(gRg9, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'E2:E5000');
        writematrix(gRg11, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'F2:F5000');
        writematrix(gRg13, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'G2:GA5000');
        writematrix(gRg15, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'H2:H5000');
        writematrix(gRg17, ExcelNameG, 'Sheet', 'Radius of Gyration', 'Range', 'I2:I5000');

        ExcelNameR = append(MainFolder, filesep, SubFolders{t}, filesep, SubsubFolders{r}, filesep, 'ResultsRed.xlsx');
        writematrix(rD0, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'A2:A5000');
        writematrix(rD3, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'B2:B5000');
        writematrix(rD5, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'C2:C5000');
        writematrix(rD7, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'D2:D5000');
        writematrix(rD9, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'E2:E5000');
        writematrix(rD11, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'F2:F5000');
        writematrix(rD13, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'G2:GA5000');
        writematrix(rD15, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'H2:H5000');
        writematrix(rD17, ExcelNameR, 'Sheet', 'Diffusion', 'Range', 'I2:I5000');
       
        writematrix(rn0, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'A2:A5000');
        writematrix(rn3, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'B2:B5000');
        writematrix(rn5, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'C2:C5000');
        writematrix(rn7, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'D2:D5000');
        writematrix(rn9, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'E2:E5000');
        writematrix(rn11, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'F2:F5000');
        writematrix(rn13, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'G2:GA5000');
        writematrix(rn15, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'H2:H5000');
        writematrix(rn17, ExcelNameR, 'Sheet', 'Viscosity', 'Range', 'I2:I5000');

        writematrix(ra0, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'A2:A5000');
        writematrix(ra3, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'B2:B5000');
        writematrix(ra5, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'C2:C5000');
        writematrix(ra7, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'D2:D5000');
        writematrix(ra9, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'E2:E5000');
        writematrix(ra11, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'F2:F5000');
        writematrix(ra13, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'G2:GA5000');
        writematrix(ra15, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'H2:H5000');
        writematrix(ra17, ExcelNameR, 'Sheet', 'AnExp', 'Range', 'I2:I5000');
        
        writematrix(gRg0, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'A2:A5000');
        writematrix(gRg3, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'B2:B5000');
        writematrix(gRg5, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'C2:C5000');
        writematrix(gRg7, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'D2:D5000');
        writematrix(gRg9, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'E2:E5000');
        writematrix(gRg11, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'F2:F5000');
        writematrix(gRg13, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'G2:GA5000');
        writematrix(gRg15, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'H2:H5000');
        writematrix(gRg17, ExcelNameR, 'Sheet', 'Radius of Gyration', 'Range', 'I2:I5000');
        
>>>>>>> 08411330cc35a0f416504520e8c1f4351498a589
    end
end