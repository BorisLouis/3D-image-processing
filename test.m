MainFolder = 'G:\DDM_testdata\20241112_DDM_2D_validation';
DateFolder = {'PS_100nm', 'PS_200nm', 'PS_500nm', 'PS_1000nm'};
CellFolder = {'sample1', 'sample2', 'sample3'};
Category = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4];

CatTrack1 = [];
CatTrack2 = [];
CatTrack3 = [];
CatTrack4 = [];
CatDDM1 = [];
CatDDM2 = [];
CatDDM3 = [];
CatDDM4 = [];

figure()
CatIdx = 1;
for i = 1:size(DateFolder, 2)
    for j = 1:size(CellFolder,2)
        try
            %% Tracking
            FilenameTrack = append(MainFolder, filesep, DateFolder{i}, filesep, CellFolder{j}, filesep, 'msdRes.mat');
            MSDTrack = load(FilenameTrack);
            MSDTrack = MSDTrack.allRes;
            for z = 1:size(MSDTrack,2)
                Diff(z) = MSDTrack(z).DR;
                Visc(z) = MSDTrack(z).nR;
                Anom(z) = MSDTrack(z).aR;
            end  
            CatNew = [mean(Diff), mean(Visc), mean(Anom)];  
            Cat = Category(CatIdx);    
            if Cat == 1
              CatTrack1 = [CatTrack1; CatNew];
            elseif Cat == 2
              CatTrack2 = [CatTrack2; CatNew];
            elseif Cat == 3
              CatTrack3 = [CatTrack3; CatNew];
            elseif Cat == 4
              CatTrack4 = [CatTrack4; CatNew];
            end    
            

            %% DDM
            FilenameDDM1 = append(MainFolder, filesep, DateFolder{i}, filesep, CellFolder{j}, filesep, 'MSD.mat');
            DDMMSD = load(FilenameDDM1);
            DDMMSD = DDMMSD.MSDAv;

            if size(DDMMSD,2) > 100
                FitRange = round(0.06*size(DDMMSD,2));
            else 
                FitRange = round(0.15*size(DDMMSD,2));
            end
            if FitRange < 4
                FitRange = 4;
            end
            
            try
                model = 'a*x + log(b) + log(4)';
                ft = fittype(model);
                options = fitoptions(ft);
                options.Lower = [0 0];
                options.Upper = [1 10];
                options.StartPoint = [1 0.5];
                [f, gov] = fit(log(DDMMSD(2,2:FitRange)).', log(DDMMSD(1,2:FitRange)).', ft);
                coef = coeffvalues(f);
                DDMD = coef(2);
            catch
                DDMD = NaN;
            end
            
            if Cat == 1
              CatDDM1 = [CatDDM1; DDMD];
              c = 'r';
            elseif Cat == 2
              CatDDM2 = [CatDDM2; DDMD];
              c = 'y';
            elseif Cat == 3
              CatDDM3 = [CatDDM3; DDMD];
              c = 'g';
           elseif Cat == 4
              CatDDM4 = [CatDDM4; DDMD];
              c = 'b';
            end

            FilenameDDM2 = append(MainFolder, filesep, DateFolder{i}, filesep, CellFolder{j}, filesep, 'data', filesep, 'DDMOutput.mat');
            DDMSignal = load(FilenameDDM2);
            DDMSignal = DDMSignal.DDMOutput;

            [~,Idx] = min(abs(DDMSignal.QVector-1));
            Signal = DDMSignal.DDMSignalValue{Idx, 1};
            plot(DDMSignal.Time{Idx, 1}(1,1:20), DDMSignal.DDMSignalValue{Idx, 1}(1,1:20),'Color', c)
            xlabel('Time (s)')
            ylabel('ISF')
            title('MCF7 ISF - DDM')
            hold on
            CatIdx = CatIdx + 1;
            
        catch
            continue
        end
    end
end


DiffusionTrack = [[CatTrack1(:,1); nan(5,1)], [CatTrack2(:,1)], [CatTrack3(:,1); nan(4,1)]];
ViscosityTrack = [[CatTrack1(:,3); nan(5,1)], [CatTrack2(:,3)], [CatTrack3(:,3); nan(4,1)]];
DiffusionDDM = [[CatDDM1(:,1); nan(5,1)], [CatDDM2(:,1)], [CatDDM3(:,1); nan(4,1)]];

figure()
boxplot(DiffusionTrack,'Labels',{'Cat 1','Cat 2', 'Cat 3', 'Cat 4'})
ylabel('Diffusion (µm^2s^-^1)')
title('MCF7 Diffusion - SPT')

% figure()
% boxplot(ViscosityTrack,'Labels',{'Cat 1','Cat 2', 'Cat 3'})
% ylabel('Viscosity (cP)')
% title('KM12C Viscosity - SPT')

figure()
boxplot(DiffusionDDM,'Labels',{'Cat 1','Cat 2', 'Cat 3', 'Cat 4'})
ylabel('Diffusion (µm^2s^-^1)')
title('MCF7 Diffusion - DDM')