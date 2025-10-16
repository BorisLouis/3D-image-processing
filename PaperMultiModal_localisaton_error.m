Dimensionality = 3;
Exptime = 0.250;
R = 100*10^(-9);

MainFolder = 'S:\Dual Color\Localisation_error\20251016\ExpTIme_250ms';
Folders = dir(MainFolder);

%%% for channel 1
diffx = [];
diffy = [];
diffz = [];
diffr = [];
for j = 3:size(Folders, 1)
    if Folders(j).isdir == 1
        Filename = append(Folders(j).folder, filesep, Folders(j).name, filesep, 'traces3D_noSRCal1.mat');
        load(Filename);
        for i = 1:size(traces, 1)
            if size(traces{i, 1}, 1) > 10
                CurrTrace = traces{i, 1};
                diffx = [diffx; median(abs(diff(CurrTrace.row)))];
                diffy = [diffy; median(abs(diff(CurrTrace.col)))];
                diffz = [diffz; median(abs(diff(CurrTrace.z)))];
                diffr = [diffr; sqrt(mean(abs(diff(CurrTrace.row))).^2 + mean(abs(diff(CurrTrace.col))).^2 + mean(abs(diff(CurrTrace.z))).^2)];
            end
        end
    end
end

Results.Error.errorX = median(diffx)*10^(-3);
Results.Error.errorY = median(diffy)*10^(-3);
Results.Error.errorZ = median(diffz)*10^(-3);
Results.Error.errorR = median(diffr)*10^(-3);

Results.Diffusion.Dx = ((Results.Error.errorX).^2)./(2*1*Exptime);
Results.Diffusion.Dy = ((Results.Error.errorY).^2)./(2*1*Exptime);
Results.Diffusion.Dz = ((Results.Error.errorZ).^2)./(2*1*Exptime);
Results.Diffusion.Dr = ((Results.Error.errorR).^2)./(2*Dimensionality*Exptime);

Results.Viscosity.nx = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dx*10^(-12))*1000;
Results.Viscosity.ny = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dy*10^(-12))*1000;
Results.Viscosity.nz = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dz*10^(-12))*1000;
Results.Viscosity.nr = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dr*10^(-12))*1000;

SaveFile = append(MainFolder, filesep, 'ResultsLocErrorCh1.mat');
save(SaveFile, "Results")


%%% for channel 2
diffx = [];
diffy = [];
diffz = [];
diffr = [];
for j = 3:size(Folders, 1)
    if Folders(j).isdir == 1
        Filename = append(Folders(j).folder, filesep, Folders(j).name, filesep, 'traces3D_noSRCal2.mat');
        load(Filename);
        for i = 1:size(traces, 1)
            if size(traces{i, 1}, 1) > 10
                CurrTrace = traces{i, 1};
                diffx = [diffx; median(abs(diff(CurrTrace.row)))];
                diffy = [diffy; median(abs(diff(CurrTrace.col)))];
                diffz = [diffz; median(abs(diff(CurrTrace.z)))];
                diffr = [diffr; sqrt(mean(abs(diff(CurrTrace.row))).^2 + mean(abs(diff(CurrTrace.col))).^2 + mean(abs(diff(CurrTrace.z))).^2)];
            end
        end
    end
end

Results.Error.errorX = median(diffx)*10^(-3);
Results.Error.errorY = median(diffy)*10^(-3);
Results.Error.errorZ = median(diffz)*10^(-3);
Results.Error.errorR = median(diffr)*10^(-3);

Results.Diffusion.Dx = ((Results.Error.errorX).^2)./(2*1*Exptime);
Results.Diffusion.Dy = ((Results.Error.errorY).^2)./(2*1*Exptime);
Results.Diffusion.Dz = ((Results.Error.errorZ).^2)./(2*1*Exptime);
Results.Diffusion.Dr = ((Results.Error.errorR).^2)./(2*Dimensionality*Exptime);

Results.Viscosity.nx = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dx*10^(-12))*1000;
Results.Viscosity.ny = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dy*10^(-12))*1000;
Results.Viscosity.nz = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dz*10^(-12))*1000;
Results.Viscosity.nr = (1.380649*10^-(23)*296.15)./(6*pi*R*Results.Diffusion.Dr*10^(-12))*1000;

SaveFile = append(MainFolder, filesep, 'ResultsLocErrorCh2.mat');
save(SaveFile, "Results")