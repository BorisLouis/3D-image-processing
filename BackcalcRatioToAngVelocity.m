clc 
clear 
close all;
%calibration info
path2RotCal = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';

%file info
MainFolder = 'S:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';
subFolders = {'sample_1', 'sample_2'};

%% Get RotCalibration info
RotCalib = open(append(path2RotCal, filesep, 'RotCalib.mat'));
RotCalib = RotCalib.RotCalib;
Amplitude = RotCalib.I_mean;
NumRotations = 2; %estimated number of rotations

%% Calculate angular velocity from the traces
for a = 1:size(subFolders, 2)
    Traces = open(append(MainFolder, filesep, subFolders{a}, filesep, "CommonTraces.mat"));
    Traces = Traces.CommonTraces;

    angular_speed = [];
    for q = 1:size(Traces, 1)
        MSD = [];
        Diff = [];
        time = [];
        tau_values = [];
        deltas = [];
        if isnan(Traces.I(q))
            angular_speed(q, 1) = NaN;
        else
            Diff = Traces.Diff{q};
            time = Traces.Diff{q};
            
            for i = 1:size(Diff, 1)
                Angles(i,1) = 0.25*real(acos(Diff(i)./Amplitude));
            end
            tau_values = 1:floor(length(Angles));
            
            for n = 1:size(tau_values,2) - 1
                dt = tau_values(n);
                for i = 1:size(Angles,1)-dt
                    deltas(i) = Angles(i+dt) - Angles(i);  % Difference between successive angles
                    deltas(i) = (deltas(i)).^2;
                end
                deltas(deltas == 0) = [];
                value = median(deltas);
                MSD(q,dt) = value;
            end
                 
            time_lags = tau_values*0.010;
            p = fit(time_lags(1, 50:80).', MSD(q, 50:80).', 'a*x+b');
            coeff = coeffvalues(p);
            D_theta = coeff(1) / 2;
            angular_speed(q, 1) = sqrt(2*D_theta)*180/pi;
        end
    end
    Traces.AngSpeed = angular_speed;
    CommonTraces = Traces;
    Filename = append(MainFolder, filesep, subFolders{a}, filesep, "CommonTraces.mat");
    save(Filename, "CommonTraces")
end

