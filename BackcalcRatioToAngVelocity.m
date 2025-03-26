clc 
clear 
close all;
%calibration info
path2RotCal = 'D:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\100ms_exp';

%file info
MainFolder = 'D:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\100ms_exp';
subFolders = {'sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5', 'sample_6', 'sample_7', 'sample_8'};
ExpTime = 0.100; % in sec

%% Get RotCalibration info
RotCalib = open(append(path2RotCal, filesep, 'RotCalib.mat'));
RotCalib = RotCalib.RotCalib;
Amplitude = RotCalib.I_mean;
NumRotations = 2; %estimated number of rotations

%% Calculate angular velocity from the traces
AngSpeed = [];
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

            time_lags = tau_values*ExpTime;
            p = fit(time_lags(1, 70:90).', MSD(q, 70:90).', 'a*x+b'); %% take range where the slope is the steepest
            coeff = coeffvalues(p);
            D_theta = coeff(1) / 2;
            if isreal(sqrt(2*D_theta)*180/pi)
                angular_speed(q, 1) = sqrt(2*D_theta)*180/pi;
            else
                angular_speed(q, 1) = NaN;
            end

        end
    end
    Traces.AngSpeed = angular_speed;
    CommonTraces = Traces;
    Filename = append(MainFolder, filesep, subFolders{a}, filesep, "CommonTraces.mat");
    save(Filename, "CommonTraces")

    AngSpeed = [AngSpeed; angular_speed];
end

AngSpeed(isnan(AngSpeed)) = [];
AngularSpeed.All = AngSpeed;
AngularSpeed.mean = mean(AngSpeed);
AngularSpeed.median = median(AngSpeed);
AngularSpeed.std = std(AngSpeed);
save(append(path2RotCal, filesep, 'AngularSpeeds.mat'), 'AngularSpeed')