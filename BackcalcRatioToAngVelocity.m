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
        Diff = [];
        time = [];
        tau_values = [];
        deltas = [];
        if isnan(Traces.I(q))
            angular_speed(q, 1) = NaN;
        else
            Diff = Traces.Diff{q};
            time = Traces.Diff{q};

            Theta = 0.5*real(acos(Diff./Amplitude));
            tau = Traces.Time(q,:);

            [msadTheta] = MSD.Rotational.calc(Theta, tau, ExpTime);
            DTheta   = MSD.Rotational.getDiffCoeff(msadTheta,tau,4,'2D');

            if isreal(sqrt(DTheta)*180/pi)
                angular_speed(q, 1) = sqrt(DTheta)*180/pi;
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