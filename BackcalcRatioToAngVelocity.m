clc 
clear 
close all;
%calibration info
path2RotCal = 'E:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';

%file info
MainFolder = 'E:\Rotational Tracking\20250228_AuBPs_184x92_calib\2DCal_184x91_rotational\10ms_exp';
subFolders = {'sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5', 'sample_6', 'sample_7', 'sample_8'};
ExpTime = 0.010; % in sec

%% Get RotCalibration info
RotCalib = open(append(path2RotCal, filesep, 'RotCalib.mat'));
RotCalib = RotCalib.RotCalib;
Amplitude = RotCalib.I_mean;
NumRotations = 2; %estimated number of rotations

%% Calculate angular velocity from the traces
AngSp1 = [];
AngSp2 = [];
AngSp3 = [];
AngSp4 = [];
MSADPeaks = [];
f = waitbar(0, 'initializing');

for a = 1:size(subFolders, 2)
    Traces = open(append(MainFolder, filesep, subFolders{a}, filesep, "CommonTraces.mat"));
    Traces = Traces.CommonTraces;

    angular_speed = [];
    AngSpeed1 = [];
    AngSpeed2 = [];
    AngSpeed3 = [];
    AngSpeed4 = [];
    for q = 1:size(Traces, 1)
        waitbar(q./size(Traces, 1), f, append('Calculating: sample ', num2str(a), ' out of ', num2str(size(subFolders, 2))))
        Diff = [];
        time = [];
        tau_values = [];
        deltas = [];
        if isnan(Traces.I(q))
            angular_speed(q, 1) = NaN;
        else
            Diff = Traces.Diff{q};
            r = fit(Traces.Time(q, :)', Traces.Diff{q}, ['0.3*cos(0.*x)+c']);
            figure()
            plot(r, Traces.Time(q, :)', Traces.Diff{q})
            coeff = coeffvalues(r);
            Diff = Diff - coeff(3);
            time = Traces.Diff{q};

            %%% Check1
            Diffsmooth = medfilt1(Diff, 90);
            [pks1,locs1,w,p] = findpeaks(Diffsmooth);
            [pks,locs2,w,p] = findpeaks(-Diffsmooth);
            Int1 = [diff(locs1); diff(locs2)];
            AngSpeed1 = [AngSpeed1; 360./(Int1 * ExpTime)/4];

            %%% Check2
            Theta = 0.25*real(acos(Diff./coeff(1)));
            Thetasmooth = medfilt1(Theta, 90);
            [pks,locs1,w,p] = findpeaks(Thetasmooth);
            [pks,locs2,w,p] = findpeaks(-Thetasmooth);
            peaks = sort([locs1; locs2]);
            for w = 1:size(peaks, 1)-1
                try
                    segment = Theta(peaks(w)+50:peaks(w+1)-10, 1);
                    timelag = [1:size(segment, 1)]'*ExpTime;
                    g = fit(timelag, segment, 'a*x +b');
                    coeff = coeffvalues(g);
                    % figure()
                    % plot(f, timelag, segment)
                    AngSpeed2 = [AngSpeed2; abs(coeff(1))*180/pi];
                catch
                    AngSpeed2 = [AngSpeed2; NaN];
                end
            end

            tau = Traces.Time(q,:);

            [msadTheta] = MSD.Rotational.calc(Theta, tau, ExpTime);
            % figure()
            % plot(msadTheta(2,:), msadTheta(1,:))
            % xlabel('Time (s)')
            % ylabel('MSAD (rad^2)')

            %%% Check 3
            [pksTop,locs1,w,p] = findpeaks(medfilt1(msadTheta(1,:), 50));
            [pks,locs2,w,p] = findpeaks(medfilt1(-msadTheta(1,:), 50));
            Interval = [diff(locs1), diff(locs2)]';
            AngSpeed3 = [AngSpeed3; 360./(Interval*ExpTime)/4];

            %%% save peaks
            MSADPeaks{q, 1} = Diff;
            MSADPeaks{q, 2} = mean(pksTop);

            %%% Try to get diff coeff
            StartPoint = [mean(pksTop)/2, mean(Interval)/180*pi*ExpTime];
            [bier, gov] = fit(msadTheta(2, 1:500)', msadTheta(1, 1:500)', 'a*(1-cos(b*x))', 'StartPoint', StartPoint);
            % figure()
            % plot(bier, msadTheta(2, :), msadTheta(1, :))
            coeff = coeffvalues(bier);
            AngSpeed4 = [AngSpeed4; coeff(2)*180/pi/4];


        end
    end

    RotSpeedCheck.AngSpeed1 = AngSpeed1;
    RotSpeedCheck.AngSpeed2 = AngSpeed2;
    RotSpeedCheck.AngSpeed3 = AngSpeed3;
    Filename = append(MainFolder, filesep, subFolders{a}, filesep, "RotSpeedCheck.mat");
    save(Filename, "RotSpeedCheck")

    AngSp1 = [AngSp1; AngSpeed1];
    AngSp2 = [AngSp2; AngSpeed2];
    AngSp3 = [AngSp3; AngSpeed3];
    AngSp4 = [AngSp4; AngSpeed4];
end

CheckSpeed1.All = AngSp1;
CheckSpeed1.mean = nanmedian(AngSpeed1);
CheckSpeed1.std = std(AngSpeed1);
CheckSpeed2.All = AngSp2;
Checkpeed2.mean = nanmedian(AngSpeed2);
CheckSpeed2.std = std(AngSpeed2);
CheckSpeed3.All = AngSp3;
CheckSpeed3.mean = nanmedian(AngSpeed3);
CheckSpeed3.std = std(AngSpeed3);
CheckSpeed4.All = AngSp4;
CheckSpeed4.mean = nanmedian(AngSpeed4);
CheckSpeed4.std = std(AngSpeed4);
AngularSpeed.Check1 = CheckSpeed1;
AngularSpeed.Check2 = CheckSpeed2;
AngularSpeed.Check3 = CheckSpeed3;
AngularSpeed.Check4 = CheckSpeed4;
save(append(path2RotCal, filesep, 'AngularSpeeds.mat'), 'AngularSpeed')