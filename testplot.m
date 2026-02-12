%% Diffusion reference values from "Res.xlsx"

% Time points (seconds and minutes)
DiffData.time.min = [3 6];
DiffData.time.sec = [180 360];

%% bAA conditions
DiffData.bAA(1).name = '1x bAA';
DiffData.bAA(1).ch1 = [2.22 2.03];
DiffData.bAA(1).ch2 = [1.35 1.21];

DiffData.bAA(2).name = '2x bAA';
DiffData.bAA(2).ch1 = [2.62 2.45];
DiffData.bAA(2).ch2 = [1.56 1.38];

DiffData.bAA(3).name = '3x bAA';
DiffData.bAA(3).ch1 = [2.36 1.34];
DiffData.bAA(3).ch2 = [2.39 1.35];

DiffData.bAA(4).name = '4x bAA';
DiffData.bAA(4).ch1 = [2.32 2.32];
DiffData.bAA(4).ch2 = [1.23 1.24];

%% AA conditions
DiffData.AA(1).name = '2x AA';
DiffData.AA(1).ch1 = [2.04 2.04];
DiffData.AA(1).ch2 = [1.24 1.23];

DiffData.AA(2).name = '3x AA';
DiffData.AA(2).ch1 = [1.94 1.91];
DiffData.AA(2).ch2 = [1.19 1.20];


%% Color gradients
n_bAA = numel(DiffData.bAA);
n_AA  = numel(DiffData.AA);

colors_bAA = [linspace(0.2,0.8,n_bAA)', zeros(n_bAA,1), linspace(0.8,0.2,n_bAA)']; % blue→purple
colors_AA  = [linspace(0.2,0.8,n_AA)', linspace(0.8,0.2,n_AA)', zeros(n_AA,1)];   % green→yellow


fig = figure; hold on
for i = 1:n_bAA
    plot(DiffData.time.min, DiffData.bAA(i).ch1, ...
        '-o', 'LineWidth', 2, 'Color', colors_bAA(i,:), ...
        'DisplayName', DiffData.bAA(i).name);
end
xlabel('Time (min)')
ylabel('Diffusion')
title('Channel 1 – bAA')
legend show
legend('Location', 'Best')
ylim([0.5 2.75])
grid on
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch1_bAA.svg'))
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch1_bAA.png'))

fig = figure; hold on
for i = 1:n_bAA
    plot(DiffData.time.min, DiffData.bAA(i).ch2, ...
        '-o', 'LineWidth', 2, 'Color', colors_bAA(i,:), ...
        'DisplayName', DiffData.bAA(i).name);
end
xlabel('Time (min)')
ylabel('Diffusion')
title('Channel 2 – bAA')
legend show
legend('Location', 'Best')
ylim([0.5 2.75])
grid on
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch2_bAA.svg'))
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch2_bAA.png'))

fig = figure; hold on
for i = 1:n_AA
    plot(DiffData.time.min, DiffData.AA(i).ch1, ...
        '-o', 'LineWidth', 2, 'Color', colors_AA(i,:), ...
        'DisplayName', DiffData.AA(i).name);
end
xlabel('Time (min)')
ylabel('Diffusion')
title('Channel 1 – AA')
legend show
legend('Location', 'Best')
ylim([0.5 2.75])
grid on
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch1_AA.svg'))
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch1_AA.png'))

fig = figure; hold on
for i = 1:n_AA
    plot(DiffData.time.min, DiffData.AA(i).ch2, ...
        '-o', 'LineWidth', 2, 'Color', colors_AA(i,:), ...
        'DisplayName', DiffData.AA(i).name);
end
xlabel('Time (min)')
ylabel('Diffusion')
title('Channel 2 – AA')
legend show
legend('Location', 'Best')
grid on
ylim([0.5 2.75])
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch2_AA.svg'))
saveas(fig, append('D:\Polymer Dynamics', filesep, 'Ch2_AA.png'))

