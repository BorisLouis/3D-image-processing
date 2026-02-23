%% ============================================================
%  Scatter plot of viscosity vs nanoparticle size
%  3 sub-columns per size: SPT (left), TARDIS (middle), DDM (right)
%  Equal spacing of 100, 200, 500, 1000 nm
% =============================================================

clear; clc; close all;

%% Load data
filePath = 'E:\Data Steven - GEMs\data paper\NPs\Viscosity_Grouped.csv';
data = readmatrix(filePath);

%% Nanoparticle sizes (categorical spacing)
sizes = [100 200 500 1000];
baseX = 1:length(sizes);   % 1 2 3 4 → equal spacing

%% Extract LOW concentration (rows 4–9 in your sheet)
low = data(4:9, :);

low_DDM    = [low(:,1)  low(:,4)  low(:,7)  low(:,10)];
low_SPT    = [low(:,2)  low(:,5)  low(:,8)  low(:,11)];
low_TARDIS = [low(:,3)  low(:,6)  low(:,9)  low(:,12)];

%% Extract HIGH concentration (rows 16–21)
high = data(16:21, :);

high_DDM    = [high(:,1)  high(:,4)  high(:,7)  high(:,10)];
high_SPT    = [high(:,2)  high(:,5)  high(:,8)  high(:,11)];
high_TARDIS = [high(:,3)  high(:,6)  high(:,9)  high(:,12)];

%% Plot settings
offset = 0.25;        % spacing between SPT / TARDIS / DDM columns
jitterWidth = 0.06;   % horizontal jitter width

% Colors (consistent across plots)
color_DDM    = [0 0.4470 0.7410];      % blue
color_SPT    = [0.8500 0.3250 0.0980]; % orange
color_TARDIS = [0.4660 0.6740 0.1880]; % green


%% ============================================================
%% HIGH CONCENTRATION
%% ============================================================
figure; hold on;

for i = 1:length(baseX)

    % Define mini-column centers
    x_spt    = baseX(i) - offset;
    x_tardis = baseX(i);
    x_ddm    = baseX(i) + offset;

    % Scatter with jitter
    scatter(x_spt + (rand(size(high_SPT(:,i))) - 0.5)*jitterWidth, ...
            high_SPT(:,i), 60, color_SPT, 'filled');

    scatter(x_tardis + (rand(size(high_TARDIS(:,i))) - 0.5)*jitterWidth, ...
            high_TARDIS(:,i), 60, color_TARDIS, 'filled');

    scatter(x_ddm + (rand(size(high_DDM(:,i))) - 0.5)*jitterWidth, ...
            high_DDM(:,i), 60, color_DDM, 'filled');
end

xticks(baseX);
xticklabels(string(sizes));
xlabel('Nanoparticle size (nm)');
ylabel('Viscosity');
title('High Concentration');
ylim([0 2])
yline(0.953)
legend({'SPT','TARDIS','DDM'}, 'Location','best');
box on;
hold off;


%% ============================================================
%% LOW CONCENTRATION
%% ============================================================
figure; hold on;

for i = 1:length(baseX)

    x_spt    = baseX(i) - offset;
    x_tardis = baseX(i);
    x_ddm    = baseX(i) + offset;

    scatter(x_spt + (rand(size(low_SPT(:,i))) - 0.5)*jitterWidth, ...
            low_SPT(:,i), 60, color_SPT, 'filled');

    scatter(x_tardis + (rand(size(low_TARDIS(:,i))) - 0.5)*jitterWidth, ...
            low_TARDIS(:,i), 60, color_TARDIS, 'filled');

    scatter(x_ddm + (rand(size(low_DDM(:,i))) - 0.5)*jitterWidth, ...
            low_DDM(:,i), 60, color_DDM, 'filled');
end

xticks(baseX);
xticklabels(string(sizes));
xlabel('Nanoparticle size (nm)');
ylabel('Viscosity');
title('Low Concentration');
ylim([0 2])
yline(0.953)
legend({'SPT','TARDIS','DDM'}, 'Location','best');
box on;
hold off;
