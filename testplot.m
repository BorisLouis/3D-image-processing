%% --- CONSTANTS ---
T = 296.15;               % K
long_axis_nm = 182;
short_axis_nm = 91;
a = (long_axis_nm/2)*1e-9;   % m
b = (short_axis_nm/2)*1e-9;  % m
kB = 1.380649e-23;           % J/K

% Broersma coefficients
p = a/b;
f_par  = 4*pi*a ./ (log(2*p) - 0.5);
f_perp = 8*pi*a ./ (log(2*p) + 0.5);
f_avg = (f_par + 2*f_perp)/3;

%% --- INPUT: mean ± std (viscosity in cP) ---
data.trans = [ ...
    1.23 0.90;  % 3 min
    1.19 0.72;  % 5
    1.17 1.98;  % 7
    1.20 0.71;  % 9
    1.25 1.04;  % 11
    2.06 0.50;  % 13
    10.35 8.07; % 15
    16.96 11.10;% 17
    13.45 10.83 % 19
];

data.rot = [ ...
    200 80;    % 3 min (estimated from plot)
    250 60;    % 5 min (estimated from plot)
    180 70;    % 7 min (estimated from plot)
    10  8; 
    1.97 0.59;  % 11
    1.84 1.55;  % 13
    13.28 8.22; % 15
    16.72 15.36;% 17
    40.62 38.99 % 19
];

time = [3 5 7 9 11 13 15 17 19];

%% --- FUNCTION: viscosity → diffusion ---
convertViscToDiff = @(eta) kB*T ./ ((eta*1e-3) * f_avg) .* 10^(12);

%% --- RECONSTRUCT DISTRIBUTIONS AND CONVERT TO DIFFUSION ---
rng(1) % reproducible
N = 1000; % number of synthetic samples per box

D_trans = cell(numel(time),1);
D_rot   = cell(numel(time),1);

for i = 1:numel(time)
    % Translational
    mu = data.trans(i,1); sd = data.trans(i,2);
    if ~isnan(mu)
        viscSamples = max(normrnd(mu, sd, [N,1]), 1e-4);
        D_trans{i} = convertViscToDiff(viscSamples);
    end

    % Rotational
    mu = data.rot(i,1); sd = data.rot(i,2);
    if ~isnan(mu)
        viscSamples = max(normrnd(mu, sd, [N,1]), 1e-4);
        D_rot{i} = convertViscToDiff(viscSamples);
    end
end

%% --- BUILD DATA FOR BOXCHART ---
figure; hold on; box on;

% Translational (green)
for i = 1:numel(time)
    if ~isempty(D_trans{i})
        boxchart(ones(size(D_trans{i}))*time(i)-0.25, D_trans{i}, ...
                 'BoxFaceColor',[0.4 0.9 0.4], 'BoxEdgeColor','k');
    end
end

% Rotational (orange)
for i = 1:numel(time)
    if ~isempty(D_rot{i})
        boxchart(ones(size(D_rot{i}))*time(i)+0.25, D_rot{i}, ...
                 'BoxFaceColor',[1 0.8 0.3], 'BoxEdgeColor','k');
    end
end

%% --- STYLE ---
set(gca,'YScale','log');
xlabel('Polymerization time (min)');
ylabel('Diffusion coefficient (m^2/s)');
ylim([0 15])
title('Diffusion coefficients reconstructed from viscosity boxplots');
legend({'Translational','Rotational'},'Location','best');
grid on;
