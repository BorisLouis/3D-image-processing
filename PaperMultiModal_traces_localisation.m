load('S:\Dual Color\20251015_localisation_error\Long_trace_from_water.mat')

%% Prepare coordinates
Trace1 = Trace1./1000;
coords = [Trace1.row, Trace1.col, Trace1.z];

%% 1️⃣ Add random displacement (mean = 25, std = 25)
displacementx = 69.708./1000 * randn(size(coords(:,1)));
displacementy = 178.757./1000 * randn(size(coords(:,2)));
displacementz = 119.406./1000 * randn(size(coords(:,3)));

Trace_displaced = [coords(:,1) + displacementx, ...
                   coords(:,2) + displacementy, ...
                   coords(:,3) + displacementz];

Trace_displaced_tbl = array2table(Trace_displaced, 'VariableNames', {'row','col','z'});

%% 2️⃣ Apply manual scale and translation
scale = 0.9713;
translation = [727.5575, 1049.7215, -1152.536]./1000;
Trace_transformed = Trace_displaced * scale + translation;
Trace_transformed_tbl = array2table(Trace_transformed, 'VariableNames', {'row','col','z'});

%% 3️⃣ Plot settings
% Define colors
green = [0.4660 0.6740 0.1880];   % #77ac30
red   = [0.6350 0.0740 0.0740];   % #a2142f

% Define target size in millimeters
width_mm = 180;
height_mm = 110;
mm_to_inch = 1/25.4;

Fig1 = figure('Units', 'inches', ...
              'Position', [2 2 width_mm*mm_to_inch height_mm*mm_to_inch], ...
              'Color', 'w');

% Plot traces
plot3(Trace1.row, Trace1.col, Trace1.z, 'Color', green, 'LineWidth', 0.75);
hold on
plot3(Trace_displaced_tbl.row, Trace_displaced_tbl.col, Trace_displaced_tbl.z, ...
      'Color', red, 'LineWidth', 0.75);

% Axis and grid formatting
grid on
ax = gca;
axis equal

xlim([21 24])
ylim([1 4])
zlim([-3 0])
xlabel('X (µm)', 'FontName', 'Arial', 'FontSize', 8)
ylabel('Y (µm)', 'FontName', 'Arial', 'FontSize', 8)
zlabel('Z (µm)', 'FontName', 'Arial', 'FontSize', 8)

ax.FontName = 'Arial';
ax.FontSize = 6;
ax.LineWidth = 0.75;
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';

% Make the grid thinner and lighter
ax.GridColor = [0.8 0.8 0.8];
ax.MinorGridColor = [0.8 0.8 0.8];
ax.GridAlpha = 1;
ax.MinorGridAlpha = 1;
ax.GridLineStyle = '-';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.Box = 'on';
ax.Layer = 'top';

% Adjust grid line thickness (indirectly by replotting axes)
set(findall(gca, 'Type', 'Line', '-and', 'Color', [0.15 0.15 0.15]), 'LineWidth', 0.25);

% Ensure the axis has correct aspect ratio for 3D view
daspect([1 1 1]);

%% 4️⃣ Export
set(Fig1, 'Renderer', 'painters'); % vector output
saveas(Fig1, ...
    'C:\Users\steve\Downloads\TraceColocalizationDisplaced.svg');

disp('✅ Figure exported as vector PDF, ready for Illustrator (50 x 35 mm)');





Fig2 = figure('Units', 'inches', ...
              'Position', [2 2 width_mm*mm_to_inch height_mm*mm_to_inch], ...
              'Color', 'w');

% Plot traces
plot3(Trace1.row, Trace1.col, Trace1.z, 'Color', green, 'LineWidth', 0.75);
hold on
plot3(Trace_transformed_tbl.row, Trace_transformed_tbl.col, Trace_transformed_tbl.z, ...
      'Color', red, 'LineWidth', 0.75);

% Axis and grid formatting
grid on
ax = gca;
axis equal
xlim([21 24])
ylim([1 4])
zlim([-3 0])
xlabel('X', 'FontName', 'Arial', 'FontSize', 8)
ylabel('Y', 'FontName', 'Arial', 'FontSize', 8)
zlabel('Z', 'FontName', 'Arial', 'FontSize', 8)

ax.FontName = 'Arial';
ax.FontSize = 6;
ax.LineWidth = 0.75;
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';

% Make the grid thinner and lighter
ax.GridColor = [0.8 0.8 0.8];
ax.MinorGridColor = [0.8 0.8 0.8];
ax.GridAlpha = 1;
ax.MinorGridAlpha = 1;
ax.GridLineStyle = '-';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.Box = 'on';
ax.Layer = 'top';

% Adjust grid line thickness (indirectly by replotting axes)
set(findall(gca, 'Type', 'Line', '-and', 'Color', [0.15 0.15 0.15]), 'LineWidth', 0.25);

% Ensure the axis has correct aspect ratio for 3D view
daspect([1 1 1]);

%% 4️⃣ Export
set(Fig2, 'Renderer', 'painters'); % vector output
saveas(Fig2, ...
    'C:\Users\steve\Downloads\TraceColocalizationTransformed.svg');

disp('✅ Figure exported as vector PDF, ready for Illustrator (50 x 35 mm)');

