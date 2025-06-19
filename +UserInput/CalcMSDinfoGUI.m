function [FilePath, Experiment, Filename, Dimension, ExpTime, Temp, Radius, DiffFit, MinSize, Ext] = CalcMSDinfoGUI()
     % Initialize default outputs
    FilePath = '';
    Experiment = '';
    Filename = '';
    ExpTime = NaN;
    Temp = NaN;
    Radius = NaN;
    DiffFit = NaN;
    MinSize = NaN;
    Ext = '';

    % Create UI Figure
    fig = uifigure('Name', 'Parameter Input', ...
                   'Position', [100 100 450 520]);

    y = 470; dy = 40;

    % FilePath
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'File Path:');
    editFilePath = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', '');

    % Experiment
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Experiment:');
    dropdownExperiment = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'Tracking', 'Tracking-Segmentation', 'Tracking-Phase'}, ...
        'Value', 'Tracking-Segmentation');

    % Filename
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Filename:');
    dropdownFilename = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'trackResults', 'trackResults1', 'trackResults2', 'TracesWMask'}, ...
        'Value', 'TracesWMask');

    % Dimension
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Dimension:'); % This is the label
    dropdownDimension = uidropdown(fig, ...
        'Position', [130 y 280 22], ...
        'Items', {'2D', '3D'}, ...
        'Value', '2D');

    % ExpTime
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'ExpTime (s):');
    editExpTime = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.030);

    % Temp
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Temp (K):');
    editTemp = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 296.15);

    % Radius
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Radius (µm):');
    editRadius = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.020);

    % DiffFit
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'DiffFit:');
    editDiffFit = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 4);

    % MinSize
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'MinSize:');
    editMinSize = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 20);

    % Extension
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Extension:');
    dropdownExt = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'.mat', '.csc', '.xlsx'}, 'Value', '.mat');

    % OK button
    btnOK = uibutton(fig, 'push', 'Text', 'OK', ...
        'Position', [170 y - 50 100 30], ...
        'ButtonPushedFcn', @(btn,event) submitCallback());

    % Block execution until OK is pressed
    uiwait(fig);

    % --- Submit: Collect values and close UI ---
    function submitCallback()
        FilePath = editFilePath.Value;
        Experiment = dropdownExperiment.Value;
        Filename = dropdownFilename.Value;
        Dimension = dropdownDimension.Value;
        ExpTime = editExpTime.Value;
        Temp = editTemp.Value;
        Radius = editRadius.Value;
        DiffFit = editDiffFit.Value;
        MinSize = editMinSize.Value;
        Ext = dropdownExt.Value;

        uiresume(fig);
        delete(fig);
    end
end
