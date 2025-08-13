function [FilePath, Experiment, Filename, Dimension, ExpTime, Temp, Radius, DiffFit, MinSize, Ext, ParticleType, path2RotCal] = CalcMSDinfoGUI()
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
    ParticleType = '';
    path2RotCal = '';

    % Create UI Figure
    fig = uifigure('Name', 'Parameter Input', ...
                   'Position', [100 100 450 650]);

    y = 600; dy = 40;

    % FilePath
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'File Path:');
    editFilePath = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', '');

    % Experiment
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Experiment:');
    dropdownExperiment = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'Tracking', 'Tracking-Segmentation', 'Tracking-Phase', 'Rotational Tracking'}, ...
        'Value', 'Tracking-Segmentation');

    % Filename
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Filename:');
    dropdownFilename = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'trackResults', 'trackResults1', 'trackResults2', 'TracesWMask', 'TraceswPhase', 'Traces3DCommon'}, ...
        'Value', 'TracesWMask');

    % Dimension
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Dimension:');
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

    % Radius (numeric default)
    y = y - dy;
    labelRadius = uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Radius (µm):');
    editRadiusNumeric = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.020);
    editRadiusArray   = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', '[184 92]', 'Visible', 'off');

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

    % ParticleType (hidden initially)
    y = y - dy;
    labelParticleType = uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Particle Type:', 'Visible', 'off');
    dropdownParticleType = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'Bipyramid', 'ellipsoid', 'rod', 'cilinder'}, ...
        'Value', 'Bipyramid', 'Visible', 'off');

    % path2RotCal (hidden initially)
    y = y - dy;
    labelPath2RotCal = uilabel(fig, 'Position', [20 y 100 22], 'Text', 'path2RotCal:', 'Visible', 'off');
    editPath2RotCal = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', '', 'Visible', 'off');

    % OK button
    btnOK = uibutton(fig, 'push', 'Text', 'OK', ...
        'Position', [170 y - 50 100 30], ...
        'ButtonPushedFcn', @(btn,event) submitCallback());

    % Show/hide extra fields when Experiment changes
    dropdownExperiment.ValueChangedFcn = @(dd,event) toggleRotationalFields();

    % Toggle function
    function toggleRotationalFields()
        if strcmp(dropdownExperiment.Value, 'Rotational Tracking')
            labelParticleType.Visible = 'on';
            dropdownParticleType.Visible = 'on';
            labelPath2RotCal.Visible = 'on';
            editPath2RotCal.Visible = 'on';
            editRadiusNumeric.Visible = 'off';
            editRadiusArray.Visible = 'on';
            labelRadius.Text = 'Radius array:';
        else
            labelParticleType.Visible = 'off';
            dropdownParticleType.Visible = 'off';
            labelPath2RotCal.Visible = 'off';
            editPath2RotCal.Visible = 'off';
            editRadiusNumeric.Visible = 'on';
            editRadiusArray.Visible = 'off';
            labelRadius.Text = 'Radius (µm):';
        end
    end

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
        DiffFit = editDiffFit.Value;
        MinSize = editMinSize.Value;
        Ext = dropdownExt.Value;

        if strcmp(Experiment, 'Rotational Tracking')
            ParticleType = dropdownParticleType.Value;
            path2RotCal = editPath2RotCal.Value;
            Radius = str2num(editRadiusArray.Value); %#ok<ST2NM>
        else
            Radius = editRadiusNumeric.Value;
        end

        uiresume(fig);
        delete(fig);
    end
end