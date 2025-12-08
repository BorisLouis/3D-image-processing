function [FilePath, Experiment, Filename, Dimension, ExpTime, Temp, Radius, Radius2, DiffFit, MinSize, Ext, ParticleType, path2RotCal, CutTraces, ExpModal] = CalcMSDinfoGUI()
    % Initialize default outputs
    FilePath = '';
    Experiment = '';
    Filename = '';
    ExpTime = NaN;
    Temp = NaN;
    Radius = NaN;
    Radius2 = NaN;
    DiffFit = NaN;
    MinSize = NaN;
    Ext = '';
    ParticleType = '';
    path2RotCal = '';
    CutTraces = NaN;

    % Create UI Figure
    fig = uifigure('Name', 'Parameter Input', ...
                   'Position', [100 100 450 700]);

    y = 650; dy = 40;

    % FilePath
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'File Path:');
    editFilePath = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', 'S:\Dual Color\20250121_dualcolor\PS_200g_100r\20250122_PS_200_green_PS_100_red\sample2');

    % Experiment
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Experiment:');
    dropdownExperiment = uidropdown(fig, 'Position', [130 y 280 22], ...
        'Items', {'Tracking', 'Tracking-Segmentation', 'Tracking-Phase', 'Rotational Tracking', 'Dual color tracking'}, ...
        'Value', 'Dual color tracking');

    % Filename
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Filename:');
    dropdownFilename = uidropdown(fig, 'Position', [130 y 280 22],...
        'Items', {'trackResults', 'trackResults1', 'trackResults2', 'TracesWMask', 'TraceswPhase', 'Traces3DCommon', 'traces3D_', 'traces3D_noSRCal', 'Traces3D'}, ...
        'Value', 'trackResults');

    % Dimension
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Dimension:');
    dropdownDimension = uidropdown(fig, ...
        'Position', [130 y 280 22], ...
        'Items', {'2D', '3D'}, ...
        'Value', '3D');

    % ExpTime
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'ExpTime (s):');
    editExpTime = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.010);

    % Temp
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Temp (K):');
    editTemp = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 296.15);

    % Radius (numeric default)
    y = y - dy;
    labelRadius = uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Radius (µm):');
    editRadiusNumeric = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.150);
    editRadiusArray   = uieditfield(fig, 'text', 'Position', [130 y 280 22], 'Value', '[184 92]', 'Visible', 'off');

    % Radius2 (for Dual color tracking, hidden initially)
    y = y - dy;
    labelRadius2 = uilabel(fig, 'Position', [20 y 100 22], 'Text', 'Radius 2 (µm):', 'Visible', 'off');
    editRadius2 = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 0.10, 'Visible', 'off');

    % DiffFit
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'DiffFit:');
    editDiffFit = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 4);

    % MinSize
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'MinSize:');
    editMinSize = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 20);

    % CutTraces
    y = y - dy;
    uilabel(fig, 'Position', [20 y 100 22], 'Text', 'CutTraces:');
    editCutTraces = uieditfield(fig, 'numeric', 'Position', [130 y 280 22], 'Value', 50);

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

    % ExpModal (hidden initially)
    y = y - dy;
    labelExpModal = uilabel(fig, 'Position', [20 y 140 22], 'Text', 'ExpModal:', 'Visible', 'off');
    dropdownExpModal = uidropdown(fig, 'Position', [160 y 250 22], ...
        'Items', {'Test', 'Single Exponential', 'Bi Exponential', 'Stretched Exponential'}, ...
        'Value', 'Single Exponential', 'Visible', 'off');

    % OK button
    btnOK = uibutton(fig, 'push', 'Text', 'OK', ...
        'Position', [170 y - 50 100 30], ...
        'ButtonPushedFcn', @(btn,event) submitCallback());

    % Show/hide extra fields when Experiment changes
    dropdownExperiment.ValueChangedFcn = @(dd,event) toggleExperimentFields();

    % Toggle function
    function toggleExperimentFields()
        switch dropdownExperiment.Value
            case 'Rotational Tracking'
                labelParticleType.Visible = 'on';
                dropdownParticleType.Visible = 'on';
                labelPath2RotCal.Visible = 'on';
                editPath2RotCal.Visible = 'on';
                labelExpModal.Visible = 'on';
                dropdownExpModal.Visible = 'on';
                editRadiusNumeric.Visible = 'off';
                editRadiusArray.Visible = 'on';
                labelRadius.Text = 'Radius array:';
                labelRadius2.Visible = 'off';
                editRadius2.Visible = 'off';

            case 'Dual color tracking'
                labelParticleType.Visible = 'off';
                dropdownParticleType.Visible = 'off';
                labelPath2RotCal.Visible = 'off';
                editPath2RotCal.Visible = 'off';
                labelExpModal.Visible = 'off';
                editRadiusNumeric.Visible = 'on';
                editRadiusArray.Visible = 'off';
                labelRadius.Text = 'Radius (µm):';
                labelRadius2.Visible = 'on';
                editRadius2.Visible = 'on';

            otherwise
                labelParticleType.Visible = 'off';
                dropdownParticleType.Visible = 'off';
                labelPath2RotCal.Visible = 'off';
                editPath2RotCal.Visible = 'off';
                labelExpModal.Visible = 'off';
                dropdownExpModal.Visible = 'off';
                editRadiusNumeric.Visible = 'on';
                editRadiusArray.Visible = 'off';
                labelRadius.Text = 'Radius (µm):';
                labelRadius2.Visible = 'off';
                editRadius2.Visible = 'off';
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
        CutTraces = editCutTraces.Value;
        if CutTraces == 0
            CutTraces = NaN;
        end

        if strcmp(Experiment, 'Rotational Tracking')
            ParticleType = dropdownParticleType.Value;
            path2RotCal = editPath2RotCal.Value;
            Radius = str2num(editRadiusArray.Value); %#ok<ST2NM>
            Radius2 = NaN;
            ExpModal = dropdownExpModal.Value;
        elseif strcmp(Experiment, 'Dual color tracking')
            Radius = editRadiusNumeric.Value;
            Radius2 = editRadius2.Value;
            ExpModal = '';  % no value
        else
            Radius = editRadiusNumeric.Value;
            Radius2 = NaN;
            ExpModal = '';  % no value
        end

        uiresume(fig);
        delete(fig);
    end
end