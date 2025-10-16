function [info, info1, info2, file] = infoGUI(file)
%% infoGUI - Experiment Setup GUI
% Fully fixed version: all DDM fields store user input values

% Default outputs
info = struct();
info1 = struct();
info2 = struct();

% Create the UI figure
fig = uifigure('Name', 'Experiment Setup', 'Position', [100, 100, 1000, 650]);
gl = uigridlayout(fig, [2, 3]);
gl.RowHeight = {'9x', '1x'};
gl.ColumnWidth = {'1x', '1x', '1x'};

state = struct();

%% Column 1 – General Parameters
col1 = uigridlayout(gl, [14, 2]);

state.Ext = addDropdown(col1, 'Ext:', {'.ome.tif', '.tif', '.his', '.mpg', '.spe', '.lif'}, '.ome.tif');
state.Type = addDropdown(col1, 'Type:', {'normal', 'transmission'}, 'normal');
state.RunMethod = addDropdown(col1, 'Run Method:', {'run', 'load'}, 'load');
state.calibrate = addDropdown(col1, 'calibrate', {'true', 'false'}, 'false');

state.drawROI = addDropdown(col1, 'draw ROI:', {'off', 'channel1', 'channel2'}, 'off');
state.Dimension = addDropdown(col1, 'Dimension:', {'2D', '3D'}, '3D');
state.multiModal = addDropdown(col1, 'multiModal:', {'on', 'off'}, 'on');

channelTypes = {'Translational Tracking', 'Rotational Tracking', 'Segmentation', 'Phase', 'DDM'};
ch1Dropdown = addDropdown(col1, 'Channel 1:', channelTypes, 'Rotational Tracking');
ch2Dropdown = addDropdown(col1, 'Channel 2:', channelTypes, 'Rotational Tracking');

state.PxSize = addLabelField(col1, 'PxSize (nm):', '95');
state.FWHM = addLabelField(col1, 'FWHM (px):', '3');
state.Frame2Load = addLabelField(col1, 'Frame2Load:', 'all');
state.TestFrame = addLabelField(col1, 'Test Frame:', '10');
state.Rotational = addDropdown(col1, 'Rotational:', {'on', 'off'}, 'on');
state.Bipyramid = addLabelField(col1, 'Bipyramid (nm):', '[184 92]');
state.RotationalCalib = addDropdown(col1, 'Rotational Calib:', {'on', 'off'}, 'off');
state.RadTime = addLabelField(col1, 'Rad Time (°/s):', '25');

%% Channel Panels
channel1Panel = uipanel(gl, 'Title', 'Channel 1 Setup');
channel1Layout = uigridlayout(channel1Panel, [14, 2]);
channel1Layout.RowHeight = repmat({30}, 1, 14);

channel2Panel = uipanel(gl, 'Title', 'Channel 2 Setup');
channel2Layout = uigridlayout(channel2Panel, [14, 2]);
channel2Layout.RowHeight = repmat({30}, 1, 14);

state.channel1Controls = [];
state.channel2Controls = [];

%% OK Button
btn = uibutton(gl, 'Text', 'OK', 'ButtonPushedFcn', @(btn, event) onOK());
btn.Layout.Row = 2;
btn.Layout.Column = [1 3];

%% Setup logic
ch1Dropdown.ValueChangedFcn = @(src, event) updateChannel(1, src.Value);
ch2Dropdown.ValueChangedFcn = @(src, event) updateChannel(2, src.Value);
state.RotationalCalib.ValueChangedFcn = @(src, event) toggleRadTime();
state.multiModal.ValueChangedFcn = @(src, event) toggleMultiModal();
state.calibrate.ValueChangedFcn = @(src, event) toggleCalibrate();

% Initial setup
updateChannel(1, ch1Dropdown.Value);
updateChannel(2, ch2Dropdown.Value);
toggleRadTime();
toggleBipyramid();
toggleRotationalControls();
toggleMultiModal();
toggleCalibrate();

uiwait(fig);

%% --- Nested functions ---

function updateChannel(channelNum, mode)
    if channelNum == 1
        delete(channel1Layout.Children);
        state.channel1Controls = addChannelControls(channel1Layout, mode, 'channel1Controls');
    else
        delete(channel2Layout.Children);
        state.channel2Controls = addChannelControls(channel2Layout, mode, 'channel2Controls');
    end
    toggleBipyramid();
    toggleRotationalControls();
end

function toggleMultiModal()
    if strcmp(state.multiModal.Value, 'off')
        ch2Dropdown.Enable = 'off';
        channel2Panel.Visible = 'off';
        state.channel2Controls = struct();
        state.drawROI.Items = {'off', 'channel1'};
        if strcmp(state.drawROI.Value, 'channel2')
            state.drawROI.Value = 'off';
        end
    else
        ch2Dropdown.Enable = 'on';
        channel2Panel.Visible = 'on';
        state.drawROI.Items = {'off', 'channel1', 'channel2'};
    end
end

function toggleBipyramid()
    ch1 = ch1Dropdown.Value;
    ch2 = ch2Dropdown.Value;
    if strcmp(ch1, 'Rotational Tracking') && strcmp(ch2, 'Rotational Tracking')
        state.Bipyramid.Enable = 'on';
    else
        state.Bipyramid.Enable = 'off';
    end
end

function toggleRadTime()
    if strcmp(state.RotationalCalib.Value, 'on')
        state.RadTime.Enable = 'on';
    else
        state.RadTime.Enable = 'off';
    end
end

function toggleRotationalControls()
    ch1 = ch1Dropdown.Value;
    ch2 = ch2Dropdown.Value;

    if strcmp(ch1, 'Rotational Tracking') && strcmp(ch2, 'Rotational Tracking')
        state.Rotational.Value = 'on';
        state.Rotational.Enable = 'off';
        state.RotationalCalib.Enable = 'on';
    else
        state.Rotational.Value = 'off';
        state.Rotational.Enable = 'off';
        state.RotationalCalib.Value = 'off';
        state.RotationalCalib.Enable = 'off';
    end

    toggleRadTime();
end

function toggleCalibrate()
    if strcmp(state.calibrate.Value, 'false')
        state.drawROI.Value = 'off';
        state.drawROI.Enable = 'off';
    else
        state.drawROI.Enable = 'on';
    end
end

function onOK()
    % Collect info structure
    info.Dimension = state.Dimension.Value;
    info.PxSize = str2double(state.PxSize.Value);
    info.type = state.Type.Value;
    info.runMethod = state.RunMethod.Value;
    info.multiModal = strcmp(state.multiModal.Value, 'on');
    info.drawROI = state.drawROI.Value;
    info.calibrate = strcmp(state.calibrate.Value, 'true');

    if strcmp(state.Frame2Load.Value, 'all')
        info.frame2Load = state.Frame2Load.Value;
    else
        info.frame2Load = evalin('base', state.Frame2Load.Value);
    end

    info.Channel1 = ch1Dropdown.Value;
    info.Channel2 = ch2Dropdown.Value;
    info.FWHM = str2double(state.FWHM.Value);
    info.rotational = strcmp(state.Rotational.Value, 'on');
    info.rotationalCalib = strcmp(state.RotationalCalib.Value, 'on');
    info.TestFrame = str2double(state.TestFrame.Value);
    info.RadTime = str2double(state.RadTime.Value);
    info.Bipyramid = str2num(state.Bipyramid.Value);
    file.ext = state.Ext.Value;

    % Read controls for both channels; ensure they are structs (possibly empty)
    if ~isempty(fieldnames(state.channel1Controls))
        info1 = readControls(state.channel1Controls);
    else
        info1 = struct();
    end

    if strcmp(state.multiModal.Value, 'off') || isempty(fieldnames(state.channel2Controls))
        info2 = struct();    % no second-channel info => empty struct
    else
        info2 = readControls(state.channel2Controls);
    end

    % Restructure only when the control-struct contains fields
    if ~isempty(fieldnames(info1))
        info1 = restructureChannelInfo(info1, info.Channel1);
    else
        info1 = struct();
    end

    if ~isempty(fieldnames(info2))
        info2 = restructureChannelInfo(info2, info.Channel2);
    else
        info2 = struct();
    end

    uiresume(fig);
    delete(fig);
end

function out = restructureChannelInfo(inStruct, channelType)
    out = inStruct;
    switch true
        case contains(channelType, 'Tracking')
            out.detectParam.delta = inStruct.delta;
            out.detectParam.chi2 = inStruct.chi2;
            out.detectParam.consThresh = inStruct.consThresh;
            out.trackParam.radius = inStruct.track_radius;
            out.trackParam.memory = inStruct.track_memory;
            out = rmfield(out, {'delta', 'chi2', 'consThresh', 'track_radius', 'track_memory'});

        case contains(channelType, 'Phase')
            out.optics.dz = inStruct.dz;
            out.optics.NA = inStruct.NA;
            out.optics.NA_ill = inStruct.NA_ill;
            out.optics.n = inStruct.n;
            out.optics.lambda = inStruct.lambda;
            out.optics.dlambda = inStruct.dlambda;
            out.optics.alpha = inStruct.alpha;
            out.optics.kzT = inStruct.kzT;
            % mirrorX, mirrorZ and applyFourierMask are now plain values (strings 'true'/'false' or logical)
            if isfield(inStruct, 'mirrorX')
                out.proc.mirrorX = inStruct.mirrorX;
            end
            if isfield(inStruct, 'mirrorZ')
                out.proc.mirrorZ = inStruct.mirrorZ;
            end
            if isfield(inStruct, 'applyFourierMask')
                out.proc.applyFourierMask = inStruct.applyFourierMask;
            end
            out = rmfield(out, intersect(fieldnames(out), ...
                {'dz','NA','NA_ill','n','lambda','dlambda','alpha','kzT','mirrorX','mirrorZ','applyFourierMask'}));

        case contains(channelType, 'DDM')
            % Extract values (readControls now returns plain values)
            if isfield(inStruct, 'ParticleSize'),   out.ddmParam.ParticleSize = inStruct.ParticleSize;   end
            if isfield(inStruct, 'ExpTime'),       out.ddmParam.ExpTime = inStruct.ExpTime;             end
            if isfield(inStruct, 'Wavelength'),    out.ddmParam.Wavelength = inStruct.Wavelength;       end
            if isfield(inStruct, 'NA'),            out.ddmParam.NA = inStruct.NA;                       end
            if isfield(inStruct, 'Temp'),          out.ddmParam.Temp = inStruct.Temp;                   end
            if isfield(inStruct, 'Qmin'),          out.ddmParam.Qmin = inStruct.Qmin;                   end
            if isfield(inStruct, 'Qmax'),          out.ddmParam.Qmax = inStruct.Qmax;                   end
            if isfield(inStruct, 'Scanning'),      out.ddmParam.Scanning = inStruct.Scanning;           end
            if isfield(inStruct, 'CorrectBleaching'), out.ddmParam.CorrectBleaching = inStruct.CorrectBleaching; end
            if isfield(inStruct, 'AngularAnisotropy'), out.ddmParam.AngularAnisotropy = inStruct.AngularAnisotropy; end

            if isfield(inStruct, 'ROISize')
                out.ddmParam.ROISize = inStruct.ROISize;
            end
            if isfield(inStruct, 'FitRDiff')
                out.ddmParam.FitRDiff = inStruct.FitRDiff;
            end

            out = rmfield(out, intersect(fieldnames(out), ...
                {'ParticleSize','ExpTime','Wavelength','NA','Temp','Qmin','Qmax','Scanning','CorrectBleaching','AngularAnisotropy','ROISize','FitRDiff'}));
    end
end

%% --- Channel Control Creator ---
function controls = addChannelControls(layout, type, controlsName)
    controls = struct();
    switch type
        case {'Translational Tracking', 'Rotational Tracking'}
            controls.fitMethod = addDropdown(layout, 'fitMethod', {'Phasor', 'Gauss'}, 'Phasor');
            controls.zMethod = addDropdown(layout, 'zMethod', {'Intensity', '3DFit', 'PSFE'}, 'Intensity');
            controls.detectionMethod = addDropdown(layout, 'detectionMethod', {'Intensity', 'MaxLR'}, 'MaxLR');
            controls.IntCorr = addDropdown(layout, 'Intensity correction', {'on', 'off'}, 'on');
            controls.euDist = addLabelField(layout, 'euDist (nm)', '1500');
            controls.delta = addLabelField(layout, 'delta', '10');
            controls.chi2 = addLabelField(layout, 'chi2', '35');
            controls.consThresh = addLabelField(layout, 'consThresh', '4');
            controls.track_radius = addLabelField(layout, 'track.radius (nm)', '1500');
            controls.track_memory = addLabelField(layout, 'track.memory (frames)', '15');
            controls.CorrectDrift = addDropdown(layout, 'Correct Drift', {'on', 'off'}, 'off');

        case 'Segmentation'
            controls.GlobalBgCorr = addLabelField(layout, 'GlobalBgThr', '10');
            controls.ShowSegmentation = addDropdown(layout, 'ShowSegmentation:', {'on', 'off'}, 'on');
            controls.threshold = addLabelField(layout, 'Threshold:', '0.25');
            controls.diskDim = addLabelField(layout, 'Disk dim:', '2');

        case 'Phase'
            controls.dz = addLabelField(layout, 'PxSize z (µm):', '0.56');
            controls.NA = addLabelField(layout, 'NA detection:', '1.20');
            controls.NA_ill = addLabelField(layout, 'NA illumination:', '0.26');
            controls.n = addLabelField(layout, 'Refractive index:', '1.33');
            controls.lambda = addLabelField(layout, 'Central wavelength:', '0.58');
            controls.dlambda = addLabelField(layout, 'Spectrum bandwith:', '0.075');
            controls.alpha = addLabelField(layout, 'Alpha:', '3.15');
            controls.kzT = addLabelField(layout, 'axial cutoff:', '0.01');
            controls.mirrorX = addDropdown(layout, 'Mirror along x', {'true', 'false'}, 'false');
            controls.mirrorZ = addDropdown(layout, 'Mirror along z', {'true', 'false'}, 'true');
            controls.applyFourierMask = addDropdown(layout, 'denoising Fourier', {'true', 'false'}, 'true');

        case 'DDM'
            controls.ParticleSize = addLabelField(layout, 'Particle radius (nm):', '20');
            controls.ExpTime = addLabelField(layout, 'Exp time (s):', '0.03');
            controls.Wavelength = addLabelField(layout, 'Wavelength (nm):', '561');
            controls.NA = addLabelField(layout, 'NA:', '1.20');
            controls.Temp = addLabelField(layout, 'Temperature (K):', '296.15');
            controls.Qmin = addLabelField(layout, 'Qmin (µm^-1):', '3');
            controls.Qmax = addLabelField(layout, 'Qmax (µm^-1):', '10'); controls.Qmax.Tag = 'QmaxField';
            controls.Scanning = addDropdown(layout, 'Scanning reconstruction', {'on', 'off'}, 'off');
            controls.Scanning.Value = controls.Scanning.Value; 
            controls.CorrectBleaching = addDropdown(layout, 'CorrectBleaching', {'on', 'off'}, 'on');
            controls.CorrectBleaching.Value = controls.CorrectBleaching.Value;
            controls.AngularAnisotropy = addDropdown(layout, 'Angular anisotropy', {'on', 'off'}, 'off');
            controls.AngularAnisotropy.Value = controls.AngularAnisotropy.Value;
            controls.FitRDiff = addLabelField(layout, 'fitRange D:', '4');

            controls.Scanning.ValueChangedFcn = @(src,event) onScanningChange(src, layout, controlsName);
    end

    %% --- Nested scanning callback ---
    function onScanningChange(src, layout, controlsName)
        isOn = strcmp(src.Value, 'on');

        % Lock/unlock Qmax
        qmaxField = findobj(layout.Children, 'Tag', 'QmaxField');
        if ~isempty(qmaxField)
            qmaxField.Enable = ~isOn;
        end

        % Show/hide ROISize
        roiField = findobj(layout.Children, 'Tag', 'ROISizeField');
        if isOn
            if isempty(roiField)
                uilabel(layout, 'Text', 'Kernel Size (px):', 'HorizontalAlignment', 'right', 'FontSize', 11);
                roiEdit = uieditfield(layout, 'text', 'Value', '10', 'FontSize', 11, 'Tag', 'ROISizeField');
                state.(controlsName).ROISize = roiEdit; % store dynamically
            else
                roiField.Visible = 'on';
                state.(controlsName).ROISize = roiField;
            end
        else
            if ~isempty(roiField)
                roiField.Visible = 'off';
            end
            state.(controlsName).ROISize = [];
        end

        % Force RunMethod to 'run'
        parentFig = ancestor(layout, 'figure');
        runDropdown = findall(parentFig, 'Type', 'uidropdown');
        for d = 1:length(runDropdown)
            if any(ismember({'run','load'}, runDropdown(d).Items))
                runDropdown(d).Value = 'run';
            end
        end
    end
end

%% === UI Helper Functions ===
function field = addLabelField(layout, label, defaultValue)
    uilabel(layout, 'Text', label, 'HorizontalAlignment', 'right', 'FontSize', 11);
    field = uieditfield(layout, 'text', 'Value', defaultValue, 'FontSize', 11);
end

function dd = addDropdown(layout, label, items, defaultValue)
    uilabel(layout, 'Text', label, 'HorizontalAlignment', 'right', 'FontSize', 11);
    dd = uidropdown(layout, 'Items', items, 'Value', defaultValue, 'FontSize', 11);
end

function s = readControls(controls)
    s = struct();
    if isempty(fieldnames(controls))
        return;
    end
    fields = fieldnames(controls);
    for i = 1:length(fields)
        f = fields{i};
        val = controls.(f);

        % If it's a UI object with a Value property (dropdown, editfield, etc.)
        if isobject(val) && isprop(val, 'Value')
            raw = val.Value;
            % If it's an edit field (text entry), try to convert to numeric
            if isa(val, 'matlab.ui.control.EditField') || isa(val, 'matlab.ui.control.NumericEditField')
                if isempty(raw)
                    s.(f) = [];
                else
                    num = str2double(raw);
                    if isnan(num)
                        s.(f) = raw;      % keep string
                    else
                        s.(f) = num;      % numeric value
                    end
                end
            else
                % For dropdowns and other controls, just store the Value as-is
                s.(f) = raw;
            end

        % If it's a plain UI handle array / empty or already a value (char, numeric, etc.)
        else
            % If it's empty UI handle (like []), store empty
            if isempty(val)
                s.(f) = [];
            else
                % Plain values (numbers, strings, arrays) are stored directly
                s.(f) = val;
            end
        end
    end
end


end
