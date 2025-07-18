function [info, info1, info2, file] = infoGUI(file)
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
    


    state.Ext = addDropdown(col1, 'Ext:', {'.ome.tif', '.tif', '.his', '.mpg', '.spe', '.lif'}, '.his');
    state.Type = addDropdown(col1, 'Type:', {'normal', 'transmission'}, 'normal');
    state.RunMethod = addDropdown(col1, 'Run Method:', {'run', 'load'}, 'run');
    state.calibrate = addDropdown(col1, 'calibrate', {'true', 'false'}, 'false');
    state.drawROI = addDropdown(col1, 'draw ROI', {'on', 'off'}, 'on');

    state.Dimension = addDropdown(col1, 'Dimension:', {'2D', '3D'}, '2D');
    state.multiModal = addDropdown(col1, 'multiModal:', {'on', 'off'}, 'on');
    ch1Dropdown = addDropdown(col1, 'Channel 1:', ...
        {'Translational Tracking', 'Rotational Tracking', 'Segmentation', 'Phase'}, ...
        'Translational Tracking');
    ch2Dropdown = addDropdown(col1, 'Channel 2:', ...
        {'Translational Tracking', 'Rotational Tracking', 'Segmentation', 'Phase'}, ...
        'Segmentation');

    state.PxSize = addLabelField(col1, 'PxSize (nm):', '95');
    state.FWHM = addLabelField(col1, 'FWHM (px):', '3');
    
    state.Frame2Load = addLabelField(col1, 'Frame2Load:', 'all');
    state.TestFrame = addLabelField(col1, 'Test Frame:', '1');

    state.Rotational = addDropdown(col1, 'Rotational:', {'on', 'off'}, 'on');
    state.Bipyramid = addLabelField(col1, 'Bipyramid (nm):', '[184 92]');
    state.RotationalCalib = addDropdown(col1, 'Rotational Calib:', {'on', 'off'}, 'off');
    state.RadTime = addLabelField(col1, 'Rad Time (°/s):', '25');
    

    %% Channel Panels
    channel1Panel = uipanel(gl, 'Title', 'Channel 1 Setup');
    channel1Layout = uigridlayout(channel1Panel, [14, 2]);
    channel1Layout.RowHeight = {'fit'};
    state.channel1Controls = [];

    channel2Panel = uipanel(gl, 'Title', 'Channel 2 Setup');
    channel2Layout = uigridlayout(channel2Panel, [14, 2]);
    channel2Layout.RowHeight = {'fit'};
    state.channel2Controls = [];

    channel1Layout.RowHeight = repmat({30}, 1, 14);
    channel2Layout.RowHeight = repmat({30}, 1, 14);

    %% OK Button
    btn = uibutton(gl, 'Text', 'OK', ...
        'ButtonPushedFcn', @(btn, event) onOK());
    btn.Layout.Row = 2;
    btn.Layout.Column = [1 3];

    %% Setup logic
    ch1Dropdown.ValueChangedFcn = @(src, event) updateChannel(1, src.Value);
    ch2Dropdown.ValueChangedFcn = @(src, event) updateChannel(2, src.Value);
    state.RotationalCalib.ValueChangedFcn = @(src, event) toggleRadTime();

    % Initial setup
    updateChannel(1, ch1Dropdown.Value);
    updateChannel(2, ch2Dropdown.Value);
    toggleRadTime();
    toggleBipyramid();
    toggleRotationalControls();

    % Wait for user to press OK
    uiwait(fig);

    %% --- Nested functions ---

    function updateChannel(channelNum, mode)
        if channelNum == 1
            delete(channel1Layout.Children);
            state.channel1Controls = addChannelControls(channel1Layout, mode);
        else
            delete(channel2Layout.Children);
            state.channel2Controls = addChannelControls(channel2Layout, mode);
        end
        toggleBipyramid();
        toggleRotationalControls();
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
            state.Rotational.Enable = 'off';  % Fixed to on
            state.RotationalCalib.Enable = 'on';  % Let user change it
        else
            state.Rotational.Value = 'off';
            state.Rotational.Enable = 'off';  % Fixed to 0
            state.RotationalCalib.Value = 'off';
            state.RotationalCalib.Enable = 'off';
        end
    
        toggleRadTime();  % Always refresh dependent RadTime state
    end

    function toggleSegmentDiameter()
        ctrlSets = {state.channel1Controls, state.channel2Controls};
        for i = 1:2
            ctrlSet = ctrlSets{i};
            if isfield(ctrlSet, 'CheckSize') && isfield(ctrlSet, 'SegmentDiameter')
                if strcmp(ctrlSet.CheckSize.Value, 'on')
                    ctrlSet.SegmentDiameter.Value = '[]';
                    ctrlSet.SegmentDiameter.Enable = 'off';
                else
                    ctrlSet.SegmentDiameter.Enable = 'on';
                end
            end
        end
    end

    function toggleTestFrame()
        ctrlSets = {state.channel1Controls, state.channel2Controls};
        for i = 1:2
            ctrlSet = ctrlSets{i};
            if isfield(ctrlSet, 'ShowSegmentation') && isfield(ctrlSet, 'CheckSize') && isfield(ctrlSet, 'TestFrame')
                checkOn = strcmp(ctrlSet.CheckSize.Value, 'on');
                showSegOn = strcmp(ctrlSet.ShowSegmentation.Value, 'on');
    
                if checkOn || showSegOn
                    ctrlSet.TestFrame.Enable = 'on';
                else
                    ctrlSet.TestFrame.Value = '[]';
                    ctrlSet.TestFrame.Enable = 'off';
                end
            end
        end
    end


    function onOK()
        % Collect info structure
        info.Dimension = state.Dimension.Value;
        info.PxSize = str2double(state.PxSize.Value);
        info.type = state.Type.Value;
        info.runMethod = state.RunMethod.Value;
        if strcmp(state.multiModal.Value, 'on')
            info.multiModal = 1;
        else 
            info.multiModal = 0;
        end
        if strcmp(state.drawROI.Value, 'on')
            info.drawROI = 1;
        else 
            info.drawROI = 0;
        end

        if strcmp(state.calibrate.Value, 'true')
            info.calibrate = 1;
        else
            info.calibrate = 0;
        end
        if strcmp(state.Frame2Load.Value, 'all')
            info.frame2Load = state.Frame2Load.Value;
        else
            info.frame2Load = evalin('base', state.Frame2Load.Value);  %#ok<EVLC>
        end
        info.Channel1 = ch1Dropdown.Value;
        info.Channel2 = ch2Dropdown.Value;
        info.FWHM = str2double(state.FWHM.Value);
        if strcmp(state.Rotational.Value, 'on')
            info.rotational = 1;
        else 
            info.rotational = 0;
        end
        if strcmp(state.Rotational.Value, 'on')
            info.rotationalCalib = 1;
        else 
            info.rotationalCalib = 0;
        end
        info.TestFrame = str2double(state.TestFrame.Value);
        info.RadTime = str2double(state.RadTime.Value);
        info.Bipyramid = str2num(state.Bipyramid.Value); %#ok<ST2NM>
        file.ext = state.Ext.Value;

        % Substructures
        info1 = readControls(state.channel1Controls);
        info2 = readControls(state.channel2Controls);

        if contains(info.Channel1, 'Tracking')
            info1.detectParam.delta = info1.delta;
            info1.detectParam.chi2 = info1.chi2;
            info1.detectParam.consThresh = info1.consThresh;
            info1.trackParam.radius = info1.track_radius;
            info1.trackParam.memory = info1.track_memory;

            info1 = rmfield(info1, {'delta', 'chi2', 'consThresh', 'track_radius', 'track_memory'});
        elseif contains(info.Channel1, 'Phase')
            info1.optics.dz = info1.dz;
            info1.optics.NA = info1.NA;
            info1.optics.NA_ill = info1.NA_ill;
            info1.optics.n = info1.n;
            info1.optics.lambda = info1.lambda;
            info1.optics.dlambda = info1.dlambda;
            info1.optics.alpha = info1.alpha;
            info1.optics.kzT = info1.kzT;

            info1.proc.mirrorX = info1.mirrorX;
            info1.proc.mirrorZ = info1.mirrorZ;
            info1.proc.applyFourierMask = info1.applyFourierMask;

            info1 = rmfield(info1, {'dz', 'NA', 'NA_ill', 'n', 'lambda', 'dlambda', 'alpha', 'kzT',...
                        'mirrorX', 'mirrorZ', 'applyFourierMask'});
        end

        if contains(info.Channel2, 'Tracking')
            info2.detectParam.delta = info2.delta;
            info2.detectParam.chi2 = info2.chi2;
            info2.detectParam.consThresh = info2.consThresh;
            info2.trackParam.radius = info2.track_radius;
            info2.trackParam.memory = info2.track_memory;

            info2 = rmfield(info2, {'delta', 'chi2', 'consThresh', 'track_radius', 'track_memory'}); 
        elseif contains(info.Channel1, 'Phase')
            info2.optics.dz = info2.dz;
            info2.optics.NA = info2.NA;
            info2.optics.NA_ill = info2.NA_ill;
            info2.optics.n = info2.n;
            info2.optics.lambda = info2.lambda;
            info2.optics.dlambda = info2.dlambda;
            info2.optics.alpha = info2.alpha;
            info2.optics.kzT = info2.kzT;

            info2.proc.mirrorX = info2.mirrorX;
            info2.proc.mirrorZ = info2.mirrorZ;
            info2.proc.applyFourierMask = info2.applyFourierMask;

            info2 = rmfield(info2, {'dz', 'NA', 'NA_ill', 'n', 'lambda', 'dlambda', 'alpha', 'kzT',...
                        'mirrorX', 'mirrorZ', 'applyFourierMask'});
        end

        
        % Resume and close
        uiresume(fig);
        delete(fig);
    end
function controls = addChannelControls(layout, type)
    controls = struct();
    switch type
        case {'Translational Tracking', 'Rotational Tracking'}
            controls.fitMethod = addDropdown(layout, 'fitMethod', {'Phasor', 'Gauss'}, 'Phasor');
            controls.zMethod = addDropdown(layout, 'zMethod', {'Intensity', '3DFit', 'PSFE'}, 'Intensity');
            controls.detectionMethod = addDropdown(layout, 'detectionMethod', {'Intensity', 'MaxLR'}, 'MaxLR');
            controls.IntCorr = addDropdown(layout, 'Intensity correction', {'on', 'off'}, 'off');
            controls.euDist = addLabelField(layout, 'euDist (nm)', '1000');
            controls.delta = addLabelField(layout, 'delta', '6');
            controls.chi2 = addLabelField(layout, 'chi2', '50');
            controls.consThresh = addLabelField(layout, 'consThresh', '4');
            controls.track_radius = addLabelField(layout, 'track.radius (nm)', '2500');
            controls.track_memory = addLabelField(layout, 'track.memory (frames)', '3');

        case 'Segmentation'
            controls.GlobalBgCorr = addLabelField(layout, ...
                'GlobalBgThr', '10');
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
            controls.mirrorZ = addDropdown(layout, 'Mirror along z', {'true', 'false'}, 'true');             % mirror the input stack along Z
            controls.applyFourierMask = addDropdown(layout, 'denoising Fourier', {'true', 'false'}, 'true'); 
    end
end
end

%% === UI Helper Functions ===

function field = addLabelField(layout, label, defaultValue)
    uilabel(layout, 'Text', label, ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 11);
    field = uieditfield(layout, 'text', ...
        'Value', defaultValue, ...
        'FontSize', 11);
end

function dd = addDropdown(layout, label, items, defaultValue)
    uilabel(layout, 'Text', label, ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 11);
    dd = uidropdown(layout, ...
        'Items', items, ...
        'Value', defaultValue, ...
        'FontSize', 11);
end

function s = readControls(controls)
    s = struct();
    fields = fieldnames(controls);
    for i = 1:length(fields)
        f = fields{i};
        val = controls.(f);
        if isa(val, 'uidropdown')
            s.(f) = val.Value;
        else
            valStr = val.Value;
            valNum = str2num(valStr); %#ok<ST2NM>
            if isempty(valNum)
                s.(f) = valStr;
            else
                s.(f) = valNum;
            end
        end
    end
end