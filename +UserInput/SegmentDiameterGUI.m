function [val] = SegmentDiameterInput()
    % Create a UIFigure window
    fig = uifigure('Name', 'Enter Segment Diameter', 'Position', [500 500 300 150]);

    % Create a label
    lbl = uilabel(fig, ...
        'Position', [30 90 240 22], ...
        'Text', 'Enter Segment Diameter (in Px):');

    % Create an edit field for input
    edt = uieditfield(fig, 'numeric', ...
        'Position', [30 60 240 22]);

    % Create OK button
    btn = uibutton(fig, 'push', ...
        'Position', [100 20 100 30], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @(btn,event) okButtonCallback());

    % Wait for user input and button press before closing
    uiwait(fig);

    function okButtonCallback()
        val = edt.Value;
        if isempty(val) || ~isnumeric(val)
            uialert(fig, 'Please enter a valid numeric value.', 'Invalid Input');
            return
        end
        % Save input value to obj.info.SegmentDiameter
        uiresume(fig);
        close(fig);
    end
end

