function set_listener_enabled(hListener, value)
%SET_LISTENER_ENABLED Helper function to set 'Enabled' property for listeners.
%
% Usage:
%
%   SET_LISTENER_ENABLED(H, VALUE)
%
% Inputs:
%
%   H <1x1 listener handle>
%     - Handle to listener object to set 'Enabled' property for
%
%   VALUE <numeric, logical, or character array>
%     - Numeric: 0 or 1
%     - Logical: false or true
%     - Character array: 'off' or 'on'
%
% Description:
%
%   At some point, the input for setting a listener's 'Enabled' property
%   changed from a numeric/logical (0,1) to a character array ('on','off').
%   This helper function abtracts away that difference by allowing a
%   numeric, logical, or character array set value to be specified and then
%   setting the listener 'Enabled' property with the expected type.

% Author: Todd Baxter
% Date: 07-03-2018

    % determine if 'Enabled' needs to be set as a string or logical
    logicalSetType = true;
    if ischar(hListener.Enabled)
        logicalSetType = false;
    end
    % handle numeric, string, and logical set value types
    if ischar(value)
        value = lower(value);
    end
    switch value
        case {0,'off'}
            if logicalSetType
                hListener.Enabled = false;
            else
                hListener.Enabled = 'off';
            end
        case {1,'on'}
            if logicalSetType
                hListener.Enabled = true;
            else
                hListener.Enabled = 'on';
            end
        otherwise
            error('SET_LISTENER_ENABLED:invalidSetValue', 'Invalid listener ''Enabled'' set value');
    end
end
