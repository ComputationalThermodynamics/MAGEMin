function exe = MAGEMin_exe(Computation)
% Returns the executable name of MAGEMin, which is slightly different on
% different OS's

if isunix
    switch  Computation.BinaryMethod
        case 'Default'
            exe = 'MAGEMin';    % use julia version
        case 'Local'
            exe = './MAGEMin';  % use local version
        otherwise
            error('Unknown method')
    end
else
    exe = 'MAGEMin.exe';    % windows
end
