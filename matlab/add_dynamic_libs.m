function command = add_dynamic_libs(command, Computation)
% Adds dynamic libraries if required
if isfield(Computation,'path_bin')
    
    if isunix
        command_dyn = ['export PATH=', Computation.path_bin, ' ;  export DYLD_LIBRARY_PATH=',Computation.path_dylib ,'; ' ];
        
    else
        % windows system
        command_dyn = ['set PATH=', Computation.path_bin, ';',Computation.path_dylib,';%PATH% && ' ];
        
    end
    command = [command_dyn, command];
end

end