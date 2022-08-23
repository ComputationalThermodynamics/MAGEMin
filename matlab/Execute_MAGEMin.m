function Execute_MAGEMin(Computation, VerboseLevel, n_points, Test, OxProp, sys_in)
% This executes MAGEMin, either on the local machine (on 1 or more MPI
% threads), or on a remote, more powerful, server.
% 
% Locally, there are two 


% Run simulation
NumRanks     =   Computation.NumRanks;
MPI_dir      =   Computation.MPI_dir;
RemoteServer =   Computation.RemoteServer;
MatlabOut    =   Computation.MatlabOut;


% Retrieve name of executable
exe = MAGEMin_exe(Computation);

% Setup the general execution command
command = [exe,' --out_matlab=',num2str(MatlabOut),' --Verb=',num2str(VerboseLevel),' --sys_in=',sys_in,' --File=MAGEMin_input.dat --n_points=',num2str(n_points)];

if ~isnan(Test)
    % employ a prededined test
    command = [command, ' --test=',num2str(Test)];
else
    % employ specified chemistry
    command = [command, ' --Bulk='];
    for iTable=1:size(OxProp,1)
        command = [command,num2str(table2array(OxProp(iTable,2))),','];
    end
end



if ~RemoteServer
    % In case we run it locally:    
    if NumRanks>1 
        if isempty(MPI_dir) || strcmp(Computation.BinaryMethod,'Default')==true
            command_MPI = ['mpiexec -n ',num2str(NumRanks),' '];
        else
            command_MPI = [MPI_dir,'/mpiexec -n ',num2str(NumRanks),' '];
        end
    else
        command_MPI = '';
    end
    command = [command_MPI, command];
    
    % In case we do NOT use the locally compiled MAGEMin, but instead the
    % automatically downloaded MAGEMin binaries, we need to load the correct
    % directories first:
     
    disp(command)
    if Computation.Julia_MAGEMin_binary==true
        command = add_dynamic_libs(command, Computation)
    end
    if isunix
        system('killall MAGEMin 2>&1')
    end
    system(command);
    
    
    
else
    
   % We want to run it remotely, we need to 
   %  1) copy the MAGEMin_input.dat to the remote machine, 
   %  2) run MAGEMin on the requested # of MPI cores
   %  3) copy the results back
   RemoteServer_setup % load all data (server name etc.)
   
   % 1) Copy local file
   start = tic;
   copy_to_remote = ['scp -q -i ',ssh_key_file,' MAGEMin_input.dat ',remote_ssh_login,':',remote_execution_dir];
   system(copy_to_remote);
   tot_time_step1 = toc(start);
   disp(['Copy local file took ',num2str(tot_time_step1),'s'])
   
   
   % 2) Perform remote calculation
   start = tic;
   command_zip      = ['; zip -q -r output OUTPUT/'];   % zip folder remotely
   command_remote   = ['ssh -i ',ssh_key_file,' ',remote_ssh_login ,' "cd ',remote_execution_dir,'; mpiexec -n ',num2str(NumRanks),' ', command,command_zip, '"']
   command_remote   = [command_remote, command_zip];
   disp(command_remote)
   system(command_remote);
   disp(command_remote)
   tot_time_step2 = toc(start);
   disp(['Performing remote computations took ',num2str(tot_time_step2),'s'])
   
   
   % 3) Copy remote back to local:
   start = tic;
   copy_from_remote = ['scp -q -i ',ssh_key_file,' ',remote_ssh_login,':',remote_execution_dir,'/output.zip .'];
   disp(copy_from_remote)
   system(copy_from_remote);
   system('unzip -q -o output.zip');
   tot_time_step3 = toc(start);
   disp(['Retrieving remote computational results took ',num2str(tot_time_step3),'s'])
   
   
end

end

