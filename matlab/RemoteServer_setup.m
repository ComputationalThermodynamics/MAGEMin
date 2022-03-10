% This file is used if you want to perform MAGEMin calculations on a remote
% linux server
%
% You need several pre-requisites:
%
% 1) Make sure that MAGEMin runs fine on that machine in the directory that
%       will be used for the remote calculations. You will have to compile MAGEMin 
%       remotely beforehand (just as on your local machine, except that you may have 
%       to copy it to a different directory)   
%
% 2) Make sure that the required libraries (NLopt, lapacke) are linked correctly 
%       at the BEGINNING of your .bashrc file. Before a comments like: 
%           # If not running interactively, don't do anything
%
%    In my case, this looked like below. Note that I also added "/usr/bin/"
%    which is where "scp" is located on that machine (won't work otherwise)
%           # ~/.bashrc: executed by bash(1) for non-login shells.
%           # see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
%           # for examples
%     
%           # Add the required libraries for MAGEMin before the command below
%           export LD_LIBRARY_PATH=/opt/mpich3/lib/:/local/home/boris/Software/NLopt/install/lib
%           export PATH=PATH:/opt/mpich3/bin:/usr/bin/
%     
%           # If not running interactively, don't do anything
%           case $- in
%               *i*) ;;
%               *) return;;
%           esac
%
% 3) We will use ssh to remotely (and non-interactively) login to the
%       remote machine. Turns out that by default this sources
%       .bash_profile and not .bashrc. As a result, you may have to add
%       this to your .bash_profile file:
%
%           # Get the aliases and functions
%           if [ -f ~/.bashrc ]; then
%                   . ~/.bashrc
%           fi
%
% 4) You may want to setup keyless login to the remote machine; if not you
%       will be prompted multiple times for your password. Here a short
%       summary of how to set that up. 
%       We follow the steps as indicated here:
%
%          http://www.linuxproblem.org/art_9.html
%
%       Yet, importantly, we do NOT save it to the default file (id_rsa), 
%       but instead create a new file in ~/.ssh/id_rsa_MAGEMin
%       We copy that public key (id_rsa_MAGEMin.pub) over to the remote
%       machine as indicated
%  
%
% This has been tested on a mac with a remote linux machine; it'll likely
%   work the same on a linux machine, but is untested on windows

% Here your ssh login (as you do from the terminal):
remote_ssh_login     = 'boris@apollo.geo.uni-mainz.de';

% The remote directory where we'll run MAGEMin (which must be present
% there, so you must create it beforehand!)
remote_execution_dir = '~/MAGEMin_run';

% This is the file that contains the passwordless ssh key to login to the
% remote machine
ssh_key_file = '~/.ssh/id_rsa_MAGEMin';

% In case you want to test that this all works, run this in your terminal:
% ssh -i ~/.ssh/id_rsa_MAGEMin boris@apollo.geo.uni-mainz.de "cd ~/MAGEMin_run; mpiexec -n 1 ./MAGEMin"
% 
% and, provided that MAGEMin_input.dat exists locally, you can do:
% scp -i ~/.ssh/id_rsa_MAGEMin MAGEMin_input.dat boris@apollo.geo.uni-mainz.de:~/MAGEMin_run


