function [TrainingData] = TrainingData_SolutionModel(TrainingData, SS, MPI_dir, MPI_cores, n_pc, varargin)
% This generates training data for a SS model, for use with a neural
% network. 
%
% There are 2 modes:
%   1) Predefined Gamma, T,P 
%         [TrainingData]    = TrainingData_SolutionModel(TrainingData, SS, MPI_dir, MPI_cores, n_pc, T, P, Gamma)
%           TrainingData    =   Previous data (optional), or empty
%           SS              =   solid solution (e.g. 'spl')
%           MPI_dir         =   directory that contains the MPI executable (empty if no MPI)
%           MPI_cores       =   # of cores (optional)
%           n_pc            =   # points/dimension
%           T               =   vector with T points
%           P               =   vector with P points
%           Gamma           =   11xn vector with Gamma points
%
%   2) Create random selection
%           [TrainingData]  =   TrainingData_SolutionModel(TrainingData, SS, MPI_dir, MPI_cores, n_pc, nPoints, Bounds_T, Bounds_P, Bounds_Gamma)
%           TrainingData    =   Previous data (optional), or empty
%           SS              =   solid solution (e.g. 'spl')
%           MPI_dir         =   directory that contains the MPI executable (empty if no MPI)
%           MPI_cores       =   # of cores (optional)
%           nPoints         =   # of points
%           Bounds_T        =   lower and upper bounds of T
%           Bounds_P        =   lower and upper bounds of P
%           Bounds_Gamma    =   lower and upper bounds of Gamm


if nargin==8
   T_vec            =   varargin{1};
   P_vec            =   varargin{2};
   Gamma_vec        =   varargin{3};
   
   nPoints          =   length(T_vec);
   
   RandomPoints     =   false;
   
elseif nargin==9
   nPoints          =   varargin{1};
   Bounds_T         =   varargin{2};
   Bounds_P         =   varargin{3};
   Bounds_Gamma     =   varargin{4};
   
   RandomPoints     =   true;
   
end

if ~isempty(TrainingData)
    numTraining = length(TrainingData.Gamma);
else
    numTraining = 0;
end

% Create directory based on SS
curdir = pwd;
dirname = [SS,'_training'];
mkdir(dirname)
cd(dirname)
system('cp ../MAGEMin .')


if ismac
    fig = uifigure;
    progressDialog = uiprogressdlg(fig,'Title',['Training points for ',SS] );
end
for iPoint=1:nPoints
    i       = numTraining + iPoint
    
    if RandomPoints
        % Take random sample
        T       = rand(1)*diff(Bounds_T)        + Bounds_T(1);
        P       = rand(1)*diff(Bounds_P)        + Bounds_P(1);
        Gamma   = rand(size(Bounds_Gamma)).*diff(Bounds_Gamma,[],2)  + Bounds_Gamma(:,1);
        
    else
        T       = T_vec(iPoint);
        P       = P_vec(iPoint);
        Gamma   = Gamma_vec(iPoint,:);
    end
    
    % Do calculation
    [Local_best]    = ScanParameterSpace_SolutionModel(P,T,Gamma,SS, n_pc, MPI_dir , MPI_cores, logical(0));
    
    % Store data
    TrainingData.Gamma(i,:)     = Gamma(:)';
    TrainingData.T(i,1)         = T;
    TrainingData.P(i,1)         = P;
    TrainingData.df(i,1)        = Local_best.df;
    TrainingData.xEOS_end(i,:)  = Local_best.xEOS_end;
    TrainingData.Prop_end(i,:)  = Local_best.Prop_end;
    
    str = ['Processing the points ; ',num2str(iPoint/nPoints*100),'% Total:  ',num2str(nPoints),' points'];
    if ismac
        % update progress bar 
        progressDialog.Value = iPoint/nPoints;
        progressDialog.Message = str;
        drawnow
        
    else
        disp(str);

    end
    
    if mod(iPoint,10)==0
    
       % Save intermediate data to not loose anything if something goes
       % wrong
       save(['Training_debug_',SS],'TrainingData');
       
    end
    
    

end

if ismac
    delete(fig)
end

cd(curdir)





