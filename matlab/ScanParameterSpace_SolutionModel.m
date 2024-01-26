function [Local_best, Local_cluster, Local] = ScanParameterSpace_SolutionModel(P,T,Gamma,SS, n_pc, MPI_dir, MPI_cores, SampleG)
% This scans the local parameter space of a solution model and returns the
% best-fit model
%

% Prepare string
str = sprintf('./MAGEMin --Mode=2 --Temp=%f --Pres=%f --n_pc=%i --Phase=%s --Gam=',T, P, n_pc, SS);


for i=1:length(Gamma)
    str = [str, sprintf('%4.12f,',Gamma(i));];
end


if SampleG
    str = [str, ' --maxeval=1 '];        % in case we want to sample the G-surface rather than do optimizations
end

% Perform computation
if isempty(MPI_dir);
    disp(str)
    system(str);
else
    str = [MPI_dir,'/mpiexec -np ',num2str(MPI_cores),' ',str];
    disp(str)
    system(str);
end


% Load data
Local               =   Read_LocalOptimizations_MAGEMin('./OUTPUT');
Local.SolutionModel =   SS;     % add solution model

% Remove points where the SF are not OK
Local.all           =   Local;
ind_SF_OK           =   find(Local.SF_ok==1);
Local.df            =   Local.df(ind_SF_OK);
Local.Prop_start    =   Local.Prop_start(ind_SF_OK,:);
Local.xEOS_start    =   Local.xEOS_start(ind_SF_OK,:);
Local.xEOS_end      =   Local.xEOS_end(ind_SF_OK,:);
Local.Prop_end      =   Local.Prop_end(ind_SF_OK,:);
Local.Status        =   Local.Status(ind_SF_OK,:);
Local.SF_ok         =   Local.SF_ok(ind_SF_OK,:);


% 1) 'Best' model (aka, one with lowest df)
[~,ind_best]                =   min(Local.df);
Local_best.SolutionModel	=   SS;
Local_best.df               =   Local.df(ind_best);
Local_best.Prop_start       =   Local.Prop_start(ind_best,:);
Local_best.Prop_end         =   Local.Prop_end(ind_best,:);
Local_best.xEOS_start       =   Local.xEOS_start(ind_best,:);
Local_best.xEOS_end         =   Local.xEOS_end(ind_best,:);
Local_best.Status           =   Local.Status(ind_best,:);
Local_best.Prop_end         =   Local.Prop_end(ind_best,:);


Local_cluster=[];
if nargout>1
    % 2) cluster the data and determine optimal # of clusters
    
    % We need to find the optimal number of clusters, which is in fact the one
    % that gives the overall lowest variance in the data. In many cases we have
    % 1-2 clusters, which is why we check up to 15 here
    for iCluster=1:15;
        [clusters] = kmeans(Local.xEOS_end',iCluster);
        for i=1:max(clusters)
            ind         =   find(clusters==i);
            variance(i) =   sum(abs(var(Local.xEOS_end(ind,:),1)));
        end
        VarCluster(iCluster) = sum(variance);
        
    end
    [~,optNumClusters] = min(VarCluster);    % optimal # of clusters
    
    
    [clusters] = kmeans(Local.xEOS_end',optNumClusters);
    
    Local_cluster.SolutionModel	=   SS;
    
    for i=1:optNumClusters
        ind                                 =   find(clusters==i);
        Local_cluster.df(i)                 =   mean(Local.df(ind),1);
        Local_cluster.Prop_start(i,:)       =   mean(Local.Prop_start(ind,:),1);
        Local_cluster.Prop_end(i,:)         =   mean(Local.Prop_end(ind,:),1);
        Local_cluster.xEOS_start(i,:)       =   mean(Local.xEOS_start(ind,:),1);
        Local_cluster.xEOS_end(i,:)         =   mean(Local.xEOS_end(ind,:),1);
        Local_cluster.Prop_end(i,:)         =   mean(Local.Prop_end(ind,:),1);
        
        Local_cluster.Status(i,:)           =   mean(Local.Status(ind,:),1);
    end
    
end
