function [Data_new, Data_TC_in] = Compute_MAGEMin_SinglePoint(varargin)
% This performs a MAGEMin computation for a single point
%
% It is essentially used to debug the code, and test the effect of different 
% Gamma's, solution models and computational options on the MAGEMin solution
% 

Data                =   varargin{1};
UseGammaEstimation  =   varargin{2};
Use_xEOS            =   varargin{3};
VerboseLevel        =   varargin{4};

if UseGammaEstimation
    disp('Using Gamma estimation')
    if nargin>4
        % optionally provide Gamma to override the one in the Data
        % structure
        Gamma           =   varargin{5};
        if length(Gamma)==10
           Gamma(end+1)=0;      % in case we use TC input w/out H2O
        end
        if ~isempty(Gamma)
            Data.Gamma = Gamma;
        end
    end
else
    disp('Not using Gamma estimation')
end    

if Use_xEOS
    disp('Using xEOS estimation')
    if nargin>5
         SolutionModels =   varargin{6};
         Data.StartingValues_xEOS.SolutionModels = SolutionModels;
    end
    
    str = ' Solution Models considered : ';
    for i=1:length( Data.StartingValues_xEOS.SolutionModels)
        str = [str,' ',Data.StartingValues_xEOS.SolutionModels{i}];
    end
    disp(str)
    disp('')
    
    % We can either supply the xEOS as input for all solution models
    % considered above, or use the values that exist in the structure
    % already
    if nargin>6
        CompositionalVar =   varargin{7};
        Data.StartingValues_xEOS.CompositionalVar = CompositionalVar;
    end
    if ~isfield(Data.StartingValues_xEOS,'CompositionalVar')
        Data.StartingValues_xEOS.CompositionalVar = [];
    end
    Data.StartingValues_xEOS.StableFractions = Data.StableFractions;

    % Check if size(CompositionalVar)==size(SolutionModels)
    % If not, we expect the last ones to be 'addional'
    if length(Data.StartingValues_xEOS.SolutionModels)> length(Data.StartingValues_xEOS.CompositionalVar)
        
        for iSS=length(Data.StartingValues_xEOS.CompositionalVar)+1:length(Data.StartingValues_xEOS.SolutionModels)
            SS                          = Data.StartingValues_xEOS.SolutionModels{iSS};
            [VarName,  InitialValues]   = SolidSolution_VarNames(SS);
            Data.StartingValues_xEOS.CompositionalVar{iSS} = InitialValues;
            
            disp([' Added starting initial guess for ',Data.StartingValues_xEOS.SolutionModels{iSS}])

        end
    end
    if length(Data.StartingValues_xEOS.StableFractions)< length(Data.StartingValues_xEOS.SolutionModels)
        for iSS=length(Data.StartingValues_xEOS.StableFractions)+1:length(Data.StartingValues_xEOS.SolutionModels)
             Data.StartingValues_xEOS.StableFractions(iSS)  = 1e-3;
        end
    end
    
    

else
    disp('Not using xEOS estimation')
end    

% Default values
newPoints           =   1;
XY_vec              =   [Data.T, Data.P];
Chemistry           =   Data.Chemistry;
dlg.CancelRequested =   false;
ComputeAllPoints    =   logical(0);

Computation.Use_xEOS = Use_xEOS;
Computation.NumRanks = 1;
Computation.MPI_dir  = '';
Computation.MinPhaseFraction = 0;

PhaseData_in{1} = Data;

% Perform computation and retrieve result 
[Data_new, XY_vec, FailedSimulations]  = 	PerformMAGEMin_Simulation(PhaseData_in, newPoints, XY_vec, VerboseLevel, Chemistry, dlg, ComputeAllPoints, UseGammaEstimation, Computation);

Data_new = Data_new{1};
Data_new.Chemistry = Data.Chemistry;



% Also output the "Data" structure in a format that can be used directly by TC
if Use_xEOS
    Data_TC_in                  =   Data;
    Data_TC_in.StableFractions  =   Data.StartingValues_xEOS.StableFractions;
    Data_TC_in.StableSolutions  =   Data.StartingValues_xEOS.SolutionModels;
    Data_TC_in.CompositionalVar =   Data.StartingValues_xEOS.CompositionalVar;
else
    Data_TC_in=[];
end







