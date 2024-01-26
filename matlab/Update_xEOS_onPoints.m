function PhaseData = Update_xEOS_onPoints(PhaseData,XY_vec, newPoints,PrPts,Computation)
% This updates the xEOS for new points based on previous calculations.
% Four different methods are available:
%       1) Use a previously trained neural network to estimate x-eos for
%           the main phases on the phase diagram. NOTE: this is under
%           development
%       2) Average from coarser mesh, which is known if we employ AMR
%       3) Use tha average of all SS models in the full domain
%       4) Average from a certain # of closests surrounding points



switch Computation.EOS_Method
    case 'Neural Network'
        % Use a previously trained NN to estimate x-eos for the most
        % important phases here
        
        % Load networks (for now we employ the deep learning TB; can be changed to scripts later)
        if length(PhaseData)==0
            PreviousSolutions = false;
        else
            PreviousSolutions = true;
        end
        
        % Use the NN to evaluate the data for the points as a function of T,P (and later chemistry)
        SolutionModels = {'spn';};
%         SolutionModels = {'spn';'ol';'cpx'; 'opx'; 'pli';'g'; 'liq'};
        
        for iSS = 1:length(SolutionModels)
           CompositionalVar_NN{iSS}  = xEOS_From_NN(SolutionModels{iSS},XY_vec(newPoints,:));
        end
        
        % Loop over points
        for i=1:length(newPoints)
            iPoint = newPoints(i);
            if PreviousSolutions 
                Data        = PhaseData{newPoints(i)};
            else
                Data        = [];
            end
            
            Data.T               = XY_vec(iPoint,1);
            Data.P               = XY_vec(iPoint,2);
            
            % Specify stable solutions from NN
            Data.StableSolutions = SolutionModels;
            
            if ~isfield(Data,'Gamma')
                PhaseData{iPoint}.Gamma= zeros(1,11);
            end
            for iSS = 1:length(SolutionModels)
                Data.CompositionalVar{iSS} = CompositionalVar_NN{iSS}(i,:);
            end
            
            PhaseData{iPoint}.StartingValues_xEOS.T = Data.T;
            PhaseData{iPoint}.StartingValues_xEOS.P = Data.P;
            PhaseData{iPoint}.StartingValues_xEOS.CompositionalVar =  Data.CompositionalVar;
            PhaseData{iPoint}.StartingValues_xEOS.SolutionModels   =  Data.StableSolutions;
            
        end
        
    
    case 'Inherit from coarse mesh'
        
        if ~isempty(PrPts)
            for i=1:length(newPoints)
                SolidSolution_all           = [];
                CompositionalVariable_all   = [];
                Prev = PrPts(i,:); Prev=Prev(Prev>0);
                for j=1:length(Prev)
                    SolidSolution_all         = [SolidSolution_all; PhaseData{Prev(j)}.StableSolutions];
                    CompositionalVariable_all = [CompositionalVariable_all;  PhaseData{Prev(j)}.CompositionalVar'];
                end
                


                % get unique xEOS values for duplicate ones
                [SolidSolution, CompositionalVariable] = ObtainUnique_xEOS(SolidSolution_all, CompositionalVariable_all);
                PhaseData{newPoints(i)}.StartingValues_xEOS.CompositionalVar = CompositionalVariable;
                PhaseData{newPoints(i)}.StartingValues_xEOS.SolutionModels   = SolidSolution;
            end
        end
        
    case 'Average surrounding points'
     
        TP_norm = max(XY_vec)-min(XY_vec);
        
        oldPoints = 1:length(PhaseData);
        rmPoints = [];
        for i=1:length(PhaseData)
            if ~isfield(PhaseData{i},'StableSolutions')
                rmPoints(end+1) = i;
            end
        end
        oldPoints(rmPoints)=[];

        
        % Set this to all newPoints
        for i=1:length(newPoints)
            id = newPoints(i);
            
            if isfield(Computation,'num_Neighbours')
                num_Neighbours  =   Computation.num_Neighbours;
            else
                num_Neighbours  =   1;
            end
            
            % Get closest points to this one
            idx = [];
            r   = 0.05;
            while length(idx)<num_Neighbours
                [idx]   =   nearestneighbour([XY_vec(id,:)./TP_norm]', [XY_vec(oldPoints(:),:)./TP_norm]','n',num_Neighbours,'r',r);
                r       =   r+0.05;     % increase the search radius until we have sufficient points
            end
            
            SolidSolution_all           = [];
            CompositionalVariable_all   = [];
            for j=1:length(idx)
                SolidSolution_all         = [SolidSolution_all; PhaseData{idx(j)}.StableSolutions];
                CompositionalVariable_all = [CompositionalVariable_all;  PhaseData{idx(j)}.CompositionalVar'];
            end
            
            
            % get unique xEOS values for duplicate ones
            [SolidSolution, CompositionalVariable] = ObtainUnique_xEOS(SolidSolution_all, CompositionalVariable_all);
         
%             % testing: don't do this for liq 
%             ind_Liq                         = find(strcmp(SolidSolution,'liq'));
%             CompositionalVariable(ind_Liq)  = [];
%             SolidSolution(ind_Liq)          = [];
% 
            PhaseData{newPoints(i)}.StartingValues_xEOS.CompositionalVar = CompositionalVariable;
            PhaseData{newPoints(i)}.StartingValues_xEOS.SolutionModels   = SolidSolution;
        end
            
    case 'Average whole domain'
        % Make a list of all SS models
        SolidSolution_all           = [];
        CompositionalVariable_all   = [];
        for j=1:length(PhaseData)
            if isfield(PhaseData{j},'StableSolutions')
                % PhaseData can also include new points (just refined), for
                % which we obviously don't have a solution yet
                SolidSolution_all         = [SolidSolution_all; PhaseData{j}.StableSolutions];
                CompositionalVariable_all = [CompositionalVariable_all;   PhaseData{j}.CompositionalVar'];
            end
        end
        
        % get unique xEOS values for duplicate ones
        [SolidSolution, CompositionalVariable] = ObtainUnique_xEOS(SolidSolution_all, CompositionalVariable_all);
        
        % Set this to all newPoints  
        for i=1:length(newPoints)
            PhaseData{newPoints(i)}.StartingValues_xEOS.CompositionalVar = CompositionalVariable;
            PhaseData{newPoints(i)}.StartingValues_xEOS.SolutionModels   = SolidSolution;
        end
        
    otherwise
        error('Method to compute xEOS not yet implemented')
end


function [SolidSolution, CompositionalVariable] = ObtainUnique_xEOS(SolidSolution_all, CompositionalVariable_all)
% Get unique values from a whole list of solid solutions & 

SolidSolution           = [];
CompositionalVariable   = [];
[SolidSolution,~,iC] = unique(SolidSolution_all); % unique values
for iSol=1:length(SolidSolution)
    ind                         =   find(iC==iSol);
    CompositionalVariable{iSol} =   reshape([CompositionalVariable_all{ind}], [size(CompositionalVariable_all{ind(1)},2) length(ind)])';
    CompositionalVariable{iSol} =   mean(CompositionalVariable{iSol},1); % average all values
end