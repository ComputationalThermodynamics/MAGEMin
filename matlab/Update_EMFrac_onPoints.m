function PhaseData = Update_EMFrac_onPoints(PhaseData,XY_vec, newPoints,PrPts,Computation)
% This updates the EM fraction for new points based on previous calculations.
% Four different methods are available:
%       1) Use a previously trained neural network to estimate x-eos for
%           the main phases on the phase diagram. NOTE: this is under
%           development
%       2) Average from coarser mesh, which is known if we employ AMR
%       3) Use tha average of all SS models in the full domain
%       4) Average from a certain # of closests surrounding points



switch Computation.EMFrac_Method
    case 'Neural Network'
        % Use a previously trained NN to estimate the EM fraction for the
        % ph
        
        if length(PhaseData)==0
            PreviousSolutions = false;
        else
            PreviousSolutions = true;
        end
        
        PreviousSolutions = false;
      
        % Loop over points
        [EM_Fractions,SolidSolutions]  = EM_Fraction_From_NN(XY_vec(newPoints,2), XY_vec(newPoints,1));
        
        for i=1:length(newPoints)
            iPoint = newPoints(i);
        
            PhaseData{iPoint}.StartingValues_EMFraction.T = XY_vec(iPoint,1);
            PhaseData{iPoint}.StartingValues_EMFraction.P = XY_vec(iPoint,2);
            PhaseData{newPoints(i)}.StartingValues_EMFraction.EMFractions       = EM_Fractions{i};
            PhaseData{newPoints(i)}.StartingValues_EMFraction.SolutionModels    = SolidSolutions{i};
            
        end
        
    
    case 'Inherit from coarse mesh'
        
        if ~isempty(PrPts)
            for i=1:length(newPoints)
                SolidSolution_all   = [];
                EMFractions_all     = [];
                Prev = PrPts(i,:); Prev=Prev(Prev>0);
                for j=1:length(Prev)
                    SolidSolution_all   = [SolidSolution_all;   PhaseData{Prev(j)}.StableSolutions];
                    EMFractions_all     = [EMFractions_all;     PhaseData{Prev(j)}.EMFractions'];
                end

                % get unique Em_Fractions values for duplicate ones
                [SolidSolution, EMFractions] = ObtainUnique_EMFractions(SolidSolution_all, EMFractions_all);
                
                PhaseData{newPoints(i)}.StartingValues_EMFraction.EMFractions       = EMFractions;
                PhaseData{newPoints(i)}.StartingValues_EMFraction.SolutionModels    = SolidSolution;
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
            EMFractions_all             = [];
            for j=1:length(idx)
                SolidSolution_all   	= [SolidSolution_all;   PhaseData{idx(j)}.StableSolutions];
                EMFractions_all         = [EMFractions_all;     PhaseData{idx(j)}.EMFractions'];
            end
            
            
            % get unique Em_Fractions values for duplicate ones
            [SolidSolution, EMFractions] = ObtainUnique_EMFractions(SolidSolution_all, EMFractions_all);
            
            PhaseData{newPoints(i)}.StartingValues_EMFraction.EMFractions       = EMFractions;
            PhaseData{newPoints(i)}.StartingValues_EMFraction.SolutionModels    = SolidSolution;
        end
            
    case 'Average whole domain'
        % Make a list of all SS models
        SolidSolution_all 	= [];
        EMFractions_all     = [];
        for j=1:length(PhaseData)
            if isfield(PhaseData{j},'StableSolutions')
                % PhaseData can also include new points (just refined), for
                % which we obviously don't have a solution yet
                SolidSolution_all   = [SolidSolution_all; PhaseData{j}.StableSolutions];
                EMFractions_all     = [EMFractions_all;   PhaseData{j}.EMFractions'];
            end
        end
        
        % get unique Em_Fractions values for duplicate ones
        [SolidSolution, EMFractions] = ObtainUnique_EMFractions(SolidSolution_all, EMFractions_all);
        
        PhaseData{newPoints(i)}.StartingValues_EMFraction.EMFractions       = EMFractions;
        PhaseData{newPoints(i)}.StartingValues_EMFraction.SolutionModels    = SolidSolution;
        
    otherwise
        error('Method to compute xEOS not yet implemented')
end


function [SolidSolution, EMFraction] = ObtainUnique_EMFractions(SolidSolution_all, EMFractions_all)
% Get unique values from a whole list of solid solutions & 

SolidSolution           = [];
CompositionalVariable   = [];
[SolidSolution,~,iC] = unique(SolidSolution_all); % unique values
for iSol=1:length(SolidSolution)
    ind                         =   find(iC==iSol);
   
    % Average values
    EMFraction_local    =   reshape([EMFractions_all{ind}], [size(EMFractions_all{ind(1)},2) length(ind)])';
    average_EMFraction  =   mean(EMFraction_local,1);
    
    dist                =   EMFraction_local-repmat(average_EMFraction, [size(EMFraction_local,1),1]);
    dist                =   sum(dist.^2,2);
    max_dist            =   max(abs(dist));
    if max(abs(dist))<0.1
        % Seems safe to use the average value
        EMFraction{iSol}    =   average_EMFraction; % average all values
    else
        
        disp('EMFraction_local is: ')
        EMFraction_local
        disp(['We should really be using >1 compositions for: ',SolidSolution_all{ind(1)}, ' as the compositional distance is ',num2str(max_dist)])
        
        EMFraction{iSol}    =   EMFraction_local(end,:);
%         pause
    end
    
    
    
    
end