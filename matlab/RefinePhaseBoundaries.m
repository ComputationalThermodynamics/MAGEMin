function [refine_Elements, React_output, PhaseData] = RefinePhaseBoundaries(PhaseData, TP_vec, elements, NumAssemblage)
% This computes finds elements in which a phase boundary occurs


% Reactions are when 2 points have a different assemblage number
el_ext          =   elements(:,[1:4, 1]);         % circular extended elements
Reactions       =   NumAssemblage(el_ext(:,2:end)) - NumAssemblage(el_ext(:,1:end-1));


% Grid refinement where the reactions happen
refine_Elements = find(abs(Reactions)>0,2);



% 3) Compute the locations of all reactions
% find all reactions on PhaseDiagram:
% TotalReactions  =   unique((Reactions(:)));
% TotalReactions  =   TotalReactions(TotalReactions~=0);

for iReact = 1:length(Reactions_number_CELL)
    Reac = Reactions_number_CELL(iReact);
    ind  = find(Reactions_number==Reac);
    id   = elements4(ind);      % points
    Loc  = TP_vec(id,:);
    
    
    % NOTE: there is somewhere a bug, so I do not compute mid-points for
    % now
    
    % next, we re-order the points to make a continuous curve
    React_output{iReact}.TP     = Loc;
    React_output{iReact}.Loc    = Reac;
end





% Reactions=[];

