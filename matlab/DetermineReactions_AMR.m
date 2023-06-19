function [refine_Elements] = DetermineReactions_AMR(PseudoSectionData)
% This computes reactions from the phase data

PhaseData   =   PseudoSectionData.PhaseData;

for i=1:length(PhaseData)
    array(i).StableSolutions = PhaseData{i}.StableSolutions;
    % disp(array(i).StableSolutions);
    % disp('\n');
end
    

% 1) Give every solid solution and pure phase a number
% 1a) Create a list with all stable phases/solutions in the diagram
StableSol = PhaseData{1}.StableSolutions;
for i=1:length(PhaseData)
    sol = unique(PhaseData{i}.StableSolutions);
    for j=1:length(sol)
        if ~(any(strcmp(StableSol,sol(j))))
            StableSol(end+1)  = sol(j);
        end
    end
    
end
StableSol = unique(StableSol);

% 1b) Give them all a number
% StablePhaseNumber = [1:length(StableSol)];  % NOTE: we may need to assign random numbers
rng(0,'twister');   % reproducible random numbers
StablePhaseNumber = randi([0 1e4],1,length(StableSol));


% 1c)
for i=1:length(array)
    NumAssemblage = 0;
    sol = array(i).StableSolutions;
    sol_unique = unique(sol);
    for j=1:length(sol_unique)
        id              =   find(strcmp(StableSol,sol_unique{j}));
        NumAssemblage   =   NumAssemblage + StablePhaseNumber(id);
        
    end
    NumAssemblage_vec(i) = NumAssemblage;
end

% Generate a new refined mesh
XY_vec      =   PseudoSectionData.XY_vec;
elements    =   PseudoSectionData.elements;
newPoints   =   [];
PrPts0      =   [];

% Reactions are when 2 points have a different assemblage number
el_ext          =   elements(:,[1:4, 1]);         % circular extended elements
Reactions       =   NumAssemblage_vec(el_ext(:,2:end)) - NumAssemblage_vec(el_ext(:,1:end-1));

% Grid refinement where the reactions happen
refine_Elements = find(any(abs(Reactions)>0,2));


