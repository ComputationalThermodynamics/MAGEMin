function PseudoSectionData = ReprocessPointsMinPhaseFraction(PseudoSectionData,MinPhaseFraction)
% This goes through all points and adapts them based on a newly selected
% minimum phase fraction. In this manner, we don't have to redo the phase
% diagram calculations



PhaseData = PseudoSectionData.PhaseData;

for iPoint = 1:length(PhaseData)
    
    PD =  PhaseData{iPoint};
    
    FullInfo                =   PD.FullInfo;
    StableSolutions         =   FullInfo.StableSolutions;
    StableFractions         =   FullInfo.StableFractions;
    CompositionalVar        =   FullInfo.CompositionalVar;
    Density                 =   FullInfo.Density;
    
    % Reprocess:
    ind                  	=   find(StableFractions>MinPhaseFraction);
    StableFractions      	=   StableFractions(ind);
    StableSolutions     	=   StableSolutions(ind);
    CompositionalVar      	=   CompositionalVar(ind);
    Density              	=   Density(ind);
    
     % Extract melt fraction
    ind_liq = find(ismember(StableSolutions,'liq'));
    if ~isempty(ind_liq)>0
        liq = StableFractions(ind_liq);
    else
        liq = 0;
    end
    
    % Compute average Density of full assemblage
    Density_total           =   sum(StableFractions.*Density);    
    
    % Compute Density of liq
    if liq>0
        Density_liq         =   Density(ind_liq);
        
        ind_sol             =   find(~ismember(StableSolutions,'liq'));
        if length(ind_sol)>0
            Density_sol     	=   sum(StableFractions(ind_sol).*Density(ind_sol))/sum(StableFractions(ind_sol));
        else
            Density_sol         =   0;
        end
        
    else
        Density_liq         =   0;
        Density_sol         =   Density_total;
    end
    
    
     % Save output to structure
    PD.Density_total    =   Density_total;
    PD.Density_sol      =   Density_sol;
    PD.Density_liq      =   Density_liq;
    PD.StableSolutions  =   StableSolutions;
    PD.StableFractions  =   StableFractions;
    PD.Density          =   Density;
    PD.CompositionalVar =   CompositionalVar;
    PD.numStablePhases  =   length(StableFractions);
    
    PhaseData{iPoint}   =   PD; % set back
end

PseudoSectionData.PhaseData = PhaseData;

            