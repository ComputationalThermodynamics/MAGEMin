function [PseudoSectionData, CancelCalculation] = RecomputePoints_TC(dlg,newPoints,PseudoSectionData)
% Recomputes points using Thermocalc

% Extract data
PhaseData               =   PseudoSectionData.PhaseData;
TP_vec                  =	PseudoSectionData.TP_vec;
Computation             =   PseudoSectionData.Computation;
Chemistry               =   PseudoSectionData.Chemistry;

dlg.Title               =  'Redoing pointwise calculations using ThermoCalc';
dlg.Indeterminate       =   'off';

CancelCalculation       =   false;

% Rerun every point with ThermoCalc. Mainly for debugging and
% development, as it i
for iPoint=1:length(newPoints)
   
    save debug
    
    % Read point
    id              =   newPoints(iPoint);
    Data            =   PhaseData{id};
    
    % Obtain x-EOS and stable phases in vicinity
    Computation.Method          =   'Average surrounding points';
    Computation.num_Neighbours  =   300;
    Computation.EOS_Method      =   Computation.Method;
    PhaseData                   =   Update_xEOS_onPoints(PhaseData,TP_vec, newPoints(iPoint),[],Computation);
    xEOS                        =   PhaseData{id}.StartingValues_xEOS;
    
    % Set solution models & x-EOS as TC input
    Data_TC_in                  =   Data;
    if length(Data_TC_in.StableFractions)<length(xEOS.SolutionModels)
        ind = length(Data_TC_in.StableFractions)+1: length(xEOS.SolutionModels);
        Data_TC_in.StableFractions(ind) = 1e-3;
    end
    Data_TC_in.StableSolutions  =   xEOS.SolutionModels;
    Data_TC_in.CompositionalVar =   xEOS.CompositionalVar;
    Data_TC_in.Chemistry        =   Chemistry;
    
    % Perform TC optimization
    Data_TC_best                =   dogmin_TC(Data_TC_in, './ThermoCalc');
    
    if Data_TC_best.success
        % If success, transfer data
        PhaseData{id}.Gibbs                 =   Data_TC_best.Gibbs;
        PhaseData{id}.StableSolutions       =   Data_TC_best.StableSolutions;
        PhaseData{id}.StableFractions       =   Data_TC_best.StableFractions;
        PhaseData{id}.Gamma                 =   Data_TC_best.Gamma;
        PhaseData{id}.CompositionalVar      =   Data_TC_best.CompositionalVar;
        PhaseData{id}.Density_total         =   Data_TC_best.rho;
        PhaseData{id}.liq                   =   Data_TC_best.liq;
        PhaseData{id}.numStablePhases       =   Data_TC_best.numStablePhases;
        PhaseData{id}.TC_recompute.success  =   true;
        PhaseData{id}.Status                =   0;  % all is fine
    else
        PhaseData{id}.TC_recompute.success  =   false;
    end
    
    % update box
    dlg.Value   = iPoint/length(newPoints);
    dlg.Message = [num2str(length(newPoints)),' points; ', num2str(round(iPoint/length(newPoints)*1000)/10),'% done.. '];
    pause(1/1000)
    
    if dlg.CancelRequested
        CancelCalculation = true;
        break
    end
    disp(['Recalculating points with Thermocalc:  ',num2str(iPoint/length(newPoints)*100),' % finished'])
    
end

% Set back data
if ~CancelCalculation
    PseudoSectionData.PhaseData = PhaseData;
end


