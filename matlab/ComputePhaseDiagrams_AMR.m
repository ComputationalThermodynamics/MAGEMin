    function [PseudoSectionData, CancelComputation] = ComputePhaseDiagrams_AMR(PseudoSectionData, DisplayPlots)
% This employs an adaptive mesh refinement approach to compute  
% pseudosections, using the MAGEMin code
%
%

% Specify the parameters to compute a pseudosection
VerboseLevel                    =   PseudoSectionData.Computation.Verbose;                %   How much info will be shown during computations?

Chemistry                    	=   PseudoSectionData.Chemistry;               % predefined chemical composition
PseudoSectionData0              =   PseudoSectionData;

if isfield(PseudoSectionData,'fig')
    fig = PseudoSectionData.fig;
else
    fig = uifigure;
end


% Generate initial (coarse) mesh
P_dat                           =   PseudoSectionData.CoarseGrid.P;
T_dat                           =   PseudoSectionData.CoarseGrid.T;

P_vec1D                         =   P_dat(1):P_dat(3):P_dat(2);
T_vec1D                         =   T_dat(1):T_dat(3):T_dat(2);
numPoints                       =   0;
[TP_vec, elements, irregular]   =   GenerateRegularMesh_AMR(T_vec1D,P_vec1D);
newPoints                       =   numPoints+1:size(TP_vec,1);             % newly added points

PhaseData                       =   [];
FailedSimulations               =   [];

% Start refinement loop
for iter=1:PseudoSectionData.Computation.RefinementLevels
    
    % refine existing grid
    if iter>1
        
        % Determine which cells to refine
        [refine_Elements]                   =   DetermineReactions_AMR(PseudoSectionData);
        
        % refine mesh
        numPoints                           =   length(TP_vec);
        elements_old                        =   elements;
        [TP_vec,elements,irregular, PrPts]  =   QrefineR(TP_vec,elements,irregular,refine_Elements);
        newPoints                           =   numPoints+1:size(TP_vec,1); % newly added points
        
        % Update Gamma on the new Points
        if PseudoSectionData.Computation.UseGammaEstimation
            PhaseData = Update_Gamma_onPoints(PhaseData,TP_vec, elements_old, newPoints,PrPts,PseudoSectionData.Computation);
        end
        
        % Update x-EOS on new points from coarse mesh:
        if PseudoSectionData.Computation.Use_xEOS
            PhaseData = Update_xEOS_onPoints(PhaseData,TP_vec, newPoints,PrPts,PseudoSectionData.Computation);
        end
        
        % Update EM Fractions on new points from coarse mesh:
        if PseudoSectionData.Computation.Use_EMFrac
            PhaseData = Update_EMFrac_onPoints(PhaseData,TP_vec, newPoints,PrPts,PseudoSectionData.Computation);
        end
        
        
        if DisplayPlots
            figure(2), clf
            patch('Faces', elements, 'Vertices', TP_vec, 'Facecolor','none')
            hold on
            plot(TP_vec(newPoints,1),TP_vec(newPoints,2),'ro')
            title(['requires computation of ',num2str(length(newPoints)),' new points (red)'])
            drawnow
        end
        
    else
        
        % Perhaps we can already make an estimation of Gamma
        if PseudoSectionData.Computation.UseGammaEstimation
            switch PseudoSectionData.Computation.Gamma_Method
                case {'Inherit from coarse mesh','AMR interpolation','Average surrounding points'}
                    % do nothing, as previous Gamma is not yet known
                    
                otherwise
                    elements_old                        =   elements;
                    PrPts=[];
                    PhaseData = Update_Gamma_onPoints(PhaseData,TP_vec, elements_old, newPoints,PrPts,PseudoSectionData.Computation);
                    
            end
        end
        
        % Perhaps we can also already make an estimation of x-eos
        if PseudoSectionData.Computation.Use_xEOS
            switch PseudoSectionData.Computation.EOS_Method
                case 'Neural Network'
                    PrPts     = [];
                    PhaseData = Update_xEOS_onPoints(PhaseData,TP_vec, newPoints,PrPts,PseudoSectionData.Computation);
                    
                otherwise
                    % nope, no info about xEOS yet
            end
                
        end
        
        % Perhaps we can also already make an estimation of x-eos
        if PseudoSectionData.Computation.Use_EMFrac
            switch PseudoSectionData.Computation.EMFrac_Method
                case 'Neural Network'
                    PrPts     = [];
                    PhaseData = Update_EMFrac_onPoints(PhaseData,TP_vec, newPoints,PrPts,PseudoSectionData.Computation);
                    
                otherwise
                    % nope, no info about xEOS yet
            end
                
        end
        
        
    end
    
    % Perform a MAGEMin simulation for all new points
    dlg   = uiprogressdlg(fig,'Title',['Performing computations; refinement level ',num2str(iter),' of ',num2str(PseudoSectionData.Computation.RefinementLevels)],'Cancelable','on');
    ComputeAllPoints    =   PseudoSectionData.Computation.AllPoints;
    UseGammaEstimation  =   PseudoSectionData.Computation.UseGammaEstimation; 
    Computation         =   PseudoSectionData.Computation;
    
    [PhaseData, TP_vec, FailedSimulations, CancelComputation] = PerformMAGEMin_Simulation(PhaseData, newPoints, TP_vec, VerboseLevel, Chemistry, dlg, ComputeAllPoints, UseGammaEstimation, Computation);
    
    if CancelComputation
        break
    end
    
%     array       =   cell2mat(PhaseData);
    
    % liquid fraction at every point
    liq         = zeros(length(PhaseData),1);
    numPhase    = zeros(length(PhaseData),1);
    for i=1:length(PhaseData)
       
        
        liq(i)      =   PhaseData{i}.liq(1);
        numPhase(i) =   PhaseData{i}.numStablePhases;
    end
    
    
    if DisplayPlots
        % Plot the results
        % color of mesh:
        col         =   liq;
        
        figure(1), clf
        h = patch('Faces', elements, 'Vertices', TP_vec,'FaceVertexCData',col(:),'FaceColor','interp');
        xlabel('Temperature [C]','Fontsize',15)
        ylabel('Pressure [kbar]','Fontsize',15)
        colorbar, title('melt fraction')
        
        
        col         =   numPhase;
        figure(3), clf
        h = patch('Faces', elements, 'Vertices', TP_vec,'FaceVertexCData',col(:),'FaceColor','interp');
        xlabel('Temperature [C]','Fontsize',15)
        ylabel('Pressure [kbar]','Fontsize',15)
        colorbar, title('# of stable phases')
    end
    
    
    % Store data in structure, which we will later use in the GUI
    PseudoSectionData.TP_vec            =   TP_vec;
    PseudoSectionData.elements          =   elements;
    PseudoSectionData.PhaseData         =   PhaseData;
    PseudoSectionData.FailedSimulations =   FailedSimulations;
    
    pause(1)
end


if CancelComputation % cancel button was pushed; set back data
    PseudoSectionData = PseudoSectionData0;
end


