function [PhaseData, XY_vec, FailedSimulations, CancelComputation]  = 	PerformMAGEMin_Simulation(PseudoSectionData,PhaseData, newPoints, XY_vec, VerboseLevel, Chemistry, dlg, ComputeAllPoints, UseGammaEstimation, Computation);
% This performs a MAGEMin simulation for a bunch of points

DiagramType = Computation.DiagramType;
Mode        = {'AllPoints','SinglePoint'};
if ComputeAllPoints
    Mode = Mode{1};
else
    Mode = Mode{2};
end
sys_in              =   Chemistry.sys_in;               % predefined chemical composition
Test                =   Chemistry.Predefined;
if isnan(Test)
    OxProp         =   Chemistry.OxProp;    % we do not employ a predefined test, but specify mol proportions insteadOxProp
else
    OxProp         =   [];
end

Use_xEOS            =   Computation.Use_xEOS;
Use_EMFrac          =   Computation.Use_EMFrac;

% Delete GRID files before starting a new run
cur_dir = pwd;
cd('output')
system('rm -rf _pseudosection_output.*.*');
cd(cur_dir);

CancelComputation   =   false;
PhaseData0          =   PhaseData;
switch Mode
    case 'AllPoints'
        
        % Write all points to be processed to a single file
        n_points = Write_MAGEMin_InputFile(newPoints, XY_vec, PhaseData, UseGammaEstimation, Use_xEOS, Use_EMFrac, OxProp, DiagramType,PseudoSectionData);
        if ~isempty(PhaseData)
            for i=1:length(newPoints)
                
                if length(PhaseData)>newPoints(i)
                    if isfield(PhaseData{newPoints(i)},'Gamma')
                        Gamma = PhaseData{newPoints(i)}.Gamma;
                    else
                        Gamma = [];
                    end
                    PhaseData{newPoints(i)}.StartingValue.Gamma = Gamma;
                    
                end
                
                % store this data for debugging later (so we can reconstruct how
                % the point was calculated)
                PhaseData{newPoints(i)}.StartingValue.Use_xEOS      = Use_xEOS;
                PhaseData{newPoints(i)}.StartingValue.Use_Gamma     = UseGammaEstimation;
                PhaseData{newPoints(i)}.StartingValue.Computation   = Computation;
                
            end
        end
        
        dlg.Message         = [num2str(length(newPoints)),' points; computing them all @ once. '];
        dlg.Indeterminate   = 'on';
        
        % Perform the MAGEMin simulation for all the points
%         % Run simulation
%         NumRanks    =   Computation.NumRanks;
%         MPI_dir     =   Computation.MPI_dir;
%         tic
%         if NumRanks>1
%             command = [MPI_dir,'/mpiexec -n ',num2str(NumRanks),' '];
%         else
%             command = '';
%         end
%         command = [command, './MAGEMin --Verb=',num2str(VerboseLevel),' --File=MAGEMin_input.dat --n_points=',num2str(n_points)];
%         
%         if ~isnan(Test)
%             % employ a prededined test
%             command = [command, ' --test=',num2str(Test)];
%         else
%             % employ specified chemistry
%             command = [command, ' --Bulk='];
%             for iTable=1:size(OxProp,1)
%                 command = [command,num2str(table2array(OxProp(iTable,2))),','];
%             end
%         end
%         system('killall MAGEMin 2>&1')
%         disp(command)
%         system(command);
%         
        
        tic
        Execute_MAGEMin(Computation, VerboseLevel, n_points, Test, OxProp, sys_in)
        ForwardSimulation_Time = toc
        
        if dlg.CancelRequested
            CancelComputation = true;
        end
        
        % Read data of all points
        tic
        if Computation.MatlabOut
            [PhaseData, Status] = ReadEquilibriumPathData_MAGEMin(newPoints, PhaseData);
        else
            [PhaseData, Status] = ReadPseudoSectionData_MAGEMin(newPoints, PhaseData, Computation.MinPhaseFraction,Computation.labelstyle,Computation.db);
        end
        ReadData_Time = toc
        
    case 'SinglePoint'
        % we do computations point-by-point
        
        ForwardSimulation_Time = 0;
        
        for iPoint=1:length(newPoints)
            disp(['SinglePoint;  iPoint=',num2str(iPoint)])
            id  = newPoints(iPoint);
            
            T   = XY_vec(id,1);
            P   = XY_vec(id,2);
            
            tic;
            
            % Run simulation
            command = ['./MAGEMin --Verb=',num2str(VerboseLevel)];
            if ~isnan(Test)
                % employ a prededined test
                command = [command, ' --test=',num2str(Test)];
            else
                % employ specified chemistry
                command = [command, ' --Bulk='];
                for iTable=1:size(OxProp,1)
                    command = [command,num2str(table2array(OxProp(iTable,2))),','];
                end
            end
            % store this data for debugging later (so we can reconstruct how
            % the point was calculated)
            StartingValue.Use_xEOS      = Use_xEOS;
            StartingValue.Use_Gamma     = UseGammaEstimation;
            StartingValue.Computation   = Computation;
            
            % if we use an xEOS, we need to write a file and use that
            if iPoint>1
                id = newPoints(iPoint-1:iPoint);
            end
            
            n_points = Write_MAGEMin_InputFile(id, XY_vec, PhaseData, UseGammaEstimation, Use_xEOS, Use_EMFrac, OxProp, DiagramType, PseudoSectionData);
           
            Computation.NumRanks=1; % makes no sense to do a single calculation on >1 core 
            Execute_MAGEMin(Computation, VerboseLevel, n_points, Test, OxProp, sys_in)
            
            ForwardSimulation_Time = ForwardSimulation_Time + toc;
            disp(['Percentage finished ',num2str(iPoint/length(newPoints)*100),'% ; total forward computational time = ',num2str(ForwardSimulation_Time),'s'])
            disp(' ')
            
            % Read data for this point
            if iPoint>1
                PhaseData_in1   = ReadPseudoSectionData_MAGEMin(1:2,[],Computation.MinPhaseFraction,Computation.labelstyle,Computation.db);
                PhaseData_in{1} = PhaseData_in1{2};
                id = id(2);
            else
                PhaseData_in   = ReadPseudoSectionData_MAGEMin(1,[],Computation.MinPhaseFraction,Computation.labelstyle,Computation.db);
            end

            PhaseData{id}  = PhaseData_in{:};
            PhaseData{id}.StartingValue = StartingValue; % save starting values
            Status(iPoint) = PhaseData_in{:}.Status;
            
            % Update progress bar
            if mod(iPoint,10)==0
                dlg.Value   = iPoint/length(newPoints);
                dlg.Message = [num2str(length(newPoints)),' points; ', num2str(round(iPoint/length(newPoints)*1000)/10),'% done.. '];
                pause(1/1000)
            end
            
            if dlg.CancelRequested
                CancelComputation = true;
                break
            end
            
        end
        
        
end
FailedSimulations=[];


FailedSimulations   =   [];
indFailed           =   find( Status==3 );      % failed to converge
indSuccess          =   find( Status <3 );      % did work

if CancelComputation
    indFailed = [];
end



%--------------------------------------------------------------------------
% Write input file for MAGEMin
function n_points = Write_MAGEMin_InputFile(newPoints, XY_vec, PhaseData, Use_Gamma, Use_xEOS, Use_EMFrac, OxProp, DiagramType, PseudoSectionData)

fid         = fopen('MAGEMin_input.dat','w');
n_points    = length(newPoints);
for i=1:n_points
    % Read line with P/T and Gamma

    % in case PX or TX diagram is selected the bulk is set
    nox    = size(PseudoSectionData.Chemistry.OxProp(:,2),1);
    Bulk   =   zeros(1,nox);

    if DiagramType == 'PT'
        id      =   newPoints(i);
        T       =   XY_vec(id,1);
        P       =   XY_vec(id,2);
    elseif DiagramType == 'PX'
        id      =   newPoints(i)
        P       =   XY_vec(id,2)
        T       =   (PseudoSectionData.FixedTemperature)
        Bulk    =   table2array(PseudoSectionData.Chemistry.OxProp(:,2))*(1.-XY_vec(id,1)) + table2array(PseudoSectionData.Chemistry.OxProp(:,3))*(XY_vec(id,1));
    elseif DiagramType == 'TX'
        id      =   newPoints(i);
        T       =   XY_vec(id,2);
        P       =   PseudoSectionData.FixedPressure;
        Bulk    =   table2array(PseudoSectionData.Chemistry.OxProp(:,2))*(1.-XY_vec(id,1)) + table2array(PseudoSectionData.Chemistry.OxProp(:,3))*(XY_vec(id,1));
    end
    

    % if we do not use the previous Gamma
    Gamma   =   zeros(1,11);
    if Use_Gamma
        if length(PhaseData)>=id
            % check if a previous Gamma estimatione exist
            if isfield(PhaseData{id},'Gamma')
                Gamma = PhaseData{id}.Gamma;
            end
        end
    end
    
    
    % for now assume that we don't use a previous guess of phases & initial guess of compostional variables.
    % note: this could be changed in the future!
    if (~Use_xEOS & ~Use_EMFrac)
        n_phases = 0;
    else
        if length(PhaseData)>=id
            if isfield(PhaseData{id},'StartingValues_xEOS')
                % We only do this if we have specified StartingValues for a specific point
                % They may not exists in the first loop, for example.
                
                StartingValues_xEOS     =   PhaseData{id}.StartingValues_xEOS;
                SolutionModels_Comp     =   StartingValues_xEOS.SolutionModels;
                InitialGuess_Num_Comp   =   zeros(size(SolutionModels_Comp));            % Change that if we want >1
                CompositionalVar        =   StartingValues_xEOS.CompositionalVar;
                n_phases_Comp           =   length(SolutionModels_Comp);                 % solution models
            else
                n_phases_Comp           =   0;
                SolutionModels_Comp     =   [];
            end
            
            if isfield(PhaseData{id},'StartingValues_EMFraction')
                % We only do this if we have specified StartingValues for a specific point
                % They may not exists in the first loop, for example.
                
                StartingValues_EMFraction   =   PhaseData{id}.StartingValues_EMFraction;
                SolutionModels_EM           =   StartingValues_EMFraction.SolutionModels;
                InitialGuess_Num_EM         =   zeros(size(SolutionModels_EM));            % Change that if we want >1
                EM_Fractions                =   StartingValues_EMFraction.EMFractions;
                n_phases_EM                 =   length(SolutionModels_EM);                 % solution models
            else
                n_phases_EM                 =   0;
                SolutionModels_EM           =   [];
            end
            
            n_phases = max([n_phases_EM n_phases_Comp]);
            
            if n_phases>0
                if n_phases_Comp>=n_phases_EM
                    SolutionModels   = SolutionModels_Comp;
                    InitialGuess_Num = InitialGuess_Num_Comp;
                else
                    SolutionModels   = SolutionModels_EM;
                    InitialGuess_Num = InitialGuess_Num_EM;
                end
            end
        else
            n_phases = 0;
        end
    end

    if nox == 11
        fprintf(fid,'%i %f %f %f %f %f %f %f %f %f %f %f %f %f\n',n_phases, P,T,Bulk);
    elseif nox == 10
        fprintf(fid,'%i %f %f %f %f %f %f %f %f %f %f %f %f\n',n_phases, P,T,Bulk); 
    elseif nox == 7
        fprintf(fid,'%i %f %f %f %f %f %f %f %f %f\n',n_phases, P,T,Bulk); 
    elseif nox == 8
        fprintf(fid,'%i %f %f %f %f %f %f %f %f %f %f\n',n_phases, P,T,Bulk); 
    end
    if (Use_xEOS | Use_EMFrac) & n_phases>0
        % Add compositional variable guesses
        for i=1:n_phases
            SS = SolutionModels{i};
            
            fprintf(fid,'%s ',SS  );
%             fprintf(fid,'%i ',InitialGuess_Num(i));     % we can expand this later to take
           
            CompVar_vec	= zeros(1,11);
            EM_Frac_vec = zeros(1,12);
            if Use_xEOS
                if length(SolutionModels_Comp)>0
                    id                  =   find(contains(SolutionModels_Comp,SS));
                    if length(id)>0
                        CompVar             =   CompositionalVar{id};
                        n                   =   length(CompVar);
                        CompVar_vec(1:n)    =   CompVar;
                    end
                end
            end
            
            if Use_EMFrac
                if length(SolutionModels_EM)>0
                    id                  =   find(contains(SolutionModels_EM,SS));
                    if length(id)>0
                        EM_Frac             =   EM_Fractions{id};
                        n                   =   length(EM_Frac);
                        EM_Frac_vec(1:n)    =   EM_Frac;
                    end
                end
            end
            
            % Write to file. Note that we have a fxed # of values we write
            for iVar=1:11
                fprintf(fid,'%f ',CompVar_vec(iVar));
            end
            
            for iVar=1:12
                fprintf(fid,'%f ',EM_Frac_vec(iVar));
            end
            
            
            fprintf(fid,' \n');
        end
        
    end
    
end
fclose(fid);
