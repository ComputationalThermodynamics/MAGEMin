function [PhaseData, Status] = ReadData_MAGEMin(varargin)
% This reads the output of a MAGEMin simulation (a file) & stores it in a structure

if nargin==0
    newPoints           =   1;
    PhaseData           =   [];
    numPoints           =   1;
    MinPhaseFraction    =   0;
elseif nargin==1
    newPoints           =   varargin{1};
    PhaseData           =   [];
    MinPhaseFraction    =   0;
elseif nargin==2
    newPoints           =	varargin{1};
    PhaseData           =	varargin{2};
    MinPhaseFraction    =   0;
elseif nargin==3
    newPoints           =	varargin{1};
    PhaseData           =	varargin{2};
    MinPhaseFraction    =   varargin{3};
else
    error('wrong number of input parameters')
    
end


% open file
fid     = fopen('./output/_pseudosection_output.txt');
fgetl(fid);     % skip comment line

for iPoint=1:length(newPoints)
    % Read line with P/T and Gamma
    line = fgetl(fid);
    A       = sscanf(line,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    
    % Retrieve info from first (numeric) line
    numPoint= A(1); % the number of the calculation point (in parallel, this may get mixed up)
    STATUS  = A(2);
    P       = A(3);
    T       = A(4);
    Gibbs   = A(5);
    br_norm = A(6);
    Gamma   = A(7:end);
    
    % Read stable assemblage
    StableSolutions =   [];
    StableFractions =   [];
    Density         =   [];
    CompositionalVar=   [];
    EMFractions     =   [];
    
    i               =   1;
    line            =   fgetl(fid);
    
    
    while ~isempty(line)
        out                   = split(line);
        data                  = str2double(out); 
        StableSolutions{i,1}  = out{1};
        StableFractions(i)    = data(2);
        Density(i)            = data(3);
%         StableFractions(i)    = str2num(out{2});
%         Density(i)            = str2num(out{3});
%         n_xeos                = str2num(out{4});
        
        CompVar               = [];
        EM_Frac                = [];
        if length(out)>4
            % We have a solution model; read in the compositional variables
            % and their proportions
            
            n_xeos                = data(4);
            for j=1:n_xeos
                CompVar(j) =  str2num(out{j+4});
                %CompVar(j) =  data(j+7);
            end
            for j=1:n_xeos+1
                EM_Frac(j) =  str2num(out{j+4+n_xeos});
                %EM_Frac(j) =  data(j+7+n_xeos);
            end
        else
            n_xeos=0;
        end
         
        CompositionalVar{i} = CompVar; % the compositional variables for this file
        EMFractions{i} = EM_Frac; % the compositional variables for this file
        
        line	=   fgetl(fid);
        i       =   i+1;
    end
    
    % Depending on MinPhaseFraction, we may not take minor phases into account
    % while plotting the diagram, even when they are taken into account in
    % computing the Gibbs energy and Gamma, etc.
    % We also store the full info in the structure below, such that the GUI
    % can later adapt the value of MinPhaseFraction again, without need to
    % recompute all again
    FullInfo.StableSolutions    = StableSolutions;
    FullInfo.StableFractions    = StableFractions;
    FullInfo.CompositionalVar	= CompositionalVar;
    FullInfo.EMFractions        = EMFractions;
    FullInfo.Density            = Density;
    
    
    ind                         =   find(StableFractions>MinPhaseFraction);
    StableFractions             =   StableFractions(ind);
    StableSolutions             =   StableSolutions(ind);
    CompositionalVar            =   CompositionalVar(ind);
    EMFractions                 =   EMFractions(ind);
    Density                     =   Density(ind);
    
    
    % Extract melt fraction
    ind_liq = find(ismember(StableSolutions,'liq'));
    if ~isempty(ind_liq)>0
        liq = StableFractions(ind_liq(1));
    else
        liq = 0;
    end
    
    % Compute average Density of full assemblage
    Density_total           =   sum(StableFractions.*Density);
    
    % Compute Density of liq
    if liq>0
        Density_liq         =   Density(ind_liq(1));
        
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
    
    
    % A complication with the MPI-parallel computations is that the results are
    % not in the same order as the original list of points
    %
    
    % Sort according to names (makes things easier to compare later)
    [StableSolutions,iS]                            =   sort(StableSolutions);
    CompositionalVar                                =   CompositionalVar(iS);
    EMFractions                                     =   EMFractions(iS);
    StableFractions                                 =   StableFractions(iS);
    Density                                         =   Density(iS);

    % Save output to structure
    PhaseData{newPoints(numPoint)}.Status           =   STATUS;
    PhaseData{newPoints(numPoint)}.P                =   P;
    PhaseData{newPoints(numPoint)}.T                =   T;
    PhaseData{newPoints(numPoint)}.Gibbs            =   Gibbs;
    PhaseData{newPoints(numPoint)}.Density_total    =   Density_total;
    PhaseData{newPoints(numPoint)}.Density_sol      =   Density_sol;
    PhaseData{newPoints(numPoint)}.Density_liq      =   Density_liq;
    PhaseData{newPoints(numPoint)}.br_norm          =   br_norm;
    PhaseData{newPoints(numPoint)}.Gamma            =   Gamma;
    PhaseData{newPoints(numPoint)}.StableSolutions  =   StableSolutions;
    PhaseData{newPoints(numPoint)}.StableFractions  =   StableFractions;
    PhaseData{newPoints(numPoint)}.Density          =   Density;
    PhaseData{newPoints(numPoint)}.CompositionalVar =   CompositionalVar;
    PhaseData{newPoints(numPoint)}.EMFractions      =   EMFractions;
    PhaseData{newPoints(numPoint)}.liq              =   liq;
    PhaseData{newPoints(numPoint)}.numStablePhases  =   length(StableFractions);
    
    % Store info of all phases, included the ones that are discarded
    % because of small mass fraction:
    PhaseData{newPoints(numPoint)}.FullInfo         =   FullInfo;
    
    Status(numPoint)                                =   STATUS;
end

fclose(fid);




