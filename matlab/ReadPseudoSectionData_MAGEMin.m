function [PhaseData, Status] = ReadPseudoSectionData_MAGEMin(varargin)
% This reads the output of a MAGEMin simulation (a file) & stores it in a structure

if nargin==0
    newPoints           =   1;
    PhaseData           =   [];
    numPoints           =   1;
    MinPhaseFraction    =   0;
    labelstyle          =   0;
    db                  =   'ig';
elseif nargin==1
    newPoints           =   varargin{1};
    PhaseData           =   [];
    MinPhaseFraction    =   0;
elseif nargin==2
    newPoints           =	varargin{1};
    PhaseData           =	varargin{2};
    MinPhaseFraction    =   0;
elseif nargin==4
    newPoints           =	varargin{1};
    PhaseData           =	varargin{2};
    MinPhaseFraction    =   varargin{3};
elseif nargin==5
    newPoints           =	varargin{1};
    PhaseData           =	varargin{2};
    MinPhaseFraction    =   varargin{3};
    labelstyle          =   varargin{4};
    db                  =   varargin{5};
else
    error('wrong number of input parameters')
    
end

lbl_exist = 0;
if labelstyle == 0 %then this is default TC style
    % this part load the labelling option files
    if strcmp(db,'ig') | strcmp(db,'igd')
        in = readcell('./matlab/labels/ig.txt','delimiter','');
        lbl_exist = 1;
    end
    if strcmp(db,'alk')
        in = readcell('./matlab/labels/alk.txt','delimiter','');
        lbl_exist = 1;
    end

    if lbl_exist == 1
        for i=1:length(in)
            ph{i} = split(in(i),';');
        end

        j = 1;
        for i=2:size(ph,2)
            if str2double(ph{i}(2)) > 1
                ph_new{j}  = ph{i};
                ph_name(j) = string(ph{i}(1));
                j = j + 1;
            end
        end
    else
        disp(' ');
        disp(' LABEL ERROR: No file is available for this database, if the file has been created (in ./matlab/labels) you need to add a loading line in the matlab script in ./matlab/ReadPseudoSectionData_MAGEMin.m around line 37');
        disp(' -> using Only Phase labeling mode...');
        disp(' ');
    end

end

% open file
fid     = fopen('./output/_pseudosection_output.txt');
fgetl(fid);     % skip comment line

for iPoint=1:length(newPoints)
    % Read line with P/T and Gamma
    line    = fgetl(fid);
    A       = sscanf(line,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    
    % Retrieve info from first (numeric) line
    numPoint= A(1); % the number of the calculation point (in parallel, this may get mixed up)
    STATUS  = A(2);
    P       = A(3);
    T       = A(4);
    Gibbs   = A(5);
    br_norm = A(6);
    Gamma   = A(7:17);
    Vp      = A(18);
    Vs      = A(19);
    entropy = A(20);

    % Read stable assemblage
    StableSolutions =   [];
    StableFractions =   [];
    Density         =   [];
    CompositionalVar=   [];
    EMFractions     =   [];
    EMlist		    =   [];  	%list of end-members name
    SScomp_wt       =   [];
    PHfractions_wt  =   [];
    SSoxide         =   [];  
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
        EM_Frac               = [];
        EM_list               = [];
        SS_comp_wt            = [];
        SS_oxide              = [];

        if str2num(out{4}) > 0
            % We have a solution model; read in the compositional variables
            % and their proportions
            
            n_xeos                = data(4);
            for j=1:n_xeos
                CompVar(j) =  str2num(out{j+4});
            end

            % load endmember fraction and name
            for j=1:n_xeos+1
                EM_Frac(j) =  str2num(out{j*2+4+n_xeos});
            end
            for j=1:n_xeos+1
                EM_list{j} =  out{j*2+4+n_xeos - 1};
            end

            % load solution phase composition and oxide names
            len_ox = str2num(out{n_xeos+(n_xeos+1)*2+5});
            start  = n_xeos+(n_xeos+1)*2 + 5;
            k = 1;
            for j=2:2:len_ox*2
                SS_comp_wt(k)  =  str2num(out{start+j});
                SS_oxide{k} =  out{start+j - 1};
                k = k + 1;
            end

            PH_fractions_wt = str2num(out{start+len_ox*2 + 1});

        else
            len_ox  = str2num(out{5});
            k = 1;
            for j=2:2:len_ox*2
                SS_comp_wt(k)  =  str2num(out{5+j});
                SS_oxide{k} =  out{5+j - 1};
                k = k + 1;
            end

            id      = 5 + len_ox*2 + 1;
            PH_fractions_wt = str2num(out{id});
            n_xeos  = 0;
        end

        CompositionalVar{i} =  CompVar; % the compositional variables for this file
        EMFractions{i} 		=  EM_Frac; % the compositional variables for this file
        EMlist{i}  			=  EM_list;
        SScomp_wt{i}        =  SS_comp_wt;
        PHfractions_wt{i}   =  PH_fractions_wt;
        SSoxide{i}          =  SS_oxide;
        line	            =  fgetl(fid);
        i                   =  i+1;
    end
    
    % routine to change the name of the phase according to secondary parameters
    if labelstyle == 0 & lbl_exist == 1     %default thermocalc style of labeling
        for i=1:size(StableSolutions,1)
            ps_ph_name = string(StableSolutions{i});
            ps_ph_name = split(ps_ph_name,'_');
            % disp(ps_ph_name)
            if length(ps_ph_name) > 1
                ps_ph_name = ps_ph_name(1);
            end

            if any(strcmp(ph_name,ps_ph_name))
                id = find(strcmp(ph_name,ps_ph_name));

                x = CompositionalVar{i};
                islabeled = 0;

                maxv = round(str2double(ph_new{id}{2}));
                k = 0;
                while islabeled == 0 & k < maxv-1
                    val = eval(ph_new{id}{3+k});
                    if val > 0.0
                        lbl = ph_new{id}{3+maxv+k};
                        islabeled = 1;
                    end
                    k = k + 1;

                end
                if islabeled == 0
                    lbl = ph_new{id}{3+maxv*2-1};
                end 
                StableSolutions{i} = lbl;

            end
        end
    end

    if labelstyle == 2     %label using dominant end-member fraction
        for i=1:size(StableSolutions,1)
            % first get the phase name without the suffix for solvii
            ps_ph_name = string(StableSolutions{i});
            ps_ph_name = split(ps_ph_name,'_');

            if length(ps_ph_name) > 1
                ps_ph_name = ps_ph_name(1);
            end

            p = EMFractions{i};

            if length(p) > 1                                % is this a solution phase
                id      = find(p == max(p));                % find the id of the dominant endmember
                em_name = string(EMlist{i}{id});            % get endmember name
                ph_lbl  = strcat(ps_ph_name,':',em_name);   % create label
                StableSolutions{i} = ph_lbl{1};                % attribute label
            end

        end
    end

    % Depending on MinPhaseFraction, we may not take minor phases into account
    % while plotting the diagram, even when they are taken into account in
    % computing the Gibbs energy and Gamma, etc.
    % We also store the full info in the structure below, such that the GUI
    % can later adapt the value of MinPhaseFraction again, without need to
    % recompute all again
    FullInfo.StableSolutions    = StableSolutions;
    FullInfo.StableFractions    = StableFractions;
    FullInfo.EMlist    			= EMlist;
    FullInfo.CompositionalVar	= CompositionalVar;
    FullInfo.EMFractions        = EMFractions;
    FullInfo.Density            = Density;
    FullInfo.SScomp_wt    		= SScomp_wt;
    FullInfo.SSoxide    	    = SSoxide;
    FullInfo.PHfractions_wt     = PHfractions_wt;


    ind                         =   find(StableFractions>MinPhaseFraction);
    StableFractions             =   StableFractions(ind);
    EMlist             			=   EMlist(ind);
    StableSolutions             =   StableSolutions(ind);
    CompositionalVar            =   CompositionalVar(ind);
    EMFractions                 =   EMFractions(ind);
    Density                     =   Density(ind);
    SScomp_wt             		=   SScomp_wt(ind);
    SSoxide             		=   SSoxide(ind);
    PHfractions_wt              =   PHfractions_wt(ind);

    liq                         =   0;
    Density_liq                 =   0;

    % Extract melt fraction
    ind_liq = find(ismember(StableSolutions,'liq'));
    if ~isempty(ind_liq) > 0
        liq = StableFractions(ind_liq(1));
        solvus = 0;
    end

    ind_liq2 = find(ismember(StableSolutions,'liq_2'));
    if ~isempty(ind_liq2) > 0
        liq = sum(StableFractions(ind_liq2));
        solvus = 1;
    end
    


    % Compute average Density of full assemblage
    Density_total           =   sum(StableFractions.*Density);
    
    % Compute Density of liq
    if liq>0

        if solvus == 0
            Density_liq         =   Density(ind_liq(1));
            ind_sol             =   find(~ismember(StableSolutions,'liq'));
        elseif solvus == 1
            Density_liq         =   mean(Density(ind_liq2));
            ind_sol             =   find(~ismember(StableSolutions,'liq_2'));
        end

        if length(ind_sol)>0
            Density_sol     	=   sum(StableFractions(ind_sol).*Density(ind_sol))/sum(StableFractions(ind_sol));
        else
            Density_sol         =   0;
        end
        
    else
        Density_liq         =   2000;
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
    EMlist                                          =   EMlist(iS);
    SScomp_wt                                       =   SScomp_wt(iS);
    SSoxide                                         =   SSoxide(iS);
    PHfractions_wt                                  =   PHfractions_wt(iS);

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
    PhaseData{newPoints(numPoint)}.Vp               =   Vp;
    PhaseData{newPoints(numPoint)}.Vs               =   Vs;
    PhaseData{newPoints(numPoint)}.entropy          =   entropy;

    PhaseData{newPoints(numPoint)}.StableSolutions  =   StableSolutions;
    PhaseData{newPoints(numPoint)}.StableFractions  =   StableFractions;
    PhaseData{newPoints(numPoint)}.EMlist  			=   EMlist;
    PhaseData{newPoints(numPoint)}.Density          =   Density;
    PhaseData{newPoints(numPoint)}.CompositionalVar =   CompositionalVar;
    PhaseData{newPoints(numPoint)}.EMFractions      =   EMFractions;
    PhaseData{newPoints(numPoint)}.liq              =   liq;
    PhaseData{newPoints(numPoint)}.numStablePhases  =   length(StableFractions);
    
    PhaseData{newPoints(numPoint)}.SScomp_wt        =   SScomp_wt;
    PhaseData{newPoints(numPoint)}.PHfractions_wt   =   PHfractions_wt;
    PhaseData{newPoints(numPoint)}.SSoxide          =   SSoxide;
    % Store info of all phases, included the ones that are discarded
    % because of small mass fraction:
    PhaseData{newPoints(numPoint)}.FullInfo         =   FullInfo;
    Status(numPoint)                                =   STATUS;
end

fclose(fid);




