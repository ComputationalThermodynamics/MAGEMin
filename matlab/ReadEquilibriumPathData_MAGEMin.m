function [OUT,Status] = ReadEquilibriumPathData_MAGEMin(varargin)
% This reads the output of a MAGEMin simulation (a file) & stores it in a structure

Status = 1; % failed simulation tracking not added yet

if nargin==0
    newPoints     =   1;
    OUT           =   [];
elseif nargin==1
    newPoints     =   varargin{1};
    OUT           =   [];
elseif nargin==2
    newPoints     =	varargin{1};
    OUT           =	varargin{2};
else
    error('wrong number of input parameters')
end


% open file
fid     = fopen('./output/_matlab_output.txt');
fgetl(fid);     % skip comment line

for iPoint=1:length(newPoints)
    numPoint = newPoints(iPoint);

    % Read line with P/T and Gamma
    fgetl(fid);       % skip line
    A = split(fgetl(fid));   % read and split line

     % retrieve stable phases
    i = 1;
    while isempty(A{1}); A = A(2:end); end
    while ~strcmp(A{i},'{')
        if i>1; if any(strcmp(StablePhases,A{i})); A{i} = [A{i} '2']; end; end
        StablePhases{i} = A{i}; 
        i = i+1;
    end
    i = i+1;

    % extract solution phases
    StableSolutions  = StablePhases(~ismember(StablePhases,{'q' 'crst' 'trd' 'coe' 'stv' 'ky' 'sill' 'and' 'ru' 'sph'}));
    StablePurePhases = StablePhases( ismember(StablePhases,{'q' 'crst' 'trd' 'coe' 'stv' 'ky' 'sill' 'and' 'ru' 'sph'}));

    % retrieve info from first (numeric) line
    P       = str2double(A{i}(1:end-1)); i = i+1;  % read pressure
    T       = str2double(A{i}(1:end-1)); i = i+1;  % read temperature

    % retrieve endmember fractions
    fgetl(fid); fgetl(fid); fgetl(fid);  % skip section header
    for iph = 1:length(StableSolutions)
        A = split(fgetl(fid));    % read and split line
        EMList.(A{2}) = [];
        i = 3;                    % read phase data
        while ~strcmp(A{i},'-') && ~isempty(A{i})
            EMList.(A{2}){i-2} = A{i};
            i = i+1; if i>length(A); break; end
        end
        B = split(fgetl(fid));    % read and split line
        if iph>1; if isfield(EMFractions,B{2}); B{2} = [B{2} '2']; end; end
        EMFractions.(B{2}) = [];
        i = 3;                    % read phase data
        while ~strcmp(B{i},'-') && ~isempty(B{i})
            EMFractions.(B{2})(i-2) = str2double(B{i});
            i = i+1; if i>length(B); break; end
        end
        EMList.(B{2}) = [A{2} EMList.(A{2})];
        EMList = rmfield(EMList,A{2});
    end

    % retrieve site fractions
    fgetl(fid); fgetl(fid); fgetl(fid);  % skip section header
    for iph = 1:length(StableSolutions)
        A = split(fgetl(fid));    % read and split line
        if iph>1; if isfield(SiteFractions,A{2}); A{2} = [A{2} '2']; end; end
        SiteFractions.(A{2}) = [];
        i = 3;                    % read phase data
        while ~strcmp(A{i},'-') && ~isempty(A{i})
            SiteFractions.(A{2})(i-2) = str2double(A{i});
            i = i+1; if i>length(A); break; end
        end
    end

    % retrieve phase compositions
    fgetl(fid); fgetl(fid); fgetl(fid);  % skip section header
    A = split(fgetl(fid));    % read and split line
    OxideList = [];
    i = 2;                    % read phase data
    while ~strcmp(A{i},'-') && ~isempty(A{i}) 
        OxideList{i-1} = A{i};
        i = i+1; if i>length(A); break; end
    end
    for iph = 1:length(StablePhases)+1
        B = split(fgetl(fid));    % read and split line
        if iph>1; if isfield(OxideFractions,B{2}); B{2} = [B{2} '2']; end; end
        OxideFractions.(B{2}) = [];
        i = 3;                    % read phase data
        while ~strcmp(B{i},'-') && ~isempty(B{i})
            OxideFractions.(B{2})(i-2) = str2double(B{i});
            i = i+1; if i>length(B); break; end
        end
    end

    % retrieve phase, system properties
    fgetl(fid); fgetl(fid); fgetl(fid);  % skip section header
    A = split(fgetl(fid));    % read and split line
    PhasePropsList = [];
    i = 3;                    % read phase data
    while ~strcmp(A{i},'-') && ~isempty(A{i}) 
        PhasePropsList{i-2} = A{i};
        i = i+1; if i>length(A); break; end
    end
    for iph = 1:length(StablePhases)
        B = split(fgetl(fid));    % read and split line
        if iph>1; if isfield(PhaseProps,B{2}); B{2} = [B{2} '2']; end; end
        PhaseProps.(B{2}) = [];
        i = 3;                    % read phase data
        while ~strcmp(B{i},'-') && ~isempty(B{i})
            PhaseProps.(B{2})(i-2) = str2double(B{i});
            i = i+1; if i>length(B); break; end
        end
    end
    B = split(fgetl(fid));    % read and split line
    SYSProps = [];
    i = 3;                    % read phase data
    while ~strcmp(B{i},'-') && ~isempty(B{i})
        SYSProps(i-2) = str2double(B{i});
        i = i+1; if i>length(B); break; end
    end
    SYSPropsList = PhasePropsList([2 5 7:end]);

    % retrieve chemical potential of oxides
    fgetl(fid); fgetl(fid); fgetl(fid);  % skip section header
    GammaList = [];
    Gamma     = [];
    i = 1;
    A = split(fgetl(fid));       % read phase data
    while length(A)>1
        GammaList{i} = A{2};
        Gamma(i)     = str2double(A{3});
        A = split(fgetl(fid));
        i = i+1;
    end

    % retrieve delta Gibbs energy, append to phase properties struct
    fgetl(fid); fgetl(fid);  % skip section header
    DeltaGibbsList = [];
    DeltaGibbs     = [];
    i = 1;
    A = split(fgetl(fid));       % read phase data
    while length(A)>1
        DeltaGibbsList{i} = A{2};
        DeltaGibbs(i)     = str2double(A{3});
        A = split(fgetl(fid));
        i = i+1;
    end
    fgetl(fid);

    % normalise phase fractions to unity sum
    sumphs = 0;
    for iph = 1:length(StablePhases)
        sumphs = sumphs+PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'fraction[wt]'));
    end
    for iph = 1:length(StablePhases)
        PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'fraction[wt]')) = PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'fraction[wt]'))./sumphs;
    end

    % Extract phase fractions in wt, vol
    if isfield(PhaseProps,'liq')
        PhaseFractions.liq_wt  = PhaseProps.liq(strcmp(PhasePropsList,'fraction[wt]'));
        PhaseFractions.liq_vol = PhaseProps.liq(strcmp(PhasePropsList,'fraction[wt]')) .* SYSProps(strcmp(SYSPropsList,'Rho[kg/m3]'))./PhaseProps.liq(strcmp(PhasePropsList,'Rho[kg/m3]'));
    else
        PhaseFractions.liq_wt  = 0;
        PhaseFractions.liq_vol = 0;
    end
    PhaseFractions.sol_wt  = 1-PhaseFractions.liq_wt;
    PhaseFractions.sol_vol = 1-PhaseFractions.liq_vol;

    % Extract mixture, melt, solid densities
    if isfield(PhaseProps,'liq')
        Density.liq = PhaseProps.liq(strcmp(PhasePropsList,'Rho[kg/m3]'));
    else
        Density.liq = nan;
        OxideFractions.liq = nan(size(OxideFractions.(StablePhases{1})));
    end
    if PhaseFractions.sol_wt>1e-6
        Density.sol = 0;
        OxideFractions.sol = zeros(size(OxideFractions.(StablePhases{1})));
        for iph = 1:length(StablePhases)
            if ~strcmp(StablePhases{iph},'liq')
                Density.sol = Density.sol + PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'fraction[wt]')) ./ PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'Rho[kg/m3]'));
                OxideFractions.sol = OxideFractions.sol + PhaseProps.(StablePhases{iph})(strcmp(PhasePropsList,'fraction[wt]')) .* OxideFractions.(StablePhases{iph});
            end
        end
        Density.sol        = PhaseFractions.sol_wt/Density.sol;
        OxideFractions.sol = OxideFractions.sol./PhaseFractions.sol_wt;
    else
        Density.sol = nan;
        OxideFractions.sol(:) = nan(size(OxideFractions.(StablePhases{1})));
    end
    Density.SYS = SYSProps(strcmp(SYSPropsList,'Rho[kg/m3]'));

    % Calculate melt viscosity in [Pas] using Giordano et al., 2008
    % sort oxides from/to
    % SiO2 Al2O3 CaO MgO FeOt K2O Na2O TiO2 O Cr2O3 H2O
    % SiO2 TiO2 Al2O3 FeOt MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
    if isfield(PhaseProps,'liq')
        liq_wt = [OxideFractions.liq(1),  OxideFractions.liq(8),  OxideFractions.liq(2), ...
                  OxideFractions.liq(5),0,OxideFractions.liq(4),  OxideFractions.liq(3), ...
                  OxideFractions.liq(7),  OxideFractions.liq(6),0,OxideFractions.liq(11),0].*100;
        Viscosity.liq = grdmodel08(liq_wt,T);
    else
        Viscosity.liq = nan;
    end

    % Calculate mol from wt fractions
    mw  = [ 60.0843, 101.961276, 56.0774, 40.3044, 71.8444, 94.1960, 61.97894, 79.8658, 15.999, 151.99, 18.01528]; % molar weights
    phs = [StablePhases{:} {'sol' 'SYS'}];
    for iph = 1:length(phs)
        OxideFractions_mol.(phs{iph}) = OxideFractions.(phs{iph}) ./ mw;
        OxideFractions_mol.(phs{iph}) = OxideFractions_mol.(phs{iph}) ./ sum(OxideFractions_mol.(phs{iph}));
    end

    % Calculate Mg# for liq and SYS 
    phs = [StablePhases{:} {'SYS'}];
    for iph = 1:length(StablePhases)
        PhaseProps.(phs{iph}) = [PhaseProps.(phs{iph}),OxideFractions_mol.(phs{iph})(4) ./ (OxideFractions_mol.(phs{iph})(4)+OxideFractions_mol.(phs{iph})(5))];
    end
    PhasePropsList = [PhasePropsList 'Mg#'];
    SYSProps = [SYSProps,OxideFractions_mol.SYS(4) ./ (OxideFractions_mol.SYS(4)+OxideFractions_mol.SYS(5))];
    SYSPropsList = [SYSPropsList 'Mg#'];

    % Save output to structure
    OUT.numStablePhases(numPoint)   =  length(StablePhases);
    OUT.StablePhases(numPoint,1:length(StablePhases))  =  StablePhases;  clear StablePhases;
    OUT.StableSolutions(numPoint,1:length(StableSolutions))  =  StableSolutions;
    OUT.StablePurePhases(numPoint,1:length(StablePurePhases))  =  StablePurePhases; clear StablePurePhases;
    OUT.P(numPoint)                 =  P;              clear P;
    OUT.T(numPoint)                 =  T;              clear T;
    OUT.OxideList(numPoint,:)       =  OxideList;      clear OxideList;
    OUT.PhasePropsList(numPoint,:)  =  PhasePropsList; clear PhasePropsList;
    OUT.SYSPropsList(numPoint,:)    =  SYSPropsList;   clear SYSPropsList;
    OUT.SYSProps(numPoint,:)        =  SYSProps;       clear SYSPropsList;
    OUT.GammaList(numPoint,:)       =  GammaList;      clear GammaList;
    OUT.Gamma(numPoint,:)           =  Gamma;          clear Gamma;
    OUT.DeltaGibbsList(numPoint,1:length(StableSolutions))  =  StableSolutions;
    OUT.DeltaGibbs(numPoint,1:length(StableSolutions))      =  DeltaGibbs;     clear DeltaGibbs; clear StableSolutions;

    flds = fieldnames(PhaseFractions    ); for ifld = 1:length(flds); OUT.PhaseFractions    .(flds{ifld})(numPoint,:)  =  PhaseFractions    .(flds{ifld}); end; clear PhaseFractions;
    flds = fieldnames(Density           ); for ifld = 1:length(flds); OUT.Density           .(flds{ifld})(numPoint,:)  =  Density           .(flds{ifld}); end; clear Density;
    flds = fieldnames(Viscosity         ); for ifld = 1:length(flds); OUT.Viscosity         .(flds{ifld})(numPoint,:)  =  Viscosity         .(flds{ifld}); end; clear Viscosity;
    flds = fieldnames(EMList            ); for ifld = 1:length(flds); OUT.EMList            .(flds{ifld})(numPoint,:)  =  EMList            .(flds{ifld}); end; clear EMList;
    flds = fieldnames(EMFractions       ); for ifld = 1:length(flds); OUT.EMFractions       .(flds{ifld})(numPoint,:)  =  EMFractions       .(flds{ifld}); end; clear EMFractions;
    flds = fieldnames(SiteFractions     ); for ifld = 1:length(flds); OUT.SiteFractions     .(flds{ifld})(numPoint,:)  =  SiteFractions     .(flds{ifld}); end; clear SiteFractions;
    flds = fieldnames(OxideFractions    ); for ifld = 1:length(flds); OUT.OxideFractions    .(flds{ifld})(numPoint,:)  =  OxideFractions    .(flds{ifld}); end; clear OxideFractions;
    flds = fieldnames(OxideFractions_mol); for ifld = 1:length(flds); OUT.OxideFractions_mol.(flds{ifld})(numPoint,:)  =  OxideFractions_mol.(flds{ifld}); end; clear OxideFractions_mol;
    flds = fieldnames(PhaseProps        ); for ifld = 1:length(flds); OUT.PhaseProps        .(flds{ifld})(numPoint,:)  =  PhaseProps        .(flds{ifld}); end; clear PhaseProps;


%     PhaseData{newPoints(numPoint)}.P                    =   P; clear P;
%     PhaseData{newPoints(numPoint)}.T                    =   T; clear T;
%     PhaseData{newPoints(numPoint)}.PhaseFractions       =   PhaseFractions; clear PhaseFractions;
%     PhaseData{newPoints(numPoint)}.Density              =   Density; clear Density;
%     PhaseData{newPoints(numPoint)}.Viscosity            =   Viscosity; clear Viscosity
% 
%     PhaseData{newPoints(numPoint)}.numStablePhases      =   length(StablePhases);
%     PhaseData{newPoints(numPoint)}.StablePhases         =   StablePhases; clear StablePhases
%     PhaseData{newPoints(numPoint)}.EMList  			    =   EMList; clear EMList
%     PhaseData{newPoints(numPoint)}.EMFractions		    =   EMFractions; clear EMFractions
%     PhaseData{newPoints(numPoint)}.SiteFractions        =   SiteFractions; clear SiteFractions
%     PhaseData{newPoints(numPoint)}.OxideList  		    =   OxideList; clear OxideList
%     PhaseData{newPoints(numPoint)}.OxideFractions	    =   OxideFractions; clear OxideFractions
%     PhaseData{newPoints(numPoint)}.OxideFractions_mol   =   OxideFractions_mol; clear OxideFractions_mol
%     PhaseData{newPoints(numPoint)}.PhasePropsList  	    =   PhasePropsList; clear PhasePropsList
%     PhaseData{newPoints(numPoint)}.PhaseProps	        =   PhaseProps; clear PhaseProps
%     PhaseData{newPoints(numPoint)}.SYSPropsList  	    =   SYSPropsList; clear SYSPropsList
%     PhaseData{newPoints(numPoint)}.SYSProps	            =   SYSProps; clear SYSProps
%     PhaseData{newPoints(numPoint)}.GammaList        	=   GammaList; clear GammaList
%     PhaseData{newPoints(numPoint)}.Gamma	            =   Gamma; clear Gamma
    
    
end

fclose(fid);




