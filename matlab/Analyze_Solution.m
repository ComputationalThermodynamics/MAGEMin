% This analyzes the solution obtained with MAGEMin @ a certain point using
% our MATLAB implementation of the database


addpath ../../matlab/11_components



if 1==1
    
    % load datapoint from the command-line
    Point           =   PseudoSectionData.PointInfo;
    StableSolutions =   Point.StableSolutions;
    StableFractions =   Point.StableFractions;
    CompositionalVar=   Point.CompositionalVar;
    CompositionalVar=   CompositionalVar(1:length(StableSolutions));        % because of a bug in reading MAGEMin data
    
    P               =   Point.P;
    T_C             =   Point.T;
    
    
    RemovePhase     =   [];
   
    if ~isempty(RemovePhase)
        StableSolutions(RemovePhase)  = [];
        CompositionalVar(RemovePhase) = [];
        StableFractions(RemovePhase)  = [];
    end
     
  
else
    % specify the data point manually
    P               =   11;
    T_C             =   1237.5;
    StableSolutions =   {'liq', 'ol', 'opx', 'cpx', 'spn'}
    
%     StableFractions =  [ 0.01285460  0.61095585  0.21773249  0.15365494  0.00480212]
    StableFractions =  [  0.160000    0.180000    0.200000    0.220000    0.240000]
    
    % liq
%     CompositionalVar{1} = [ 0.29622   0.15592   0.09264   0.03293   0.34977   0.01685   0.00130   0.01477   0.02681   0.03810 0];
    CompositionalVar{1} = [ 0.310814    0.188560   0.0905651   0.0348030    0.249634   0.0154405  0.00155156   0.0153660   0.0217085   0.0483637 0];
    
    % ol
%     CompositionalVar{2} = [  0.10167   0.00537   0.00027];
    CompositionalVar{2} = [ 0.102440  0.00576000 0.000260000];
    
    
    % opx
%     CompositionalVar{3} = [ 0.09450   0.20384   0.06854  -0.04317   0.02077   0.00812   0.02464   0.00912];
    CompositionalVar{3} = [ 0.0956900    0.192170   0.0689300 -0.0475700   0.0177700  0.00825000   0.0205800  0.00907000];
    
    % cpx 
%     CompositionalVar{4} =  [0.10231   0.13327   0.31718   0.08564  -0.03163   0.01879   0.01809   0.01335   0.00269];
    CompositionalVar{4} =  [0.0963000    0.122410    0.308960   0.0848000  -0.0282600   0.0179900   0.0154900   0.0127600  0.00280000];
    
    % spn 
%     CompositionalVar{5} =  [0.14332   0.01426   0.05443   0.00734   0.47939   0.08111   0.02259];
    CompositionalVar{5} =  [0.100000   0.0500000   0.0500000   0.0500000    0.289768  0.08111   0.02259];
    
    
    
end




T_K             =   T_C + 273.15;

%% -------------------------------------------------------------------------
disp('====================================================================')
disp(' Analyzing solution for the following parameters ')
disp('====================================================================')
disp(['Pressure [kbar] : ',num2str(P)]);
disp(['Temperature [C] : ',num2str(T_C)]);
str_name = 'Stable phases   : ';
str_val  = 'Fractions       : ';
nMin = length(StableSolutions);
for iMin=1:nMin
    str_name = [str_name, sprintf('%8s\t'  ,StableSolutions{iMin})];
    str_val  = [str_val,  sprintf('%8.6f\t',StableFractions(iMin))];
end
disp(str_name)
disp(str_val)

% Print endmember composition
disp(' ')
for iMin=1:nMin
    disp(StableSolutions{iMin})
    str_val  = sprintf('  CompVar:\t');
    for j=1:length(CompositionalVar{iMin})
        str_val  = [str_val,  sprintf('%8.6f\t',CompositionalVar{iMin}(j))];
    end
    disp(str_val)
    disp(' ')
end




%% -------------------------------------------------------------------------
disp('====================================================================')
disp(' STEP 1 - comparing with MATLAB database implementation')
disp('====================================================================')

NameString  = [];
DatabaseName  	=   'Holland2018_Database';                 % database employed here
Minerals    =  feval(DatabaseName,'ListEntries');
xStart      =  [];

br              =   [ 38.451   1.774    2.821   50.510    5.880    0.010    0.250    0.100    0.096    0.109  0]';
                    
   
br              =   br/sum(br);

for iMin=1:length(StableSolutions)       %% This should all be moved to a database
    Minerals.(StableSolutions{iMin})      =  feval(DatabaseName,'Initialize',Minerals,StableSolutions{iMin},P,T_K); % Initialize solution models & endmembers
    xStart = [xStart, CompositionalVar{iMin}];    % compositional variables
    switch  Minerals.(StableSolutions{iMin}).type
        case 'SolutionModel'
            NameString  	=   [NameString, Minerals.(StableSolutions{iMin}).variables];
    end
end
%
% -------------------------------------------------------------------------
PrintOutput = true;
x0          =   [StableFractions, xStart];
[G_system, dGsys, X_system, X_sol, sf, muTotal, alpha, apf, fbc, G_sol, prop, gbase, IDM, mu_Gex, gradients, sf_struct, CompEndmembers] = G_System(x0,br,DatabaseName, Minerals, P, T_K, PrintOutput);



% Print endmember composition
disp('________________________________________________________')
for iMin=1:nMin
    disp(StableSolutions{iMin})
    str_name = sprintf('     \t');
    str_val  = sprintf('Init:\t');
    for j=1:length(Minerals.(StableSolutions{iMin}).variables)
        str_name = [str_name, sprintf('%8s\t'  ,Minerals.(StableSolutions{iMin}).variables{j})];
        str_val  = [str_val,  sprintf('%8.6f\t',CompositionalVar{iMin}(j))];
    end
    disp(str_name)
    disp(str_val)
    disp(' ')
end



disp(['G_MATLAB: ',num2str(G_system)])
disp(['G_MAGEMin: ',num2str(Point.Gibbs)])





%% -------------------------------------------------------------------------
disp(' ');
disp(' ');
disp('====================================================================')
disp(' STEP 2 - check how well we fit the constraints using the MAGEMin solution')
disp('====================================================================')


% All potential reactions are given by the nullspace
EndmemberComposition    =   [];
MineralsInvolved        =   [];
Range                   =   [];
StartGuess              =   [];
nMin                    =   length(StableSolutions); % # of minerals considered
for j = 1:nMin
    X_endmember             =   Minerals.(StableSolutions{j}).Composition_Endmember;
    EndmemberComposition    =   [EndmemberComposition; X_endmember];                                  % Keep track of chemical composition of the endmembers
    MineralsInvolved        =   [MineralsInvolved; ones(size(X_endmember,1),1)*j];                        % keep track of which endmembers are involved
    if strcmp(Minerals.(strip_char(StableSolutions{j})).type,'SolutionModel')
        Range                   = [Range; Minerals.(strip_char(StableSolutions{j})).Range];
        StartGuess              = [StartGuess; Minerals.(strip_char(StableSolutions{j})).StartingGuess];
    end
end

PotentialReactions = null(EndmemberComposition','r');

% Yet, this also includes reactions involving only 1 components; filter them out
id = [];
for i=1:size(PotentialReactions,2)
    if length(find(abs(PotentialReactions(:,i))>0))>1
        id = [id, i];           % reaction involving at least 2 components
    end
end
PotentialReactions = PotentialReactions(:,id);


ComputeGradients = true;

% constraints
[c,ceq,grad_c,grad_ceq, gradients] = constraints_Thermocalc(x0, br, DatabaseName, Minerals, P, T_K, PotentialReactions, MineralsInvolved, ComputeGradients);

disp('Mass conservation constraints:')
disp(ceq(1:11)')
disp(['sum(alpha) = ',num2str(1+ceq(12))])
disp('mu equations: ')
disp(ceq(13:end))

MassCons_MAGEMin = sum(abs(ceq( 1:11)));
SumdMu_MAGEMin   = sum(abs(ceq(13:end)));





%% -------------------------------------------------------------------------
disp(' ');
disp(' ');
disp('====================================================================')
disp(' STEP 3 - solve it Thermocalc style')
disp('====================================================================')
pause
P       =   Point.P;
T       =   Point.T;

PrintOutput = false;

[obj, Success, xOpt_TC, G_TC, Fraction_TC, CompVar_TC] =   Thermocalc_minimization(T_K,P, StableSolutions,StableFractions, CompositionalVar, 1, PrintOutput);

[c,ceq,grad_c,grad_ceq, gradients] = constraints_Thermocalc(xOpt_TC, br, DatabaseName, Minerals, P, T_K, PotentialReactions, MineralsInvolved, ComputeGradients);

disp('Mass conservation constraints after TC optimization:')
disp(ceq(1:11)')
disp(['sum(alpha)-1 = ',num2str(ceq(12))])
disp('mu equations: ')
disp(ceq(13:end))

MassCons_TC = sum(abs(ceq( 1:11)));
SumdMu_TC   = sum(abs(ceq(13:end)));

disp('________________________________________________________')
disp(['Gibbs MAGEMin    : ' , num2str(Point.Gibbs,'%10.7f')])
disp(['Gibbs Matlab DB : ' , num2str(G_system,'%10.7f')])
disp(['Gibbs Matlab TC : ' , num2str(G_TC,'%10.7f')])
disp(' ')


str_name  = 'Stable phases   : ';
str_val1  = 'Frac MAGEMin     : ';
str_val2  = 'Frac Matlab TC  : ';
nMin = length(StableSolutions);
for iMin=1:nMin
    str_name = [str_name, sprintf('%8s\t'   , StableSolutions{iMin})];
    str_val1 = [str_val1,  sprintf('%8.6f\t', StableFractions(iMin))];
    str_val2 = [str_val2,  sprintf('%8.6f\t', xOpt_TC(iMin))];
    
end
disp(str_name)
disp(str_val1)
disp(str_val2)
disp(' ')
for iMin=1:nMin
    disp(StableSolutions{iMin})
    str_name = sprintf('           \t');
    str_val1  = sprintf('MAGEMin    :\t');
    str_val2  = sprintf('Matlab TC :\t');
    for j=1:length(Minerals.(StableSolutions{iMin}).variables)
        str_name  = [str_name, sprintf('%8s\t'   , Minerals.(StableSolutions{iMin}).variables{j})];
        str_val1  = [str_val1,  sprintf('%8.6f\t', CompositionalVar{iMin}(j))];
        str_val2  = [str_val2,  sprintf('%8.6f\t', CompVar_TC{iMin}(j))];
    end
    disp(str_name)
    disp(str_val1)
    disp(str_val2)
    disp(' ')
end





%% -------------------------------------------------------------------------
disp(' ');
disp(' ');
disp('====================================================================')
disp(' STEP 4- solve it with Boris-type code optimization')
disp('====================================================================')
pause
% x0=xOpt_TC;


lb = [ones(length(StableSolutions),1)*0;  Range(:,1);];
ub = [ones(length(StableSolutions),1)*1;  Range(:,2);];

options     = optimoptions('fmincon','MaxFunEvals',1e6,'MaxIter',200,'TolX',1e-18,'TolCon',1e-7,'TolFun',1e-7,'TolProjCG',1e-7,'GradObj','on', 'Display','iter','SpecifyConstraintGradient',true);
% options     = optimoptions(options,'Algorithm','interior-point');       % works much better it seems
%         options     = optimoptions(options,'Algorithm','sqp');
%         options     = optimoptions(options,'Algorithm','sqp-legacy');

%     options     = optimoptions(options,'Algorithm','active-set');
%     options     = optimoptions(options,'Display','iter');

Compute_Gradients                =   logical(1);
[xOpt_FMINCON,fval,EXITFLAG,OUTPUT]       = fmincon(@(x)ObjectiveFunction_Equilibrium(x,br,DatabaseName, Minerals, P, T_K, PotentialReactions),x0,[],[],[],[],lb,ub, @(x)mycon2(x, br, DatabaseName, Minerals, P, T_K,PotentialReactions,MineralsInvolved,Compute_Gradients),options);
[G_FMINCON, ~, X_system_optimized, X_sol, sf, muTotal] = G_System(xOpt_FMINCON,br,DatabaseName, Minerals, P, T_K,logical(0));
Fraction_FMINCON        =   xOpt_FMINCON(1:nMin);


% [c,ceq,grad_c,grad_ceq, gradients] = mycon2(x0, br, DatabaseName, Minerals, P, T_K,PotentialReactions,MineralsInvolved,Compute_Gradients)

xCompVar                =   xOpt_FMINCON(nMin+1:end);     % Compositional Variables
num                     =   1;
for i=1:nMin
    CompVar_FMINCON{i}  =   xCompVar(num:num+length(CompositionalVar{i})-1);
    num                 =   num + length(CompositionalVar{i});
end

[c,ceq,grad_c,grad_ceq, gradients] = constraints_Thermocalc(xOpt_FMINCON, br, DatabaseName, Minerals, P, T_K, PotentialReactions, MineralsInvolved, ComputeGradients);
MassCons_FMINCON = sum(abs(ceq( 1:11)));
SumdMu_FMINCON   = sum(abs(ceq(13:end)));


disp('________________________________________________________')
disp(['Pressure [kbar] : ',num2str(P)]);
disp(['Temperature [C] : ',num2str(T_C)]);
disp(' ')
str_name = 'Stable phases   : ';
str_val  = 'Fractions       : ';
nMin = length(StableSolutions);
for iMin=1:nMin
    str_name = [str_name, sprintf('%8s\t'  ,StableSolutions{iMin})];
    str_val  = [str_val,  sprintf('%8.6f\t',StableFractions(iMin))];
end
disp(str_name)
disp(str_val)
disp(' ')

disp(['Gibbs MAGEMin         : ' , num2str(Point.Gibbs,'%10.7f')])
disp(['Gibbs Matlab DB      : ' , num2str(G_system,'%10.7f')])
disp(['Gibbs Matlab TC      : ' , num2str(G_TC,'%10.7f')])
disp(['Gibbs Matlab FMINCON : ' , num2str(G_FMINCON,'%10.7f')])

disp(' ')


str_name  = 'Stable phases        : ';
str_val1  = 'Frac MAGEMin          : ';
str_val2  = 'Frac Matlab TC       : ';
str_val3  = 'Frac Matlab FMINCON  : ';

nMin = length(StableSolutions);
for iMin=1:nMin
    str_name = [str_name, sprintf('%8s\t'   , StableSolutions{iMin})];
    str_val1 = [str_val1,  sprintf('%8.6f\t', StableFractions(iMin))];
    str_val2 = [str_val2,  sprintf('%8.6f\t', xOpt_TC(iMin))];
    str_val3 = [str_val3,  sprintf('%8.6f\t', Fraction_FMINCON(iMin))];
end
disp(str_name)
disp(str_val1)
disp(str_val2)
disp(str_val3)
disp(' ')

for iMin=1:nMin
    disp(StableSolutions{iMin})
    str_name = sprintf('                 \t');
    str_val1  = sprintf('MAGEMin         :\t');
    str_val2  = sprintf('Matlab TC      :\t');
    str_val3  = sprintf('Matlab FMINCON :\t');
    for j=1:length(Minerals.(StableSolutions{iMin}).variables)
        str_name  = [str_name,  sprintf('%8s\t'   , Minerals.(StableSolutions{iMin}).variables{j})];
        str_val1  = [str_val1,  sprintf('%8.6f\t', CompositionalVar{iMin}(j))];
        str_val2  = [str_val2,  sprintf('%8.6f\t', CompVar_TC{iMin}(j))];
        str_val3  = [str_val3,  sprintf('%8.6f\t', CompVar_FMINCON{iMin}(j))];
    end
    disp(str_name)
    disp(str_val1)
    disp(str_val2)
    disp(str_val3)
    disp(' ')
end
disp(['sum(|dMu|) MAGEMin     : ' , num2str(SumdMu_MAGEMin,'%10.7e')])
disp(['sum(|dMu|) TC         : ' , num2str(SumdMu_TC,'%10.7e')])
disp(['sum(|dMu|) FMINCON    : ' , num2str(SumdMu_FMINCON,'%10.7e')])
disp(' ')
disp(['sum(Mass) MAGEMin      : ' , num2str(MassCons_MAGEMin,'%10.7e')])
disp(['sum(Mass) TC          : ' , num2str(MassCons_TC,'%10.7e')])
disp(['sum(Mass) FMINCON     : ' , num2str(MassCons_FMINCON,'%10.7e')])

%% -------------------------------------------------------------------------
disp(' ');
disp(' ');
disp('====================================================================')
disp(' STEP 5 - Compute Gamma and find out whether there are more stable minima ')
disp('====================================================================')
pause

[Pi_MAGEMin,     Gamma_MAGEMin]   = ComputeDrivingForce_Phases(x0,            br,DatabaseName, Minerals, P, T_K);
[Pi_TC,         Gamma_TC]       = ComputeDrivingForce_Phases(xOpt_TC,       br,DatabaseName, Minerals, P, T_K);
[Pi_FMINCON,    Gamma_FMINCON]  = ComputeDrivingForce_Phases(xOpt_FMINCON,  br,DatabaseName, Minerals, P, T_K);


disp('')
ThermoCalcNames = {'SiO2' 'Al2O3' 'CaO' 'MgO' 'FeO' 'K2O' 'Na2O' 'TiO2' 'O', 'Cr2O3','H2O'};

disp('---')
disp('Driving force of the optimal solution :')
disp('Gamma =')

str         =   sprintf('\t\t');
str_val1    =   sprintf('Gamma_MAGEMin :\t');
str_val2    =   sprintf('Gamma_TC     :\t');
str_val3    =   sprintf('Gamma_FMINCON:\t');

for i=1:length(ThermoCalcNames)
    str         = [str,         sprintf('%8s\t',    ThermoCalcNames{i})];
    str_val1    = [str_val1,    sprintf('%3.3f\t',  Gamma_MAGEMin(i))];
    str_val2    = [str_val2,    sprintf('%3.3f\t',  Gamma_TC(i))];
    str_val3    = [str_val3,    sprintf('%3.3f\t',  Gamma_FMINCON(i))];
    
end
disp(str)
disp(str_val1)
disp(str_val2)
disp(str_val3)


disp(' ')
disp('Driving force:')
str         =   sprintf('\t\t');
str_val1    =   sprintf('Pi_MAGEMin :\t');
str_val2    =   sprintf('Pi_TC     :\t');
str_val3    =   sprintf('Pi_FMINCON:\t');
for i=1:length(Pi_MAGEMin)
    str         = [str,         sprintf('%8s\t',    StableSolutions{i})];
    str_val1    = [str_val1,    sprintf('% 3.3e\t',  Pi_MAGEMin(i))];
    str_val2    = [str_val2,    sprintf('% 3.3e\t',  Pi_TC(i))];
    str_val3    = [str_val3,    sprintf('% 3.3e\t',  Pi_FMINCON(i))];
end
disp(str)
disp(str_val1)
disp(str_val2)
disp(str_val3)

