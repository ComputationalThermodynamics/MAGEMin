function  [variableNames, InitialGuess] = SolidSolution_VarNames(SS)

InitialGuess = [];

switch SS
    case 'liq'          % liquids
        variableNames = {'wo(L)','sl(L)','fo(L)','fa(L)','jd(L)','hm(L)','ek(L)','ti(L)','kj(L)','yct(L)','h2o(L)'};
        InitialGuess  = [0.2, 0.2, 0.1, 0.1, 0.05, 0.001, 0.001, 0.001, 0.001,  0.001 0.001];
        
    case 'pli'           % Plagioclase (1st model from Holland and Powell, Ibar1, for anorthite-rich composition)
        
        %proportion
        variableNames = {'ca(pli)','k(pli)'};
        InitialGuess  = [0.8, 0.03];
        
    case 'plc'          % Plagioclase (2nd model from Holland and Powell, Cbar1, for albite-rich composition)
        
        %proportion
        variableNames = {'ca(plc)','k(plc)'};
        
    case 'ol'           % Olivine
        
        %proportion
        variableNames = {'x(ol)','c(ol)','Q(ol)'};
        InitialGuess  = [0.1, 0.002 0.01];

    case 'ksp'          % k-feldspar
        
        % proportion
        variableNames = {'na(ksp)','ca(ksp)'};
        InitialGuess  = [0.1, 0.001];

    case 'g'            % Garnet
        
        %proportion
        variableNames = {'x(g)','c(g)','f(g)','cr(g)','t(g)'};
        InitialGuess  = [0.455851 0.286243 0.00253871 0.000841005 0.0113394];
        
    case 'opx'          % Orthopyroxene
        
        %proportion
        variableNames = {'x(opx)','y(opx)','c(opx)','Q(opx)','f(opx)','t(opx)','cr(opx)','j(opx)'};
        InitialGuess  = [0.05, 0.006, 0.025, 0.032, 0.001, 0.001, 0.001, 0.001];	

    case 'cpx'
        variableNames = {'x(cpx)','y(cpx)','o(cpx)','n(cpx)','Q(cpx)','f(cpx)','cr(cpx)','t(cpx)','k(cpx)'};
        InitialGuess  = [0.075 0.1120 0.05 0.11 -0.0005, 0.001, 0.001, 0.001, 0.001];	

    case 'pig'
        variableNames = {'x(pig)','y(pig)','o(pig)','n(pig)','Q(pig)','f(pig)','cr(pig)','t(pig)','k(pig)'};
        InitialGuess  = [0.124, 0.1120, 0.88, 0.028, -0.0115, 0.004, 0.001, 0.148, 0.001];

    case 'amp'
        %proportion
        variableNames = {'x(amp)','y(amp)','z(amp)','a(amp)','k(amp)','c(amp)','f(amp)','t(amp)','Q1(amp)','Q2(amp)'};
        InitialGuess  = [0.3, 0.2, 0.01, 0.45, 0.01, 0.8, 0.05, 0.01, -0.01, 0.1];

    case 'ilm'          % Ilmenite
        %proportion
        variableNames = {'x(ilm)','Q(ilm)'};
        InitialGuess  = [0.8, 0.055];

    case 'spl'
        variableNames = {'x(spl)','y(spl)','c(spl)','t(spl)','Q1(spl)','Q2(spl)','Q3(spl)'};
        InitialGuess  = [0.20, 0.1, 0.8, 0.05, 0.05, 0.05, 0.05];

    case 'cm'
        variableNames = {'x(cm)','y(cm)','c(cm)','t(cm)','Q1(cm)','Q2(cm)','Q3(cm)'};
        InitialGuess  = [0.20, 0.1, 0.8, 0.05, 0.05, 0.05, 0.05];

    case 'mt'
        variableNames = {'x(mt)','y(mt)','c(mt)','t(mt)','Q1(mt)','Q2(mt)','Q3(mt)'};
        
    case 'mu'
        variableNames = {'x(mu)','y(mu)','f(mu)','n(mu)','c(mu)'};
        InitialGuess  = [0.25, 0.6, 0.17, 0.06, 0.004];

    case 'bi'
        variableNames = {'x(bi)','y(bi)','f(bi)','t(bi)','Q(bi)'};
        InitialGuess  = [0.35, 0.25, 0.04, 0.17, 0.25];
        
    case {'ru','q','lc','sph','per','stv','ne','coe','ep','H2O','law','cor','ky'}
        variableNames=[];
    otherwise
        error('Not yet programmed %s',SS )
end
