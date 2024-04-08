function Print_TC_vs_MAGEMin(Data,Data_TC)
% Prints the MAGEMin solution & the Thermocalc solution, so we can compare
% 


str = sprintf('=================================================='); disp(str)
str = sprintf('P:    %- 8.3f   kbar        ', Data.P); disp(str)
str = sprintf('T:    %- 8.3f  Celsius     ', Data.T); disp(str)
disp(' ')

if ~Data_TC.success
   % TC failed to find a solution
   Data_TC.StableFractions = Data.StableFractions*NaN;
   Data_TC.StableSolutions = Data.StableSolutions;
   Data_TC.Gibbs = NaN;
end



str = sprintf('          MAGEMin      ThermoCalc     Difference '); disp(str)
str = sprintf('          ----------   ------------   ----------  '); disp(str)
str = sprintf('Gibbs:    %- 8.6f  %-8.6f    %-8.6f        ', Data.Gibbs,         Data_TC.Gibbs,  Data.Gibbs-Data_TC.Gibbs); disp(str)
if isfield(Data_TC,'rho')
    str = sprintf('Density:  %- 8.5f  %-8.5f    %-8.6f        ', Data.Density_total, Data_TC.rho,    Data.Density_total-Data_TC.rho); disp(str)
end
disp(' ')

% Print modes
str = sprintf('mode:'); disp(str)
for i=1:length(Data.StableSolutions)
    ss = Data.StableSolutions{i};
    iTC = strcmp(Data_TC.StableSolutions,ss);
    
    str = sprintf('%5s:    %- 8.6f  %-8.6f   %- 9.7f        ',ss, Data.StableFractions(i), Data_TC.StableFractions(iTC),  Data.StableFractions(i)-Data_TC.StableFractions(iTC)); disp(str)
    
end

% print Gamma
disp(' ')
disp('Gamma: ')

disp('            SiO2        Al2O3        CaO         MgO        FeO           K2O         Na2O         TiO2         O              Cr2O3        H2O')
str='MAGEMin: ';
for i=1:length(Data.Gamma)
    str = [str, sprintf('%-11.6f,',Data.Gamma(i))]; 
end
disp(str)

if isfield(Data_TC,'Gamma')
    str='TC:      ';
    for i=1:length(Data_TC.Gamma)
        str = [str, sprintf('%-11.6f,',Data_TC.Gamma(i))];
    end
    disp(str)
    
    str='del:      ';
    for i=1:length(Data_TC.Gamma)
        str = [str, sprintf('%-12.7f ',Data.Gamma(i)-Data_TC.Gamma(i))];
    end
    disp(str)
end


% print Compositional Variables
disp(' ')
disp('Compositional Variables: ')
for iComp=1:length(Data.StableSolutions)
    disp([' ', Data.StableSolutions{iComp},':'])
    str='   MAGEMin: ';
    for i=1:length(Data.CompositionalVar{iComp})
        str = [str, sprintf(' %- 11.9f ',Data.CompositionalVar{iComp}(i))];
    end
    disp(str)
    if isfield(Data_TC,'CompositionalVar')
        str='   TC:      ';
        for i=1:length(Data_TC.CompositionalVar{iComp})
            str = [str, sprintf(' %- 11.9f ',Data_TC.CompositionalVar{iComp}(i))];
        end
        disp(str)
        
        str='   del:     ';
        for i=1:length(Data_TC.CompositionalVar{iComp})
            str = [str, sprintf(' %- 11.9f ',Data.CompositionalVar{iComp}(i)-Data_TC.CompositionalVar{iComp}(i))];
        end
        disp(str)
    end 
end


