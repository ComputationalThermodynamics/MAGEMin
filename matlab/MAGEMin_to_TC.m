function [Data_TC, Data] = MAGEMin_to_TC(varargin)

File_dir = varargin{1};
File_out = varargin{2};
if nargin>=3
    % We supply all we need already
    Data = varargin{3};
else
    % Read MAGEMin Data from disk
    [PhaseData, Status] = ReadData_MAGEMin(1,[],0);
    Data =  PhaseData{1};
end
if nargin>=4
    TC_Output = varargin{4};
else
    TC_Output = logical(1);
end

if nargin>5
    ComputeLine =   varargin{5}; 
    LineData    =   varargin{6};
else
    ComputeLine = false;
    LineData = [];
end


if ~isfield(Data,'Chemistry')
    Data.Chemistry.Predefined = 0;  % test 0 is default
end

str = '     Testing solutions: ';
for i=1:length(Data.StableSolutions)
   str = [str, sprintf(' %s ',Data.StableSolutions{i})]; 
end
str = [str, ' ... '];
fprintf(str); 


% Write TC file
cur_dir = pwd;
cd(File_dir);

% Write TC output file, which can be run with TC (using MAGEMin xEOS)
has_H2O = Write_TC_file(Data,File_out, ComputeLine, LineData);

if TC_Output
    disp(['Wrote ThermoCalc input file ',File_out,' in directory ',File_dir]);
end
% run TC

if TC_Output
    system_call =  './tc350beta';
else
    system_call =  './tc350beta > /dev/null 2>&1';
end
tic;
system(system_call);

if ComputeLine
    if exist('tc-testMAGEMin.csv','file')>0
        % Read data of the line
        Data_TC.Table  = readtable('tc-testMAGEMin.csv','headerlines',1);
        Data_TC.P      = table2array(Data_TC.Table(:,1));
        Data_TC.T      = table2array(Data_TC.Table(:,2));
    else
        Data_TC.P =[];
        Data_TC.T =[];
    end
end


cd(cur_dir);

if ~ComputeLine
    % read TC data
    Data_TC=Read_TC_Data('ThermoCalc','tc-testMAGEMin.txt', has_H2O);
    
    
    if ~Data_TC.success
        str=[' failed ... '];
        fprintf(str)
    end

end

TCtime_ms = toc*1000;

str=[' took ',num2str(TCtime_ms),' ms \n'];
fprintf(str)



%==========================================================================
function has_H2O = Write_TC_file(Data,File_out, ComputeLine, Line)

Chem=table2array(Data.Chemistry.OxProp(:,2));
if Chem(end)>0
    has_H2O = true;
else
    has_H2O = false;
end

fid = fopen(File_out,'w');
fprintf(fid,'%% ==================================================================\n');
fprintf(fid,'%% Scriptfile for calcs in KNCFMASTOCr (dry mantle melting), generated from MAGEMin \n');
fprintf(fid,'%% ==================================================================\n');
fprintf(fid,'\n');
if has_H2O
    fprintf(fid,'axfile  ig50NCKFMASHTOCr \n');
else
    fprintf(fid,'axfile  ig50NCKFMASTOCr \n');
end

fprintf(fid,'\n');
fprintf(fid,'autoexit yes \n');
if ~ComputeLine
    fprintf(fid,'debuglevel 2   %% replaces infolevel \n');
end
fprintf(fid,'\n');
fprintf(fid,'c8 yes\n');
fprintf(fid,'\n');

fprintf(fid,'with ');
for i=1:length(Data.StableSolutions)
    [variableNames, InitialGuess] = SolidSolution_VarNames(Data.StableSolutions{i});
%     if length(variableNames)>0
        % if it is a solid solution and not an EM
        fprintf(fid,' %s',Data.StableSolutions{i});
%     end
end
if ~isempty(Line)
    if ~contains(Line.PhaseDisappear{1},Data.StableSolutions)
        fprintf(fid,' %s',Line.PhaseDisappear{1});
    end
end
fprintf(fid,'\n');
fprintf(fid,'inexcess %% list of excess phases \n');
fprintf(fid,'\n');

fprintf(fid,'%% =======================\n');
fprintf(fid,'%% conditions\n');
fprintf(fid,'%% =======================\n');
fprintf(fid,'pseudosection \n');
fprintf(fid,'\n');

fprintf(fid,'diagramPT 0.001 60.001 700  2000 \n');
if ~ComputeLine
    fprintf(fid,'calcT %f \n',Data.T);
    fprintf(fid,'calcP %f \n',Data.P);
else
    if Line.CalcTatP
        fprintf(fid,'calcT %f \n',Data.T);
        Val = [Line.P_val, diff(Line.P_val)/Line.numPoints];
        fprintf(fid,'calcP %f %f %f \n',Val(1), Val(2), Val(3));
        fprintf(fid,'calcTatP \n',Data.P);
    
    else
        fprintf(fid,'calcP %f \n',Data.P);
        Val = [Line.T_val, diff(Line.T_val)/Line.numPoints];
        fprintf(fid,'calcT %f %f %f \n',Val(1), Val(2), Val(3));
        fprintf(fid,'calcTatP no \n',Data.P);
    end
    
end


fprintf(fid,'\n');
fprintf(fid,'%% ------------------------------------------------  \n');

Chem=table2array(Data.Chemistry.OxProp(:,2));
if has_H2O
    % with H2O
    fprintf(fid,'bulk H2O SiO2 Al2O3   CaO    MgO      FeOt  K2O    Na2O    TiO2   O       Cr2O3  \n');
    imax = 11;
else
    % no H2O
    fprintf(fid,'bulk SiO2     Al2O3   CaO    MgO      FeOt  K2O    Na2O    TiO2   O       Cr2O3  \n');
    imax = 10;
end
fprintf(fid,'bulk ');
if ~has_H2O
    for i=1:imax
        fprintf(fid,' %f ',Chem(i));
    end
else
    % H2O should go first
    fprintf(fid,' %f ',Chem(end));
    for i=1:imax-1
        fprintf(fid,' %f ',Chem(i));
    end
end

fprintf(fid,' \n');

fprintf(fid,'%% ------------------------------------------------  \n');
fprintf(fid,'%% ======================= \n');
fprintf(fid,'%% which calcs \n');
fprintf(fid,'%% ======================= \n');
fprintf(fid,'dogmin no  %% do G minimisation. CARE!! if bad starting guesses, stable assemblage will not be found! use 0-2 as if you want to do minimization \n');
fprintf(fid,'%%maxvar 2 8     %% min/max variance of calcs to try    \n');
fprintf(fid,'modeisopleth no \n');
if ~ComputeLine
    fprintf(fid,'zeromodeisopleth no \n');
else
    fprintf(fid,'zeromodeisopleth %s  \n',Line.PhaseDisappear{1});
end

fprintf(fid,'%% zeromodeisopleth spn %% phase that disappears (e.g. @ univariant line)\n'); 
if ~ComputeLine
    fprintf(fid,'calcmu yes \n');
end
fprintf(fid,'project no \n');
fprintf(fid,'%% ======================= \n');
fprintf(fid,'%% set up a-x relations \n');
fprintf(fid,'%% ======================= \n');
fprintf(fid,'\n');
fprintf(fid,'%% --------------------------------------------------------  \n');
fprintf(fid,'%% Some sample starting guesses for the compositional  \n');
fprintf(fid,'%% variables of phases.  \n');
fprintf(fid,'%% --------------------------------------------------------  \n');

for i=1:length(Data.StableSolutions)
    
    [variableNames, InitialGuess] = SolidSolution_VarNames(Data.StableSolutions{i});

    if length(variableNames)>0
        num             = Data.CompositionalVar{i};
        
        if 1==1
           % do not pass exactly 0 to TC:
           num(find(abs(num)<1e-4)) = 1e-4;
        end
        
        for j=1:length(variableNames)
            fprintf(fid,'xyzguess %8s   %20.12f \n',variableNames{j},num(j));
        end
    end
    fprintf(fid,'%% -----------------  \n');
    
end
fprintf(fid,'%% --------------------------------------------------------  \n');

fclose(fid);





