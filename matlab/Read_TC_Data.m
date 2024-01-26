function Data_TC = Read_TC_Data(File_dir,File_out, has_H2O)
%
% This reads a TC output file to determine if the calculation was a 
% success or not; if it was a success it reads the key data

cur_dir = pwd;
cd(File_dir)

FileName        = [File_out(1:end-4),'-o.txt'];
Data_TC.success = IsSuccess(FileName);

if Data_TC.success
    % Calculation was a success; compute stable mode & proportions
    [StableSolutions,StableFractions, Gibbs, rho,Gamma_names, Gamma] = ModeGibbsFractions(FileName);
    
    if has_H2O
        Gamma = [Gamma(2:11); Gamma(1)];
        Gamma_names1 = Gamma_names(2:11);
        Gamma_names1(11)    = Gamma_names(1);
        Gamma_names = Gamma_names1;
    end
    
    Data_TC.StableSolutions =   StableSolutions';
    Data_TC.StableFractions =   StableFractions;
    Data_TC.Gibbs           =   Gibbs;
    Data_TC.rho             =   rho;
    Data_TC.Gamma           =   Gamma(:);
    Data_TC.Gamma_names     =   Gamma_names;
    
    % Extract melt fraction
    ind_liq = find(ismember(StableSolutions,'liq'));
    if ~isempty(ind_liq)>0
        liq = StableFractions(ind_liq);
    else
        liq = 0;
    end
    Data_TC.liq             =   liq;
    Data_TC.numStablePhases =   length(StableFractions);
    
end

% Next, read the compositional variables
FileName        = [File_out(1:end-4),'-ic.txt'];
if Data_TC.success
    % Calculation was a success; compute stable mode & proportions
    [CompositionalVar] = CompositionalVariables(FileName, StableSolutions);
        
    ind_Liq                         = find(strcmp(StableSolutions,'liq'));
    if ~isempty(ind_Liq)
        if length(CompositionalVar{ind_Liq})==10
            CompositionalVar{ind_Liq}  = [CompositionalVar{ind_Liq} 0];   % because of lacking H2O in our TC database
        end
    end
    
    Data_TC.CompositionalVar = CompositionalVar;
end

if Data_TC.success
    % Sort according to names (makes things easier to compare later)
    [Data_TC.StableSolutions,iS]    =   sort(Data_TC.StableSolutions);
    Data_TC.CompositionalVar        =   Data_TC.CompositionalVar(iS);
    Data_TC.StableFractions         =   Data_TC.StableFractions(iS);

end
cd(cur_dir);



%--------------------------------------------------------------------------
function success = IsSuccess(FileName)
% Determine if the TC calculation was a success or not
fid = fopen(FileName);

success = logical(1);
while ~feof(fid)
    st = fgetl(fid);
    
    if regexp(st,'No solution') 
        success = logical(0);
    end
    if regexp(st,'bombed')>0
        success = logical(0);
    end
    
end
fclose(fid);



%--------------------------------------------------------------------------
function [StableSolutions,StableFractions, Gibbs, rho, Gamma_names, Gamma ] = ModeGibbsFractions(FileName)
% Determine the stable mode

fid = fopen(FileName);
while ~feof(fid)
    st = fgetl(fid);
    
    if regexp(st,'mode ')
        % 1) Names
        str_names           =   strsplit(st);       % all names on this line
        str_start           =   find(strcmp(str_names,'mode'));
        str_end             =   find(strcmp(str_names,'G'));
        StableSolutions     =   str_names(str_start+1:str_end-1);
        
        
        % 2) Fractions
        str_next            =   fgetl(fid);
        if regexp(str_next,'#')
            ind =  strfind(str_next,'-');
            if length(ind)>0
                str_next = [str_next(1:ind(1)-1),' ',str_next(ind(1):end)];
            end
            
            str_next1 = strsplit(str_next);
            % This is a case where we have no stable solution as one of the
            % mass fractions is negative
           
            for i=2:length(str_next1)
                mode_num(i-1) = str2num(str_next1{i});
            end
            mode_num(end+1) = NaN; % gibbs is not defined in this case
            mode_num(end+1) = NaN; % neither is density
            
        else
            mode_num            =   str2num(str_next); % all numeric values on next line
        end
   
        StableFractions     =   mode_num(1:length(StableSolutions));
        Gibbs               =   mode_num(length(StableSolutions)+1);    % Gibbs free energy
        rho                 =   mode_num(end)*1e3;                      % density
    end
    
    if regexp(st,'mu')
        % determine Gamma from TC
        str_names           =   strsplit(st);       % all names on this line
        str_start           =   find(strcmp(str_names,'mu'));
        Gamma_names         =   str_names(str_start+1:end);
        Gamma               =   str2num(fgetl(fid)); % all numeric values on next line
        Gamma = [Gamma(:); 0];
        
    end

    

    
end
fclose(fid);





%--------------------------------------------------------------------------
function [CompositionalVar] = CompositionalVariables(FileName, StableSolutions)
% Determine the stable mode

iSol = 1;
fid = fopen(FileName);
% skip first lines
for i=1:4
      st = fgetl(fid);
end

CompositionalVar=[];
for i=1:length(StableSolutions)
    CompositionalVar{i} = [];
end

while ~feof(fid)
    st = fgetl(fid);
    
    if regexp(st,'site fractions')
        iSol=length(StableSolutions)+1;
    end
    
    if iSol<= length(StableSolutions)
        CompVar=[];
        if regexp(st,StableSolutions{iSol})
            % 1) Get values on next line
            str_next            =   fgetl(fid);
            
            str_next1 = strsplit(str_next);
            % This is a case where we have no stable solution as one of the
            % mass fractions is negative
            
            for i=2:length(str_next1)
                CompVar(i-1) = str2num(str_next1{i});
            end
            
            CompositionalVar{iSol}     =   CompVar;
            iSol = iSol+1;
        end
    end
    

    
end
fclose(fid);

