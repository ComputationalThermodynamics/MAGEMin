function [surface] = SampleSolutionSpace(SS,T,P,Gamma, xEOS, varargin)
%
% This runs MAGEMin for a single solution model and a single xEOS, or for a
% range of points. It samples the solution space, setting values to NaN
% that do not fullfill the SF constraints
%
%
%

if nargin==5
    ind         = [1 2];
    lb          = xEOS(1:2);
    ub          = xEOS(1:2);
    nPoints     = [25 30];
    
else
    
    
    x0          =   xEOS;               % origin (will be at (0,0)!) 
    Normal      =   varargin{1};        % normal vector that indicates in which direction we sample
    lb          =   varargin{2};        % lower bound (in normalized sense)    
    ub          =   varargin{3};        % upper bound    
    nPoints     =   varargin{4};        % number of sample points/direction
    
end

% id = find(ub==1); ub(id) = ub(id)-1e-3;
% id = find(lb==0); lb(id) = lb(id)+1e-3;


%  Normalize normal vector to ensure that values don't go beyond bounds
Normal(1,:) = Normal(1,:); %/(max(Normal(1,:))*20);
Normal(2,:) = Normal(2,:); %/(max(Normal(2,:))*20);



% Create 2D grid
db    = (ub-lb)./(nPoints-1);

[x,y] = meshgrid(lb(1):db(1):ub(1), lb(2):db(2):ub(2));
if isempty(x)
    % in case we only consider one point
    x   =   xEOS(ind(1));
    y   =   xEOS(ind(2));
end

% Create xEOS for all points
num = 1;
for i=1:size(x,1);
    for j=1:size(x,2)
        xEOS_in = x0 + x(i,j)*Normal(1,:) + y(i,j)*Normal(2,:);

        
        PhaseData{num}.StartingValues_xEOS.CompositionalVar{1} = xEOS_in;
        PhaseData{num}.StartingValues_xEOS.SolutionModels{1}   = SS;
        PhaseData{num}.Gamma  =   Gamma;
        
        xEOS_vec(num,:) = xEOS_in;
        
        num = num+1;
    end
end

% Write input file for all points
n_points = Write_MAGEMin_InputFile(PhaseData, T, P, Gamma);

% Run MAGEMin for all points (note: only on 1 core!)
system('rm output/_thermocalc_style_output.txt')
str = ['./MAGEMin --Mode=1 --File=MAGEMin_input.dat --Verb=1  --n_points=',num2str(n_points)];
disp(str)
system(str);


% Read output for all points
[df_full, prop_EM, sf] = Read_OUTPUT();

df = df_full;
ind = find(any(sf>0,2));
if ~isempty(ind)
    df(ind) = NaN;
end


% resize df
df      = reshape(df,[size(x')])';          % here "forbidden" points are set to NaN
df_full = reshape(df_full,[size(x')])';     % this includes forbidden points that have sf<0


% surface.xP = reshape(prop_EM(:,ind(1)),[size(x')])';
% surface.yP = reshape(prop_EM(:,ind(2)),[size(x')])';


% Store data
surface.prop_EM     =   prop_EM;
surface.xEOS_vec    =   xEOS_vec;
surface.sf          =   sf;
surface.x           =   x;
surface.y           =   y;
surface.z           =   df;
surface.df_full     =   df_full;





%--------------------------------------------------------------------------
% Write input file for MAGEMin and sample drving force for this point
function n_points = Write_MAGEMin_InputFile(PhaseData, T, P, Gamma)


% newPoints           =   1;
Use_Gamma           =   true;
Use_xEOS            =   true;
fid                 =   fopen('MAGEMin_input.dat','w');
n_points            =   length(PhaseData);
for i=1:n_points
    % Read line with P/T and Gamma
    id      =   i;
    
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
    if ~Use_xEOS
        n_phases = 0;
    else
        if ~isempty(PhaseData)
            if isfield(PhaseData{id},'StartingValues_xEOS')
                % We only do this if we have specified StartingValues for a specific point
                % They may not exists in the first loop, for example.
                
                StartingValues_xEOS     =   PhaseData{id}.StartingValues_xEOS;
                SolutionModels          =   StartingValues_xEOS.SolutionModels;
                InitialGuess_Num        =   zeros(size(SolutionModels));            % Change that if we want >1
                CompositionalVar        =   StartingValues_xEOS.CompositionalVar;
                n_phases                =   length(SolutionModels);                 % solution models
            else
                n_phases                =   0;
            end
        else
            n_phases = 0;
        end
    end
    fprintf(fid,'%i %f %f %f %f %f %f %f %f %f %f %f %f %f\n',n_phases, P,T,Gamma);
    
    if Use_xEOS & n_phases>0
        % Add compositional variable guesses
        for i=1:n_phases
            fprintf(fid,'%s ',SolutionModels{i}  );
            fprintf(fid,'%i ',InitialGuess_Num(i));     % we cane expand this later to take
            for iVar=1:length(CompositionalVar{i})
                fprintf(fid,'%f ',CompositionalVar{i}(iVar));
            end
            fprintf(fid,' \n');
        end
        
    end
    
end
fclose(fid);



%--------------------------------------------------------------------------
% Read out
function [df, prop, sf] = Read_OUTPUT()

curdir = pwd;
cd('./output');

num = 0;
fid = fopen('_thermocalc_style_output.txt');
while ~feof(fid)
    line = fgetl(fid);
    if strcmp(line,'============================================================')
        num = num + 1
    end
    
    if strcmp(line, 'site fractions')
        line     = fgetl(fid);   % name of ss
        line        = fgetl(fid);
        if ~strcmp(line,'oxide compositions')
            sf(num,:) 	= strread(line,'%f');   % actual sf
        else
            sf(num,:) = NaN;
        end
    end
    
    if strcmp(line, 'Driving force')
        line      = fgetl(fid);
        if ~isempty(line)
            [~,df(num)]   = strread(line,'%s %f');
        else
            df(num) = NaN;
        end
      
    end
    
    if length(line)>5
        if strcmp(line(1:5), 'ideal')
            % read eos of SS
            i = 1;
            while ~isempty(line)
                line = fgetl(fid);
                [a] 	= strread(line,'%f');
                if ~isempty(a)
                    prop(num,i)  =   a(4);
                    i        =   i+1;
                end
            end
            
        end
    end
    
    
end
fclose(fid);

cd(curdir)

