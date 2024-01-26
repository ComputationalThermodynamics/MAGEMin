function [Local] = Read_LocalOptimizations_MAGEMin(dir)

curdir = pwd;
cd(dir)

% read solution model & # of variables
fid = fopen('__LOCAL_MINIMA.txt','r');
line = fgetl(fid); % comment
line = fgetl(fid); % comment

line = fgetl(fid); % contains data
str = strsplit(line);

Local.SolutionModel = str{1};
numVar      = str2num(str{2});
numPoints   = str2num(str{3});

% Extract Gamma
for i=1:11
    Local.Gamma(i) = str2num(str{3+i});
end
fclose(fid);


% Requires >2019b
% Data = readmatrix('__LOCAL_MINIMA.txt');
% Data(1,:)=[];


data_table          =   readtable('__LOCAL_MINIMA.txt','headerlines',4);        % works better on older matlab versions
if iscell(table2array(data_table(1,end)))
    data_table          =   data_table(:,1:end-1);
end
Data                =   table2array(data_table);


n                   =   1;
Local.Num           =   Data(:,n);
n                   =   n+1;

Local.Prop_start    =   Data(:,n:n+numVar);
n                   =   n+ numVar+1;

Local.xEOS_start    =   Data(:,n:n+numVar-1);
n                   =   n+ numVar;

Local.xEOS_end      =   Data(:,n:n+numVar-1);
n                   =   n+ numVar;

Local.Prop_end      =   Data(:,n:n+numVar);
n                   =   n+ numVar+1;

Local.df            =   Data(:,n);
n                   =   n+1;

Local.Status        =   Data(:,n);
n                   =   n+1;

Local.SF_ok         =   Data(:,n);

cd(curdir)




