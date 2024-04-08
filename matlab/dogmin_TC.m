function [Data_best] = dogmin_TC(Data, File_dir)
% This takes the MAGEMin solution and runs TC multiple times to check all
% possible permutations of the MAGEMin solution to see if TC finds a more
% stable option
%


Sol_nr  = 1:length(Data.StableSolutions);
num     = 1;

time_start = cputime;
% Perform calculations with all models passed:
[Data_dogmin{num}]          =   MAGEMin_to_TC(File_dir,'tc-testMAGEMin.txt',Data, logical(0));

if ~Data_dogmin{num}.success
    Data_dogmin{num}.StableSolutions =  Data.StableSolutions;
    Data_dogmin{num}.Gibbs           =  NaN;
end

Gibbs(num)  =   Data_dogmin{num}.Gibbs;    % to sort them accordingly

num     = num+1;

% Perform a whole range of TC calculations with all kinds of permutations
for iPerm=2:length(Sol_nr)
    disp(' ')
    permutations = nchoosek(Sol_nr,iPerm);
    for j=1:size(permutations,1)
        Data_in                     =	Data;
        Data_in.StableSolutions     =	Data.StableSolutions(permutations(j,:));
        Data_in.StableFractions     =   Data.StableFractions(permutations(j,:));
        Data_in.CompositionalVar    =   Data.CompositionalVar(permutations(j,:));
        
        TC_Output                   =   logical(1);
        [Data_dogmin{num}]          =   MAGEMin_to_TC(File_dir,'tc-testMAGEMin.txt',Data_in, TC_Output);
        
        if ~Data_dogmin{num}.success
            Data_dogmin{num}.StableSolutions =  Data_in.StableSolutions;
            Data_dogmin{num}.Gibbs           =  NaN;
        end
        
        Gibbs(num)  =   Data_dogmin{num}.Gibbs;    % to sort them accordingly
        
        num         =   num+1;
    end
    
end

% Print results
id_best = DisplayResults(Gibbs, Data, Data_dogmin);
if ~isempty(id_best)
    Data_best = Data_dogmin{id_best};
else
    Data_best.success=false;
end


time_end = cputime;
disp(['dogmin calculation took ',num2str(time_end-time_start),'s'])



%--------------------------------------------------------------------------
function [id_best] = DisplayResults(Gibbs, Data, Data_dogmin)


% Display them in order of importance
num                 =   1:length(Gibbs);
[Gibbs,id]          =   sort(Gibbs);
id(isnan(Gibbs))    =   [];

if length(id)>0
    id_best             =   id(1);
    
    disp(' ')
    disp(['  ---------------------------------------'])
    disp(['  Results of TC optimizations:'])
    disp(['     P = ',num2str(Data.P),' kbar'])
    disp(['     T = ',num2str(Data.T),' Celsius'])
    disp('   ')
    str = '  mode : ';
    for i=1:length(Data.StableSolutions)
        str = [str, sprintf('%-8s ',Data.StableSolutions{i})];
    end
    str = [str, sprintf('   %-8s ','G')  ];
    str = [str, sprintf('  %-8s ','del')];
    disp(str)
    
    del = Gibbs-Gibbs(1);
    
    for i=1:length(id)
        Dat  = Data_dogmin{id(i)};
        
        str = sprintf('   # %-4i',id(i));
        for j=1:length(Data.StableSolutions)
            
            ind = find(strcmp(Data.StableSolutions(j),Dat.StableSolutions));
            if length(ind)>0
                frac = Dat.StableFractions(ind);
            else
                frac = 0;
            end
            str = [str, sprintf('%8.6f ',frac)];
        end
        str = [str, sprintf('  %-10.6f ',Gibbs(i))];
        str = [str, sprintf('%-10.6f ',del(i))];
        disp(str)
        
    end
else
    id_best=[];
end



