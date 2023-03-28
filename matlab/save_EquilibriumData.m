function save_EquilibriumData(PseudoSectionData)

    solver    	    =   PseudoSectionData.Computation.solver;
    db    	  	    =   PseudoSectionData.Computation.db;
    k               =   PseudoSectionData.activePoint;
    P               =   PseudoSectionData.PhaseData{k}.P;
    T               =   PseudoSectionData.PhaseData{k}.T;
    sys_in          =   PseudoSectionData.Chemistry.sys_in;               % predefined chemical composition
    Test            =   PseudoSectionData.Chemistry.Predefined;
    
    if isnan(Test)
        OxProp      =  PseudoSectionData.Chemistry.OxProp;    % we do not employ a predefined test, but specify mol proportions insteadOxProp
    else
        OxProp      =   [];
    end

    exe             = MAGEMin_exe(PseudoSectionData.Computation);

    command         = [exe,' --out_matlab=1 --solver=',num2str(solver),' --Verb=-1 --sys_in=',sys_in,' --db=',db, ' --Pres=',num2str(P), ' --Temp=',num2str(T)];


    if ~isnan(Test)
        % employ a prededined test
        command = [command, ' --test=',num2str(Test)];
    else
        % employ specified chemistry
        command = [command, ' --Bulk='];
        for iTable=1:size(OxProp,1)
            command = [command,num2str(table2array(OxProp(iTable,2)),'%.8f'),','];
        end
    end

    disp(command)
    if PseudoSectionData.Computation.Julia_MAGEMin_binary==true
        command = add_dynamic_libs(command, Computation);
    end
    if isunix
        system('killall MAGEMin 2>&1');
    end
    system(command);

    dataOut = string(fileread('./output/_matlab_output.txt'));

    filter = {'*.txt'};
    [file, path] = uiputfile(filter);

    fid = fopen(file,'wt');
    fprintf(fid, dataOut);
    fclose(fid);