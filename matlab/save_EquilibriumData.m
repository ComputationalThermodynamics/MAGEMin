function save_EquilibriumData(PseudoSectionData)

    DiagramType     =   PseudoSectionData.Computation.DiagramType;

    solver    	    =   PseudoSectionData.Computation.solver;
    db    	  	    =   PseudoSectionData.Computation.db;
    k               =   PseudoSectionData.activePoint;
    P               =   PseudoSectionData.PhaseData{k}.P;
    T               =   PseudoSectionData.PhaseData{k}.T;

    sys_in          =   PseudoSectionData.Chemistry.sys_in;               % predefined chemical composition
    Test            =   PseudoSectionData.Chemistry.Predefined;
    limitCaOpx      =   PseudoSectionData.Computation.limitCaOpx;
    CaOpxLim        =   PseudoSectionData.Computation.CaOpxLim;
    buffer          =   PseudoSectionData.Computation.buffer;
    buffer_n        =   PseudoSectionData.Computation.buffer_n;
    mbCpx           =   PseudoSectionData.Computation.mbCpx;

    if isnan(Test)
        OxProp      =  PseudoSectionData.Chemistry.OxProp;    % we do not employ a predefined test, but specify mol proportions insteadOxProp
    else
        OxProp      =   [];
    end

    exe             = MAGEMin_exe(PseudoSectionData.Computation);


    if DiagramType == 'PT'
        command         = [exe,' --out_matlab=2 --solver=',num2str(solver),' --limitCaOpx=',num2str(limitCaOpx), ' --CaOpxLim=',num2str(CaOpxLim),' --buffer=',num2str(buffer),' --buffer_n=',num2str(buffer_n),' --Verb=-1 --sys_in=',sys_in,' --db=',db, ' --Pres=',num2str(P), ' --Temp=',num2str(T)];

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
    else
        k               = PseudoSectionData.activePoint
        X               = PseudoSectionData.XY_vec(k,1);
 
        Bulk            =  table2array(PseudoSectionData.Chemistry.OxProp(:,2))*(1.-PseudoSectionData.XY_vec(k,1)) + table2array(PseudoSectionData.Chemistry.OxProp(:,3))*(PseudoSectionData.XY_vec(k,1));
        command         = [exe,' --out_matlab=2 --solver=',num2str(solver),' --limitCaOpx=',num2str(limitCaOpx), ' --CaOpxLim=',num2str(CaOpxLim),' --buffer=',num2str(buffer),' --buffer_n=',num2str(buffer_n),' --Verb=-1 --sys_in=',sys_in,' --db=',db, ' --Pres=',num2str(P), ' --Temp=',num2str(T)];

        % employ specified chemistry
        command = [command, ' --Bulk='];
        command = [command,sprintf('%.8f,' , Bulk)];
    end


    if strcmp(string(db),"mb") == 1
        switch mbCpx
        case 'Omph'
            command = [command, ' --mbCpx=0']
        case 'Aug'
            command = [command, ' --mbCpx=1']
        otherwise
            disp('wrong cpx, something is fishy')
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
    dataOut = strrep(dataOut,'%','%%');
    filter  = {'*.txt'};
    [file, path] = uiputfile(filter);

    fid = fopen(file,'a+');
    fprintf(fid, dataOut);
    fclose(fid);