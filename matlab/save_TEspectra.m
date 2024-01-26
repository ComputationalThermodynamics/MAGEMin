function save_TEspectra(PseudoSectionData,full)
    k               =   PseudoSectionData.activePointTE;
    P               =   PseudoSectionData.PhaseData{k}.P;
    T               =   PseudoSectionData.PhaseData{k}.T;
    StablePhases    =   PseudoSectionData.PhaseData{k}.StableSolutions;
    fraction        = cell2mat(PseudoSectionData.PhaseData{1, k}.PHfractions_wt);

    TEname          = table2cell(PseudoSectionData.TraceElement.OxProp(:,1));
    TEbulk          = table2array(PseudoSectionData.TraceElement.OxProp(:,2));
    TEout           = 'Trace Element Predictive Modelling; J. Cornet, 2018, MsC thesis, Durham University, DOI: 10.13140/RG.2.2.25669.14562\n============================================================\n';
    TEout           = strcat(TEout,' ',strjoin(StablePhases),' {',num2str(P),', ',num2str(T),'} kbar/°C\n');
    TEout           = strcat(TEout,'mod[wt]','\t',strjoin(string(fraction)),'\n\n');
    TEout           = strcat(TEout,'Phase','\t',strjoin(TEname),'\n');
    TEout           = strcat(TEout,'bulk','\t',strjoin(string(TEbulk)),'\n');

    for j=1:length(TEname)
        Data.TE.m.(TEname{j})   = PseudoSectionData.TraceElement.m.(TEname{j})(k);
        Data.TE.s.(TEname{j})   = PseudoSectionData.TraceElement.s.(TEname{j})(k);
        Data.TE.hb.(TEname{j})  = PseudoSectionData.TraceElement.hb.(TEname{j})(k);
        Data.TE.g.(TEname{j})   = PseudoSectionData.TraceElement.g.(TEname{j})(k);
        Data.TE.pl.(TEname{j})  = PseudoSectionData.TraceElement.pl.(TEname{j})(k);
        Data.TE.ol.(TEname{j})  = PseudoSectionData.TraceElement.ol.(TEname{j})(k);
        Data.TE.opx.(TEname{j}) = PseudoSectionData.TraceElement.opx.(TEname{j})(k);
        Data.TE.cpx.(TEname{j}) = PseudoSectionData.TraceElement.cpx.(TEname{j})(k);
    end
    
    phName = fieldnames(Data.TE);
    jx = 1;
    for ix=1:length(phName)
        if Data.TE.(phName{ix}).(TEname{1}) > 0
            x1 = zeros(length(TEname),1);
            for j=1:length(TEname)
                x1(j) = Data.TE.(phName{ix}).(TEname{j});
            end
            TEout = strcat(TEout,convertStringsToChars(phName{ix}),'\t',strjoin(string(x1)),'\n');
        end
    end
    TEout = strcat(TEout,'\n');

    dataOut = [];
    if full == 1
        DiagramType     =   PseudoSectionData.Computation.DiagramType;

        solver    	    =   PseudoSectionData.Computation.solver;
        db    	  	    =   PseudoSectionData.Computation.db;

        sys_in          =   PseudoSectionData.Chemistry.sys_in;               % predefined chemical composition
        Test            =   PseudoSectionData.Chemistry.Predefined;
        limitCaOpx      =   PseudoSectionData.Computation.limitCaOpx;
        CaOpxLim        =   PseudoSectionData.Computation.CaOpxLim;
        buffer             =   PseudoSectionData.Computation.buffer;
        buffer_n           =   PseudoSectionData.Computation.buffer_n;
        mbCpx           =   PseudoSectionData.Computation.mbCpx;

        if strcmp(db,'ig') == 1
            db              =   'Igneous';
        elseif strcmp(db,'igd') == 1
            db              =   'Igneous dry';
        elseif strcmp(db,'mp') == 1
            db              =   'Metapelite';
        elseif strcmp(db,'um') == 1
            db              =   'Ultramafic';
        elseif strcmp(db,'umj') == 1
            db              =   'Ultramafic Jun';
        end 
        Data                =   PseudoSectionData.PhaseData{k};
    
        if PseudoSectionData.Computation.DiagramType == 'PT'
            P               =   PseudoSectionData.PhaseData{k}.P;
            T               =   PseudoSectionData.PhaseData{k}.T;
        elseif PseudoSectionData.Computation.DiagramType == 'PX'
            P               =   PseudoSectionData.PhaseData{k}.P;
            T               =   PseudoSectionData.FixedTemperature;
            X               =   PseudoSectionData.XY_vec(k,1);
        elseif PseudoSectionData.Computation.DiagramType == 'TX'
            P               =   PseudoSectionData.FixedPressure;
            T               =   PseudoSectionData.PhaseData{k}.T;
            X               =   PseudoSectionData.XY_vec(k,1);
        end
        
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
    end

    dataOut = strcat(TEout,dataOut);

    filter  = {'*.txt'};
    [file, path] = uiputfile(filter);

    fid = fopen(file,'a+');
    fprintf(fid, dataOut);
    fclose(fid);
end