% Trace Element Predictive Modelling has been developped by Julien Cornet 
% 2018, MsC thesis, Durham University, DOI: 10.13140/RG.2.2.25669.14562.
% Julien's work accounts for olivine, amphibole, plagioclase, garnet, clinopyroxene and 
% orthopyroxene and can be applied to mafic to ultramafic systems at suprasolidus conditions.

function PseudoSectionData = TEpredictingModel(PseudoSectionData)

    [constant,ri_CN_8,ri_8,fpe,ri_CN_6,ri_6,Molar_mass,flag,order_flag,TE_name1,TE_name,clr]=constant_input;

    % working oxide list for TEPM
    wOx         = {'SiO2' 'TiO2' 'Al2O3' 'FeO' 'MnO' 'MgO' 'CaO' 'Na2O' 'K2O' 'H2O'};
    Mox         = [ 60.08  79.88 101.96 71.85 70.94 40.30 56.08 61.98 94.2 18.015];

    TEbulk      = table2array(PseudoSectionData.TraceElement.OxProp(:,2));
    TEname      = table2cell(PseudoSectionData.TraceElement.OxProp(:,1));
    n_TE        = length(TEbulk);

    % allocate memory
    np_tot      = length(PseudoSectionData.PhaseData);
    idOx        = zeros(length(wOx),1);
    fig         = 0;
    p           = 0; %dummy variable
    % load system components of the calculation by directly reading the first solution phase
    for i=1:np_tot
        if isempty(PseudoSectionData.PhaseData{1, 1}.SSoxide{1, i} ) == 0
            SSoxide = PseudoSectionData.PhaseData{1, 2}.SSoxide{1,i};
            break;
        end
    end

    % map oxides
    for i=1:length(wOx)
        if string(wOx(i)) == 'MnO'
            idOx(i) = 1;
        else
            idOx(i) = find(strcmp(SSoxide,string(wOx(i))));
        end
    end
    
    % create structure that will hold trace elements information of the system
    SSname = {'m','s','amp','pl','cpx','g','opx','ol'};
    for i=1:length(TEname)
        for j=1:length(SSname)
            str = strcat(string(SSname(j)),'.',string(TEname(i)),' = zeros(np_tot,1);');
            eval(str);
        end
    end

    % solution phase composition
    for i=1:np_tot
        c_mx        = zeros(6,28);
        Di          = zeros(6,28);
        F           = zeros(7,1);
        SScomp_wt   = PseudoSectionData.PhaseData{1, i}.SScomp_wt;
        phase       = PseudoSectionData.PhaseData{1, i}.StableSolutions;
        fraction    = cell2mat(PseudoSectionData.PhaseData{1, i}.PHfractions_wt);

        T           = PseudoSectionData.PhaseData{1, i}.T;
        P           = PseudoSectionData.PhaseData{1, i}.P;

        if isempty(find(strcmp(phase,'liq'))) == 0
            % here the melt composition is saved for partitioning trace element with the solution phases
            id_melt      = find(strcmp(phase,'liq'));
            comp(1,:)    = SScomp_wt{1, id_melt}(idOx) .*100; comp(1,5) = 0.0;
            mc(1,:)      = comp(1,:)./Mox;
            F(1)         = fraction(id_melt);
            for ph=1:length(phase)
                % here composition of the solution phase is saved next to the melt composition
                comp(2,:)    = SScomp_wt{1, ph}(idOx) .*100; comp(5) = 0.0;
                mc(2,:)      = comp(2,:)./Mox;
                Fr            = fraction(ph);
                ss           = string(phase(ph));
                
                switch ss
                    case {"amp","anth"}
                        Di(1,:)  = TE_amph(P,T,mc,comp,ri_CN_8,ri_CN_6,fpe,ri_8,ri_6,order_flag,TE_name,fig,TE_name1,p,clr);
                        c_mx(1,:)= Fr.*Di(1,:);
                        F(2)     = Fr;
                    case {"pl4tr","pl4T","fsp"}
                        Di(2,:)  = TE_pl(constant,fpe,T,P,ri_CN_8,comp,mc,ri_8,order_flag,TE_name,fig,p,clr);
                        c_mx(2,:)= Fr.*Di(2,:);
                        F(3)     = Fr;
                    case {"cpx","dio","aug"}
                        Di(3,:)  = TE_cpx(constant,T,P,ri_CN_8,ri_CN_6,comp,mc,fpe,ri_8,ri_6,order_flag,TE_name,fig,p,clr);
                        c_mx(3,:)= Fr.*Di(3,:);
                        F(4)     = Fr;
                    case "g"
                        Di(4,:)  = TE_gt(constant,T,P,ri_CN_8,ri_CN_6,mc,fpe,comp,ri_8,ri_6,order_flag,TE_name,fig,p,clr);
                        c_mx(4,:)= Fr.*Di(4,:);
                        F(5)     = Fr;
                    case "opx"
                        Di(5,:)  = TE_opx(constant,T,P,ri_CN_8,ri_CN_6,comp,mc,fpe,ri_8,ri_6,order_flag,TE_name,fig,p,clr);
                        c_mx(5,:)= Fr.*Di(5,:);
                        F(6)     = Fr;
                    case "ol"
                        Di(6,:)  = TE_ol(T,P,ri_CN_6,comp,mc,fpe,ri_6,fig,clr);
                        c_mx(6,:)= Fr.*Di(6,:);
                        F(7)     = Fr;
                end
            end
            Dwr     = sum(c_mx,1);

            % Modal Batch Melting
            c_m    = TEbulk' ./(F(1)+Dwr.*(1.0-F(1)));               % Melt composition
            c_s    = Dwr.*c_m;                                    % Solid composition
            
            for k=1:length(TEname)
                m.(TEname{k})(i) = c_m(k);
                s.(TEname{k})(i) = c_s(k);
            end
            % fns = fieldnames(amp);

            % % Modal Fractional Melting MFM   
            % c_m= (MIXB./Dwr).*(1-(F(1)).^((1./Dwr)-1));            % Melt composition
            % c_s=Dwr.*c_m;                                         % Solid composition
            % if isnan(c_s)==1;c_s=c0;end; c_m(isinf(c_m))=0;

            c_cpx(1,:)     = zeros(1,28);
            c_gt(1,:)      = zeros(1,28);
            c_amph(1,:)    = zeros(1,28);
            c_opx(1,:)     = zeros(1,28);
            c_pl(1,:)      = zeros(1,28);
            c_ol(1,:)      = zeros(1,28);

            if F(4)~= 0
                c_cpx= c_s./(F(4)+(Di(4,:).*F(5)+Di(1,:).*F(2)+Di(5,:).*F(6)+Di(2,:).*F(3)+Di(6,:).*F(7)+F(1))./(Di(3,:)));
                for k=1:length(TEname)
                    cpx.(TEname{k})(i) = c_cpx(k);
                end
            end
            if F(5)~=0
                c_gt=c_s./(F(5)+(Di(3,:).*F(4)+Di(1,:).*F(2)+Di(5,:).*F(6)+Di(2,:).*F(3)+Di(6,:).*F(7)+F(1))./(Di(4,:))); 
                for k=1:length(TEname)
                    g.(TEname{k})(i) = c_gt(k);
                end
            end
            if F(2)~=0
                c_amph=c_s./(F(2)+(Di(4,:).*F(5)+Di(3,:).*F(4)+Di(5,:).*F(6)+Di(2,:).*F(3)+Di(6,:).*F(7)+F(1))./(Di(1,:)));
                for k=1:length(TEname)
                    amp.(TEname{k})(i) = c_amph(k);
                end
            end
            if F(6)~=0
                c_opx=c_s./(F(6)+(Di(4,:).*F(5)+Di(1,:).*F(2)+Di(3,:).*F(4)+Di(2,:).*F(3)+Di(6,:).*F(7)+F(1))./(Di(5,:))); 
                for k=1:length(TEname)
                    opx.(TEname{k})(i) = c_opx(k);
                end
            end
            if F(3)~=0
                c_pl=c_s./(F(3)+(Di(4,:).*F(5)+Di(1,:).*F(2)+Di(5,:).*F(6)+Di(3,:).*F(4)+Di(6,:).*F(7)+F(1))./(Di(2,:))); 
                for k=1:length(TEname)
                    pl.(TEname{k})(i) = c_pl(k);
                end
            end
            if F(7)~=0
                c_ol=c_s./(F(7)+(Di(4,:).*F(5)+Di(1,:).*F(2)+Di(5,:).*F(6)+Di(2,:).*F(3)+Di(3,:).*F(4)+F(1))./(Di(6,:))); 
                for k=1:length(TEname)
                    ol.(TEname{k})(i) = c_ol(k);
                end
            end
        else
            for k=1:length(TEname)
                m.(TEname{k})(i)   = NaN;
                s.(TEname{k})(i)   = NaN;
                amp.(TEname{k})(i)  = NaN;
                g.(TEname{k})(i)   = NaN;
                pl.(TEname{k})(i)  = NaN;
                ol.(TEname{k})(i)  = NaN;
                opx.(TEname{k})(i) = NaN;
                cpx.(TEname{k})(i) = NaN;
            end

        end
    end

    PseudoSectionData.TraceElement.s   = s;
    PseudoSectionData.TraceElement.m   = m;
    PseudoSectionData.TraceElement.amp  = amp;
    PseudoSectionData.TraceElement.g   = g;
    PseudoSectionData.TraceElement.pl  = pl;
    PseudoSectionData.TraceElement.ol  = ol;
    PseudoSectionData.TraceElement.opx = opx;
    PseudoSectionData.TraceElement.cpx = cpx;
end