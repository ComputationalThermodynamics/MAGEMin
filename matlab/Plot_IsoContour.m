function [h,cap] = Plot_IsoContour(UIAxes,PseudoSectionData, data)

    if isfield(PseudoSectionData,'PhaseData')
        if string(data.mode) == "mode"
            np = length(PseudoSectionData.PhaseData);
            P = zeros(np,1);
            T = zeros(np,1);
            X = zeros(np,1);

            for i=1:np
                P(i) = PseudoSectionData.PhaseData{1,i}.P;
                T(i) = PseudoSectionData.PhaseData{1,i}.T;
            end
            ngp      = str2num(data.res);
            Tstep    = (max(T) - min(T))/(ngp-1);
            Pstep    = (max(P) - min(P))/(ngp-1);

            Tx       = min(T):Tstep:max(T);
            Px       = min(P):Pstep:max(P);
            [Xm,Ym]  = meshgrid(Tx,Px);

            phase_in = string(data.phase);


            for i=1:np
                for j=1:length(PseudoSectionData.PhaseData{1,i}.StableFractions)
                    if PseudoSectionData.PhaseData{1, i}.StableSolutions{j, 1} == phase_in
                        X(i) = PseudoSectionData.PhaseData{1,i}.StableFractions(j);
                    end
                end
            end

            Xx   = griddata(T,P,X,Xm,Ym);
            minX = min(Xx(:));
            maxX = max(Xx(:));
            v    = str2num(data.min) + 1e-10:str2num(data.stp):str2num(data.max) + 1e-10;

            [c,h]  = contour(UIAxes,Tx,Px,Xx,v, data.style,'ShowText','on');
            clabel(c,h,'FontSize',str2num(data.FtSize),'Color',data.lblC);
            cap = phase_in;
        elseif string(data.mode) == "em frac"
            np = length(PseudoSectionData.PhaseData);
            P = zeros(np,1);
            T = zeros(np,1);
            X = zeros(np,1);

            for i=1:np
                P(i) = PseudoSectionData.PhaseData{1,i}.P;
                T(i) = PseudoSectionData.PhaseData{1,i}.T;
            end
            ngp      = str2num(data.res);
            Tstep    = (max(T) - min(T))/(ngp-1);
            Pstep    = (max(P) - min(P))/(ngp-1);

            Tx       = min(T):Tstep:max(T);
            Px       = min(P):Pstep:max(P);
            [Xm,Ym]  = meshgrid(Tx,Px);

            phase_in = string(data.phase);
            em_in    = string(data.em);

            for i=1:np
                for j=1:length(PseudoSectionData.PhaseData{1,i}.StableFractions)
                    if PseudoSectionData.PhaseData{1, i}.StableSolutions{j, 1} == phase_in

                        for k = 1:length(PseudoSectionData.PhaseData{1, i}.EMlist{j})
                            if PseudoSectionData.PhaseData{1, i}.EMlist{j}{k} == em_in
                                X(i) = PseudoSectionData.PhaseData{1, i}.EMFractions{j}(k);
                            end
                        end

                    end
                end
            end

            Xx   = griddata(T,P,X,Xm,Ym);
            minX = min(Xx(:));
            maxX = max(Xx(:));
            v    = str2num(data.min) + 1e-10:str2num(data.stp):str2num(data.max) + 1e-10;

            [c,h]  = contour(UIAxes,Tx,Px,Xx,v, data.style,'ShowText','on');
            clabel(c,h,'FontSize',str2num(data.FtSize),'Color',data.lblC);
            cap = phase_in+'.'+em_in;
        end 
        cap = string(cap);
    else
        disp("!!! perform a simulation before trying to plot isocontours");
    end
