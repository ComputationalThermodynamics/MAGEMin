function [h,cap] = Plot_IsoContour(UIAxes,PseudoSectionData, data)
    DiagramType  =  PseudoSectionData.Computation.DiagramType;

    if isfield(PseudoSectionData,'PhaseData')
        if string(data.mode) == "mode"
            np = length(PseudoSectionData.PhaseData);
            axisY = zeros(np,1);
            axisX = zeros(np,1);
            X = zeros(np,1);


            if DiagramType == "PT";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.P;
                    axisX(i) = PseudoSectionData.PhaseData{1,i}.T;
                end
            elseif DiagramType == "PX";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.P;
                    axisX(i) = PseudoSectionData.XY_vec(i,1);
                end
            elseif DiagramType == "TX";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.T;
                    axisX(i) = PseudoSectionData.XY_vec(i,1);
                end
            end


            ngp      = str2num(data.res);
            axisXstep    = (max(axisX) - min(axisX))/(ngp-1);
            axisYstep    = (max(axisY) - min(axisY))/(ngp-1);

            axisXx       = min(axisX):axisXstep:max(axisX);
            axisYx   = min(axisY):axisYstep:max(axisY);
            [Xm,Ym]  = meshgrid(axisXx,axisYx);

            phase_in = string(data.phase);


            for i=1:np
                for j=1:length(PseudoSectionData.PhaseData{1,i}.StableFractions)
                    if PseudoSectionData.PhaseData{1, i}.StableSolutions{j, 1} == phase_in
                        X(i) = PseudoSectionData.PhaseData{1,i}.StableFractions(j);
                    end
                end
            end

            Xx   = griddata(axisX,axisY,X,Xm,Ym);
            minX = min(Xx(:));
            maxX = max(Xx(:));
            v    = str2num(data.min) + 1e-10:str2num(data.stp):str2num(data.max) + 1e-10;

            [c,h]  = contour(UIAxes,axisXx,axisYx,Xx,v, data.style,'ShowText','on');
            clabel(c,h,'FontSize',str2num(data.FtSize),'Color',data.lblC);
            cap = phase_in;
        elseif string(data.mode) == "em frac"
            np = length(PseudoSectionData.PhaseData);
            axisY = zeros(np,1);
            axisX = zeros(np,1);
            X = zeros(np,1);

            if DiagramType == "PT";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.P;
                    axisX(i) = PseudoSectionData.PhaseData{1,i}.T;
                end
            elseif DiagramType == "PX";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.P;
                    axisX(i) = PseudoSectionData.XY_vec(i,1);
                end
            elseif DiagramType == "TX";
                for i=1:np
                    axisY(i) = PseudoSectionData.PhaseData{1,i}.T;
                    axisX(i) = PseudoSectionData.XY_vec(i,1);
                end
            end
            ngp      = str2num(data.res);
            axisXstep    = (max(axisX) - min(axisX))/(ngp-1);
            axisYstep    = (max(axisY) - min(axisY))/(ngp-1);

            axisXx       = min(axisX):axisXstep:max(axisX);
            axisYx       = min(axisY):axisYstep:max(axisY);
            [Xm,Ym]  = meshgrid(axisXx,axisYx);

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

            Xx   = griddata(axisX,axisY,X,Xm,Ym);
            minX = min(Xx(:));
            maxX = max(Xx(:));
            v    = str2num(data.min) + 1e-10:str2num(data.stp):str2num(data.max) + 1e-10;

            [c,h]  = contour(UIAxes,axisXx,axisYx,Xx,v, data.style,'ShowText','on');
            clabel(c,h,'FontSize',str2num(data.FtSize),'Color',data.lblC);
            cap = phase_in+'.'+em_in;
        end 
        cap = string(cap);
    else
        disp("!!! perform a simulation before trying to plot isocontours");
    end
