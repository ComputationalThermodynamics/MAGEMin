function em_names = Get_endmembers(UIAxes,PseudoSectionData, phase_in)
    em_names = {};
    if isfield(PseudoSectionData,'PhaseData')
        found   = 0
        n       = 1
        np      = length(PseudoSectionData.PhaseData);
        while found == 0 & n <= np
            for j=1:length(PseudoSectionData.PhaseData{1, n}.StableFractions)
                if PseudoSectionData.PhaseData{1, n}.StableSolutions{j, 1} == string(phase_in)
                    em_names    = PseudoSectionData.PhaseData{1, n}.EMlist{j};
                    found       = 1;
                end
            end
            n = n + 1;
        end
    else
        disp("!!! perform a simulation before trying to plot isocontours");
    end
