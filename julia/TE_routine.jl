#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
# Trace element partitioning prediction routine
# Note that only metapelite, metabasite and igneous database can be used for trace element prediction
# NR 11/04/2023

using MAGEMin_C  

"""
    Holds the partitioning coefficient database
"""
mutable struct TE_database
    conditions      :: String
    infos           :: String
    n_element       :: Int64
    n_phase         :: Int64
    element_name    :: Vector{String}
    phase_name      :: Vector{String}
    Kds             :: Matrix{Float64}
end

"""
    This function converts mineral names from the TE database into names that can be compared with MAGEMin mineral classification
"""
function mineral_name_convertor(    phase_name      :: Vector{String}   )

    ph = copy(phase_name)
    ph = replace.(ph,"All"  => "all");      ph = replace.(ph,"Ap"   => "ap");       ph = replace.(ph,"Ttn"  => "ttn");      ph = replace.(ph,"Amp"  => "amp");
    ph = replace.(ph,"Bt"   => "bi");       ph = replace.(ph,"Crd"  => "cd");       ph = replace.(ph,"Kfs"  => "afs");      ph = replace.(ph,"Pl"   => "pl");
    ph = replace.(ph,"Qtz"  => "q");        ph = replace.(ph,"Rt"   => "ru");       ph = replace.(ph,"Grt"  => "g");        ph = replace.(ph,"Ep"   => "ep");
    ph = replace.(ph,"Ol"   => "ol");       ph = replace.(ph,"Opx"  => "opx");      ph = replace.(ph,"Mica" => "mu");       ph = replace.(ph,"Mt"   => "mt");
    ph = replace.(ph,"Zrn"  => "zrn");      ph = replace.(ph,"Cpx"  => "cpx");      ph = replace.(ph,"Spl"  => "sp");

    return ph 
end


"""
    Classify the mineral output from MAGEMin to be able to be compared with partitioning coefficient database
"""
function mineral_classification(    out             :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    dtb             :: String  )

    ph      = Array{String}(undef, out.n_SS + out.n_PP) 
    ph_wt   = Array{Float64}(undef, out.n_SS + out.n_PP) 

    # add solution phase and classify some solution phases (spl, fsp, ilm)                             
    for i = 1:out.n_SS                             
        ss      = out.ph[i]
        ph_wt[i]= out.ph_frac_wt[i]
        ph[i]   = ss
        if ss == "fsp"
            if out.SS_vec[i].compVariables[2] - 0.5 > 0
                ph[i] = "afs"
            else
                ph[i] = "pl"
            end
        end
        if ss == "spl"
            if out.SS_vec[i].compVariables[3] - 0.5 > 0
                ph[i] = "cm"        # chromite
            else
                if out.SS_vec[i].compVariables[2] - 0.5 > 0
                    ph[i] = "mt"    # magnetite
                else
                    ph[i] = "sp"    # spinel
                end
            end
        end
        if ss == "sp"
            if out.SS_vec[i].compVariables[2] + out.SS_vec[i].compVariables[3] - 0.5 > 0
                ph[i] = "mt"        # chromite
            else
                if (1 - out.SS_vec[i].compVariables[1])*(1 + out.SS_vec[i].compVariables[3]) - 0.5 > 0
                    ph[i] = "sp"    # spinel
                else
                    if out.SS_vec[i].compVariables[3] -0.5 > 0
                        ph[i] = "FeTiOx"  # uvospinel
                    else
                        ph[i] = "sp" # hercynite
                    end
                end
            end
        end
        if ss == "dio" || ss == "aug"
            ph[i] = "cpx"
        end
        if ss == "ilm" || ss == "ilmm"
            ph[i] = "FeTiOx"
        end
        # add pure phases
        for i=1:out.n_PP
            ph[i+out.n_SS]      = out.ph[i+out.n_SS]
            ph_wt[i+out.n_SS]   = out.ph_frac_wt[i+out.n_SS]
        end
        
    end

    return ph, ph_wt
end

