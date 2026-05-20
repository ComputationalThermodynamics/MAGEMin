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
    TE_database

    Mutable structure holding a fixed trace element partitioning coefficient (KD) database.

    Fields
    ------
    conditions : String
        Description of the P-T or compositional conditions for which the KDs are valid.
    infos : String
        Citation or provenance of the database.
    n_element : Int64
        Number of trace elements.
    n_phase : Int64
        Number of mineral phases.
    element_name : Vector{String}
        Names of the trace elements.
    phase_name : Vector{String}
        Names of the mineral phases (in MAGEMin-compatible notation).
    Kds : Matrix{Float64}
        Matrix of partition coefficients (phases × elements).
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
    mineral_name_convertor(phase_name)

    Convert mineral phase names from external TE database conventions to MAGEMin-compatible short names.

    Parameters
    ----------
    phase_name : Vector{String}
        Phase names using external conventions (e.g., "Bt", "Grt", "Kfs", "Pl").

    Returns
    -------
    ph : Vector{String}
        Phase names in MAGEMin notation (e.g., "bi", "g", "afs", "pl").
"""
function mineral_name_convertor(    phase_name      :: Vector{String}   )

    ph = copy(phase_name)
    ph = replace.(ph,"All"  => "all");      ph = replace.(ph,"Ap"   => "ap");       ph = replace.(ph,"Ttn"  => "ttn");      ph = replace.(ph,"Amp"  => "amp");
    ph = replace.(ph,"Bt"   => "bi");       ph = replace.(ph,"Crd"  => "cd");       ph = replace.(ph,"Kfs"  => "afs");      ph = replace.(ph,"Pl"   => "pl");
    ph = replace.(ph,"Qtz"  => "q");        ph = replace.(ph,"Rt"   => "ru");       ph = replace.(ph,"Grt"  => "g");        ph = replace.(ph,"Ep"   => "ep");
    ph = replace.(ph,"Ol"   => "ol");       ph = replace.(ph,"Opx"  => "opx");      ph = replace.(ph,"Mica" => "mu");       ph = replace.(ph,"Mt"   => "mt");
    ph = replace.(ph,"zrc"  => "zrc");      ph = replace.(ph,"Cpx"  => "cpx");      ph = replace.(ph,"Spl"  => "sp");

    return ph 
end


"""
    mineral_classification(out, dtb)

    Classify the stable phases from a MAGEMin minimization result into mineralogical names compatible with the trace element partitioning coefficient database.

    Solution phases that straddle a solvus (feldspar, spinel, ilmenite) are disambiguated using their compositional variables. Clinopyroxene variants ("dio", "aug") are unified as "cpx"; ilmenite variants ("ilm", "ilmm") are unified as "FeTiOx".

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output structure.
    dtb : String
        Database identifier (e.g., "ig", "igad", "mp", "mpe", "mb", "ume", "mbe").

    Returns
    -------
    ph : Vector{String}
        Classified phase names compatible with the TE partitioning database.
    ph_wt : Vector{Float64}
        Weight fractions corresponding to each phase.
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
                    ph[i] = "smt"    # magnetite
                else
                    ph[i] = "spl"    # spinel
                end
            end
        end
        if ss == "sp"
            if out.SS_vec[i].compVariables[2] + out.SS_vec[i].compVariables[3] - 0.5 > 0
                ph[i] = "smt"        # magnetite
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

