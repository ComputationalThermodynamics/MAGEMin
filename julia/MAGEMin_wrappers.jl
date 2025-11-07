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
# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia

import Base.show
using Base.Threads: @threads
using ProgressMeter
using DataFrames, Dates, CSV

const VecOrMat          = Union{Nothing, AbstractVector{Float64}, AbstractVector{<:AbstractVector{Float64}}}
const available_TC_ds   = [62,633,634,635,636]

export  anhydrous_renormalization, retrieve_solution_phase_information, remove_phases, get_ss_from_mineral, mineral_classification,
        init_MAGEMin, allocate_output,finalize_MAGEMin, point_wise_minimization, 
        get_all_stable_phases, convertBulk4MAGEMin, use_predefined_bulk_rock, define_bulk_rock, create_output,
        print_info, create_gmin_struct, pwm_init, pwm_run,
        point_wise_metastability,
        single_point_minimization, multi_point_minimization, AMR_minimization, MAGEMin_Data,
        MAGEMin_data2dataframe, MAGEMin_dataTE2dataframe, MAGEMin_data2dataframe_inlined,
        
        Initialize_MAGEMin, Finalize_MAGEMin

export wt2mol, mol2wt
export compute_melt_viscosity_G08

export TE_prediction, adjust_bulk_4_zircon, create_custom_KDs_database, get_TE_database, adjust_chemical_system
export zirconium_saturation, sulfur_saturation


export initialize_AMR, split_and_keep, AMR

export out_struct, out_TE_struct

include("name_solvus.jl")

"""
    Function to allocate memory for the output
"""
function allocate_output(n::Int64)
    return Vector{gmin_struct{Float64, Int64}}(undef, n)
end


"""
    bulk_dry = anhydrous_renormalization(   bulk    :: Vector{Float64},
                                            oxide   :: Vector{String})

    This function renormalizes the bulk rock composition to remove water (H2O) if present.

"""
function anhydrous_renormalization( bulk    :: Vector{Float64},
                                    oxide   :: Vector{String})

    if "H2O" in oxide
        H2O_index = findfirst(==("H2O"), oxide)
        bulk_dry = copy(bulk)
        bulk_dry[H2O_index] = 0.0
        bulk_dry ./= sum(bulk_dry)
    else
        # println("No water oxide in the system!")
        bulk_dry = bulk ./ sum(bulk)
    end

    return bulk_dry
end

"""
    structure that holds the result of the pointwise minimization
"""
struct gmin_struct{T,I}
    MAGEMin_ver :: String
    dataset     :: String
    database    :: String
    buffer      :: String
    buffer_n    :: T
    G_system    :: T             # G of system
    Gamma       :: Vector{T}        # Gamma
    P_kbar      :: T               # Pressure in kbar
    T_C         :: T                  # Temperature in Celsius
    X           :: Vector{T}
    M_sys       :: T

    # bulk rock composition:
    bulk        :: Vector{T}
    bulk_M      :: Vector{T}
    bulk_S      :: Vector{T}
    bulk_F      :: Vector{T}

    bulk_wt     :: Vector{T}
    bulk_M_wt   :: Vector{T}
    bulk_S_wt   :: Vector{T}
    bulk_F_wt   :: Vector{T}

    # Fractions:
    # Solid, melt, fluid fractions
    frac_M      :: T
    frac_S      :: T
    frac_F      :: T

    frac_M_wt   :: T
    frac_S_wt   :: T
    frac_F_wt   :: T

    frac_M_vol   :: T
    frac_S_vol   :: T
    frac_F_vol   :: T

    # Solid, melt, fluid densities
    alpha       :: Vector{T}
    V           :: T
    s_cp        :: Vector{T}
    rho         :: T
    rho_M       :: T
    rho_S       :: T
    rho_F       :: T
    eta_M       :: T

    # Oxygen fugacity
    fO2         :: T
    dQFM        :: T

    # Activities
    aH2O        :: T
    aSiO2       :: T
    aTiO2       :: T
    aAl2O3      :: T
    aMgO        :: T
    aFeO        :: T

    # Phase fractions and type:
    n_PP        :: Int64                 # number of pure phases
    n_SS        :: Int64                 # number of solid solutions
    n_mSS       :: Int64                 # number of solid solutions

    ph_frac     :: Vector{T}            # phase fractions
    ph_frac_wt  :: Vector{T}            # phase fractions
    ph_frac_1at :: Vector{T}            # phase fractions
    ph_frac_vol :: Vector{T}            # phase fractions
    ph_type     :: Vector{I}            # type of phase (SS or PP)
    ph_id       :: Vector{I}            # id of phase
    ph_id_db    :: Vector{I}            # id of phase
    ph          :: Vector{String}       # Name of phase
    sol_name    :: Vector{String}       # Name of phase

    SS_syms     :: Dict{Symbol, Int64}
    PP_syms     :: Dict{Symbol, Int64}

    SS_vec      :: Vector{LibMAGEMin.SS_data}
    mSS_vec     :: Vector{LibMAGEMin.mSS_data}
    PP_vec      :: Vector{LibMAGEMin.PP_data}

    oxides      :: Vector{String}
    elements    :: Vector{String}

    # Seismic velocity info
    Vp              :: T                # P-wave velocity
    Vs              :: T                # S-wave velocity
    Vp_S            :: T                # P-wave velocity of solid aggregate
    Vs_S            :: T                # S-wave velocity of solid aggregate
    bulkMod         :: T                # Elastic bulk modulus
    shearMod        :: T                # Elastic shear modulus
    bulkModulus_M   :: T                # Elastic bulk modulus
    bulkModulus_S   :: T                # Elastic bulk modulus
    shearModulus_S  :: T                # Elastic shear modulus

    # thermodynamic properties
    entropy         :: Vector{T}         # entropy
    enthalpy        :: Vector{T}        # enthalpy

    # Numerics:
    iter            :: I             # number of iterations required
    bulk_res_norm   :: T    # bulk residual norm
    time_ms         :: T          # computational time for this point
    status          :: I           # status of calculations
end

struct light_gmin_struct{T <: Float32, I <: Int8} 
    P_kbar      :: T                    # Pressure in kbar
    T_C         :: T                    # Temperature in Celsius
   
    ph_frac_wt  :: Vector{T}            # phase fractions
    ph_type     :: Vector{I}            # type of phase (SS or PP)
    ph_id_db    :: Vector{I}            # id of phase

    frac_S_wt   :: T
    frac_F_wt   :: T
    frac_M_wt   :: T

    bulk_S_wt   :: Vector{T}
    bulk_F_wt   :: Vector{T}
    bulk_M_wt   :: Vector{T}

    rho_S       :: T
    rho_F       :: T
    rho_M       :: T

    s_cp        :: Vector{T}
    alpha       :: Vector{T}
end

"""
    Holds general information about solution phases
"""
mutable struct ss_infos
    ss_fName:: String
    ss_name :: String
    n_em    :: Int64
    n_xeos  :: Int64
    n_sf    :: Int64
    ss_em   :: Vector{String}
    ss_xeos :: Vector{String}
    ss_sf   :: Vector{String}
end

"""
    Holds general information about the database
"""
mutable struct db_infos
    db_name :: String
    db_info :: String
    db_dataset :: Int64
    dataset_opt :: Union{Nothing, Int64, NTuple{5,Int64}, NTuple{4,Int64}}
    data_ss :: Array{ss_infos}
    ss_name :: Array{String}
    data_pp :: Array{String}
end
        

"""
    Holds the MAGEMin databases & required structures for every thread
"""
mutable struct MAGEMin_Data{TypeGV, TypeZB, TypeDB, TypeSplxData}
    db          :: String
    gv          :: TypeGV
    z_b         :: TypeZB
    DB          :: TypeDB
    splx_data   :: TypeSplxData
end

"""
    Holds the overriding Ws parameters
0 = "mp", 1 = "mb", 11 = "mbe",2 = "ig", 3 = "igad", 4 = "um", 5 = "ume", 6 = "mtl", 7 = "mpe", 8 = "sb11", 9 = "sb21"
"""
mutable struct W_data{T <: Float64,I <: Int64}
    dtb         :: I
    ss_ids      :: I
    n_Ws        :: I
    Ws          :: Matrix{T}   #S T P * n_Ws
end


"""
    Function to retrieve the general information of the databases
"""
function retrieve_solution_phase_information(dtb)

    db_inf  = db_infos[db_infos("mp", "Metapelite (White et al., 2014)", 62, (62, 633, 634, 635, 636), ss_infos[ss_infos("liq_W14", "liq", 8, 7, 10, ["none", "q4L", "abL", "kspL", "anL", "slL", "fo2L", "fa2L", "h2oL"], ["none", "q", "fsp", "na", "an", "ol", "x", "h2o"], ["none", "fac", "pq", "xab", "xksp", "pan", "psil", "pol", "xFe", "xMg", "ph2o"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("bi_W14", "bi", 7, 6, 13, ["none", "phl", "annm", "obi", "east", "tbi", "fbi", "mmbi"], ["none", "x", "m", "y", "f", "t", "Q"], ["none", "xMgM3", "xMnM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xMnM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"]), ss_infos("g_W14", "g", 5, 4, 6, ["none", "py", "alm", "spss", "gr", "kho"], ["none", "x", "z", "m", "f"], ["none", "xMgX", "xFeX", "xMnX", "xCaX", "xAlY", "xFe3Y"]), ss_infos("ep_H11", "ep", 3, 2, 4, ["none", "cz", "ep", "fep"], ["none", "f", "Q"], ["none", "xFeM1", "xAlM1", "xFeM3", "xAlM3"]), ss_infos("ma_W14", "ma", 6, 5, 10, ["none", "mut", "celt", "fcelt", "pat", "ma", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("mu_W14", "mu", 6, 5, 10, ["none", "mut", "cel", "fcel", "pat", "ma", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("opx_W14", "opx", 7, 6, 11, ["none", "en", "fs", "fm", "mgts", "fopx", "mnopx", "odi"], ["none", "x", "m", "y", "f", "c", "Q"], ["none", "xMgM1", "xFeM1", "xMnM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xMnM2", "xCaM2", "xSiT", "xAlT"]), ss_infos("sa_W14", "sa", 5, 4, 8, ["none", "spr4", "spr5", "fspm", "spro", "ospr"], ["none", "x", "y", "f", "Q"], ["none", "xMgM3", "xFeM3", "xFe3M3", "xAlM3", "xMgM456", "xFeM456", "xSiT", "xAlT"]), ss_infos("cd_W14", "cd", 4, 3, 5, ["none", "crd", "fcrd", "hcrd", "mncd"], ["none", "x", "m", "h"], ["none", "xFeX", "xMgX", "xMnX", "xH2OH", "xvH"]), ss_infos("st_W14", "st", 5, 4, 7, ["none", "mstm", "fst", "mnstm", "msto", "mstt"], ["none", "x", "m", "f", "t"], ["none", "xMgX", "xFeX", "xMnX", "xAlY", "xFe3Y", "xTiY", "xvY"]), ss_infos("chl_W14", "chl", 8, 7, 12, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin", "mmchl"], ["none", "x", "y", "f", "m", "QAl", "Q1", "Q4"], ["none", "xMgM1", "xMnM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("ctd_W14", "ctd", 4, 3, 5, ["none", "mctd", "fctd", "mnct", "ctdo"], ["none", "x", "m", "f"], ["none", "xAlM1A", "xFe3M1A", "xFeM1B", "xMgM1B", "xMnM1B"]), ss_infos("sp_W02", "sp", 4, 3, 5, ["none", "herc", "sp", "mt", "usp"], ["none", "x", "y", "z"], ["none", "xAl", "xFe3", "xTi", "xMg", "xFe2"]), ss_infos("mt_W00", "mt", 3, 2, 5, ["none", "imt", "dmt", "usp"], ["none", "x", "Q"], ["none", "xTiM", "xFe3M", "xFeM", "xFe3T", "xFeT"]), ss_infos("ilm_W00", "ilm", 3, 2, 6, ["none", "oilm", "dilm", "dhem"], ["none", "x", "Q"], ["none", "xFe2A", "xTiA", "xFe3A", "xFe2B", "xTiB", "xFe3B"]), ss_infos("ilmm_W14", "ilmm", 5, 4, 7, ["none", "oilm", "dilm", "dhem", "geik", "pnt"], ["none", "i", "g", "m", "Q"], ["none", "xFeA", "xTiA", "xMgA", "xMnA", "xFe3A", "xFeB", "xTiB"])], ["liq", "fsp", "bi", "g", "ep", "ma", "mu", "opx", "sa", "cd", "st", "chl", "ctd", "sp", "mt", "ilm", "ilmm"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "ru", "sph", "O2", "H2O", "zo", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("mb", "Metabasite (Green et al., 2016)", 62, (62, 633, 634, 635, 636), ss_infos[ss_infos("sp_W02", "sp", 4, 3, 5, ["none", "herc", "sp", "mt", "usp"], ["none", "x", "y", "z"], ["none", "xAl", "xFe3", "xTi", "xMg", "xFe2"]), ss_infos("opx_W14", "opx", 6, 5, 9, ["none", "en", "fs", "fm", "mgts", "fopx", "odi"], ["none", "x", "y", "f", "c", "Q"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xCaM2", "xAlT", "xSiT"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("liq_G16", "liq", 9, 8, 11, ["none", "q4L", "abL", "kspL", "wo1L", "sl1L", "fa2L", "fo2L", "h2oL", "anoL"], ["none", "q", "fsp", "na", "wo", "sil", "ol", "x", "yan"], ["none", "fac", "pq", "xab", "xksp", "pwo", "psil", "ph2o", "pan", "pol", "xFe", "xMg"]), ss_infos("mu_W14", "mu", 6, 5, 10, ["none", "mu", "cel", "fcel", "pa", "mam", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("ilmm_W14", "ilmm", 4, 3, 7, ["none", "oilm", "dilm", "dhem", "geik"], ["none", "c", "t", "Q"], ["none", "xFeA", "xTiA", "xMgA", "xFe3A", "xFeB", "xTiB", "xFe3B"]), ss_infos("ilm_W00", "ilm", 3, 2, 6, ["none", "oilm", "dilm", "dhem"], ["none", "x", "Q"], ["none", "xFe2A", "xTiA", "xFe3A", "xFe2B", "xTiB", "xFe3B"]), ss_infos("ol_H11", "ol", 2, 1, 2, ["none", "fo", "fa"], ["none", "x"], ["none", "xMgM", "xFeM"]), ss_infos("amp_G16", "amp", 11, 10, 18, ["none", "tr", "tsm", "prgm", "glm", "cumm", "grnm", "a", "b", "mrb", "kprg", "tts"], ["none", "x", "y", "z", "a", "k", "c", "f", "t", "Q1", "Q2"], ["none", "xvA", "xNaA", "xKA", "xMgM13", "xFeM13", "xMgM2", "xFeM2", "xAlM2", "xFe3M2", "xTiM2", "xCaM4", "xMgM4", "xFeM4", "xNaM4", "xSiT1", "xAlT1", "xOHV", "xOV"]), ss_infos("ep_H11", "ep", 3, 2, 4, ["none", "cz", "ep", "fep"], ["none", "f", "Q"], ["none", "xFeM1", "xAlM1", "xFeM3", "xAlM3"]), ss_infos("g_W14", "g", 4, 3, 5, ["none", "py", "alm", "gr", "kho"], ["none", "x", "z", "f"], ["none", "xMgX", "xFeX", "xCaX", "xAlY", "xFe3Y"]), ss_infos("chl_W14", "chl", 7, 6, 11, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin"], ["none", "x", "y", "f", "QAl", "Q1", "Q4"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("bi_W14", "bi", 6, 5, 11, ["none", "phl", "annm", "obi", "east", "tbi", "fbi"], ["none", "x", "y", "f", "t", "Q"], ["none", "xMgM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"]), ss_infos("dio_G16", "dio", 7, 6, 12, ["none", "jd", "di", "hed", "acmm", "om", "cfm", "jac"], ["none", "x", "j", "t", "c", "Qaf", "Qfm"], ["none", "xMgM1m", "xFeM1m", "xFe3M1m", "xAlM1m", "xMgM1a", "xFeM1a", "xFe3M1a", "xAlM1a", "xNaM2c", "xCaM2c", "xNaM2n", "xCaM2n"]), ss_infos("aug_G16", "aug", 8, 7, 12, ["none", "di", "cenh", "cfs", "jdm", "acmm", "ocats", "dcats", "fmc"], ["none", "x", "y", "f", "z", "j", "Qfm", "Qa1"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT1", "xAlT1", "xSiT2", "xAlT2"]), ss_infos("abc_H11", "abc", 2, 1, 2, ["none", "abm", "anm"], ["none", "ca"], ["none", "xNaA", "xCaA"]), ss_infos("spl_W02", "spl", 3, 2, 4, ["none", "herc", "sp", "usp"], ["none", "x", "y"], ["none", "xAl", "xTi", "xMg", "xFe2"])], ["sp", "opx", "fsp", "liq", "mu", "ilmm", "ilm", "ol", "amp", "ep", "g", "chl", "bi", "dio", "aug", "abc", "spl"], ["q", "crst", "trd", "coe", "law", "ky", "sill", "and", "ru", "sph", "ab", "H2O", "zo", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("mbe", "Metabasite extended (Green et al., 2016 with oamp from Diener et al., 2007 and ta from Rebay et al., 2022)", 62, (62, 633, 634, 635, 636), ss_infos[ss_infos("sp_W02", "sp", 4, 3, 5, ["none", "herc", "sp", "mt", "usp"], ["none", "x", "y", "z"], ["none", "xAl", "xFe3", "xTi", "xMg", "xFe2"]), ss_infos("opx_W14", "opx", 6, 5, 9, ["none", "en", "fs", "fm", "mgts", "fopx", "odi"], ["none", "x", "y", "f", "c", "Q"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xCaM2", "xAlT", "xSiT"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("liq_G16", "liq", 9, 8, 11, ["none", "q4L", "abL", "kspL", "wo1L", "sl1L", "fa2L", "fo2L", "h2oL", "anoL"], ["none", "q", "fsp", "na", "wo", "sil", "ol", "x", "yan"], ["none", "fac", "pq", "xab", "xksp", "pwo", "psil", "ph2o", "pan", "pol", "xFe", "xMg"]), ss_infos("mu_W14", "mu", 6, 5, 10, ["none", "mu", "cel", "fcel", "pa", "mam", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("ilmm_W14", "ilmm", 4, 3, 7, ["none", "oilm", "dilm", "dhem", "geik"], ["none", "c", "t", "Q"], ["none", "xFeA", "xTiA", "xMgA", "xFe3A", "xFeB", "xTiB", "xFe3B"]), ss_infos("ilm_W00", "ilm", 3, 2, 6, ["none", "oilm", "dilm", "dhem"], ["none", "x", "Q"], ["none", "xFe2A", "xTiA", "xFe3A", "xFe2B", "xTiB", "xFe3B"]), ss_infos("ol_H11", "ol", 2, 1, 2, ["none", "fo", "fa"], ["none", "x"], ["none", "xMgM", "xFeM"]), ss_infos("amp_G16", "amp", 11, 10, 18, ["none", "tr", "tsm", "prgm", "glm", "cumm", "grnm", "a", "b", "mrb", "kprg", "tts"], ["none", "x", "y", "z", "a", "k", "c", "f", "t", "Q1", "Q2"], ["none", "xvA", "xNaA", "xKA", "xMgM13", "xFeM13", "xMgM2", "xFeM2", "xAlM2", "xFe3M2", "xTiM2", "xCaM4", "xMgM4", "xFeM4", "xNaM4", "xSiT1", "xAlT1", "xOHV", "xOV"]), ss_infos("ep_H11", "ep", 3, 2, 4, ["none", "cz", "ep", "fep"], ["none", "f", "Q"], ["none", "xFeM1", "xAlM1", "xFeM3", "xAlM3"]), ss_infos("g_W14", "g", 4, 3, 5, ["none", "py", "alm", "gr", "kho"], ["none", "x", "z", "f"], ["none", "xMgX", "xFeX", "xCaX", "xAlY", "xFe3Y"]), ss_infos("chl_W14", "chl", 7, 6, 11, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin"], ["none", "x", "y", "f", "QAl", "Q1", "Q4"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("bi_W14", "bi", 6, 5, 11, ["none", "phl", "annm", "obi", "east", "tbi", "fbi"], ["none", "x", "y", "f", "t", "Q"], ["none", "xMgM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"]), ss_infos("dio_G16", "dio", 7, 6, 12, ["none", "jd", "di", "hed", "acmm", "om", "cfm", "jac"], ["none", "x", "j", "t", "c", "Qaf", "Qfm"], ["none", "xMgM1m", "xFeM1m", "xFe3M1m", "xAlM1m", "xMgM1a", "xFeM1a", "xFe3M1a", "xAlM1a", "xNaM2c", "xCaM2c", "xNaM2n", "xCaM2n"]), ss_infos("aug_G16", "aug", 8, 7, 12, ["none", "di", "cenh", "cfs", "jdm", "acmm", "ocats", "dcats", "fmc"], ["none", "x", "y", "f", "z", "j", "Qfm", "Qa1"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT1", "xAlT1", "xSiT2", "xAlT2"]), ss_infos("abc_H11", "abc", 2, 1, 2, ["none", "abm", "anm"], ["none", "ca"], ["none", "xNaA", "xCaA"]), ss_infos("spl_W02", "spl", 3, 2, 4, ["none", "herc", "sp", "usp"], ["none", "x", "y"], ["none", "xAl", "xTi", "xMg", "xFe2"]), ss_infos("", "ta", 5, 4, 8, ["none", "ta", "fta", "ota", "tap", "tats"], ["none", "x", "y", "z", "q"], ["none", "xvM1", "xMgM1", "xFeM1", "xMgM23", "xFeM23", "xAlM23", "xSiT1", "xAlT1"]), ss_infos("", "oamp", 9, 8, 14, ["none", "anth", "ged", "ompa", "omgl", "otr", "fanth", "omrb", "amoa", "amob"], ["none", "x", "y", "z", "a", "c", "f", "q1", "q2"], ["none", "xvA", "xNaA", "xCaM4", "xNaM4", "xMgM4", "xFeM4", "xMgM13", "xFeM13", "xAlM2", "xFe3M2", "xMgM2", "xFeM2", "xAlT1", "xSiT1"])], ["sp", "opx", "fsp", "liq", "mu", "ilmm", "ilm", "ol", "amp", "ep", "g", "chl", "bi", "dio", "aug", "abc", "spl", "ta", "oamp"], ["q", "crst", "trd", "coe", "law", "ky", "sill", "and", "ru", "sph", "ab", "H2O", "zo", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("ig", "Igneous (Green et al., 2025, corrected after Holland et al., 2018)", 636, (62, 633, 634, 635, 636), ss_infos[ss_infos("spl_T21", "spl", 8, 7, 10, ["none", "nsp", "isp", "nhc", "ihc", "nmt", "imt", "pcr", "qndm"], ["none", "x", "y", "c", "t", "Q1", "Q2", "Q3"], ["none", "xMgT", "xFeT", "xAlT", "xFe3T", "xMgM", "xFeM", "xAlM", "xFe3M", "xCrM", "xTiM"]), ss_infos("bi_G25", "bi", 6, 5, 11, ["none", "phl", "annm", "obi", "eas", "tbi", "fbi"], ["none", "x", "y", "f", "t", "Q"], ["none", "xMgM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"]), ss_infos("cd_G25", "cd", 3, 2, 4, ["none", "crd", "fcrd", "hcrd"], ["none", "x", "h"], ["none", "xFeX", "xMgX", "xH2OH", "xvH"]), ss_infos("cpx_W24", "cpx", 10, 9, 13, ["none", "di", "cfs", "cats", "crdi", "cess", "cbuf", "jd", "cen", "cfm", "kjd"], ["none", "x", "y", "o", "n", "Q", "f", "cr", "t", "k"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xCrM1", "xTiM1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xKM2", "xSiT", "xAlT"]), ss_infos("ep_H11", "ep", 3, 2, 4, ["none", "cz", "ep", "fep"], ["none", "f", "Q"], ["none", "xFeM1", "xAlM1", "xFeM3", "xAlM3"]), ss_infos("g_W24", "g", 6, 5, 8, ["none", "py", "alm", "gr", "andr", "knom", "tig"], ["none", "x", "c", "f", "cr", "t"], ["none", "xMgM1", "xFeM1", "xCaM1", "xAlM2", "xCrM2", "xFe3M2", "xMgM2", "xTiM2"]), ss_infos("amp_G16", "amp", 11, 10, 18, ["none", "tr", "tsm", "prgm", "glm", "cumm", "grnm", "a", "b", "mrb", "kprg", "tts"], ["none", "x", "y", "z", "a", "k", "c", "f", "t", "Q1", "Q2"], ["none", "xvA", "xNaA", "xKA", "xMgM13", "xFeM13", "xMgM2", "xFeM2", "xAlM2", "xFe3M2", "xTiM2", "xCaM4", "xMgM4", "xFeM4", "xNaM4", "xSiT1", "xAlT1", "xOHV", "xOV"]), ss_infos("ilm_W24", "ilm", 5, 4, 8, ["none", "oilm", "dilm", "hm", "ogk", "dgk"], ["none", "i", "m", "Q", "Qt"], ["none", "xFeA", "xTiA", "xFe3A", "xMgA", "xFeB", "xTiB", "xFe3B", "xMgB"]), ss_infos("liq_G25w", "liq", 12, 11, 18, ["none", "q4L", "slL", "wo1L", "fo2L", "fa2L", "jdL", "hmL", "ekL", "tiL", "kjL", "ctL", "h2o1L"], ["none", "wo", "sl", "fo", "fa", "jd", "hm", "ek", "ti", "kj", "yct", "h2o"], ["none", "pq", "psl", "pwo", "pjd", "phm", "pek", "pti", "pkj", "pct", "pol", "sumT", "mgM", "feM", "CaM", "AlM", "sumM", "xh", "xv"]), ss_infos("ol_H18", "ol", 4, 3, 5, ["none", "mont", "fa", "fo", "cfm"], ["none", "x", "c", "Q"], ["none", "xMgM1", "xFeM1", "xMgM2", "xFeM2", "xCaM2"]), ss_infos("opx_W24", "opx", 9, 8, 12, ["none", "en", "fs", "fm", "odi", "mgts", "cren", "obuf", "mess", "ojd"], ["none", "x", "y", "c", "Q", "f", "t", "cr", "j"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xCrM1", "xTiM1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT", "xAlT"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("fl_G25", "fl", 11, 10, 12, ["none", "qfL", "slfL", "wofL", "fofL", "fafL", "jdfL", "hmfL", "ekfL", "tifL", "kjfL", "H2O"], ["none", "wo", "sl", "fo", "fa", "jd", "hm", "ek", "ti", "kj", "h2o"], ["none", "pq", "psl", "pwo", "pfo", "pfa", "pjd", "phm", "pek", "pti", "pkj", "ph2o", "fac"]), ss_infos("mu_W14", "mu", 6, 5, 10, ["none", "mu", "cel", "fcel", "pa", "mam", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("fper", "fper", 2, 1, 2, ["none", "per", "wu"], ["none", "x"], ["none", "xFe", "xMg"]), ss_infos("chl_W14", "chl", 7, 6, 11, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin"], ["none", "x", "y", "f", "QAl", "Q1", "Q4"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"])], ["spl", "bi", "cd", "cpx", "ep", "g", "amp", "ilm", "liq", "ol", "opx", "fsp", "fl", "mu", "fper", "chl"], ["ne", "q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "ru", "sph", "O2", "H2O", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("igad", "Igneous alkaline dry (Weller et al., 2024)", 636, (62, 633, 634, 635, 636), ss_infos[ss_infos("spl_T21", "spl", 8, 7, 10, ["none", "nsp", "isp", "nhc", "ihc", "nmt", "imt", "pcr", "usp"], ["none", "x", "y", "c", "t", "Q1", "Q2", "Q3"], ["none", "xMgT", "xFeT", "xAlT", "xFe3T", "xMgM", "xFeM", "xAlM", "xFe3M", "xCrM", "xTiM"]), ss_infos("cpx_W24", "cpx", 10, 9, 13, ["none", "di", "cfs", "cats", "crdi", "cess", "cbuf", "jd", "cen", "cfm", "kjd"], ["none", "x", "y", "o", "n", "Q", "f", "cr", "t", "k"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xCrM1", "xTiM1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xKM2", "xSiT", "xAlT"]), ss_infos("g_W24", "g", 6, 5, 8, ["none", "py", "alm", "gr", "andr", "knr", "tig"], ["none", "x", "c", "f", "cr", "t"], ["none", "xMgM1", "xFeM1", "xCaM1", "xAlM2", "xCrM2", "xFe3M2", "xMgM2", "xTiM2"]), ss_infos("ilm_W24", "ilm", 5, 4, 8, ["none", "oilm", "dilm", "hm", "ogk", "dgk"], ["none", "i", "m", "Q", "Qt"], ["none", "xFeA", "xTiA", "xFe3A", "xMgA", "xFeB", "xTiB", "xFe3B", "xMgB"]), ss_infos("liq_W24d", "liq", 14, 13, 18, ["none", "q3L", "sl1L", "wo1L", "fo2L", "fa2L", "nmL", "hmL", "ekL", "tiL", "kmL", "anL", "ab1L", "enL", "kfL"], ["none", "wo", "sl", "fo", "fa", "ns", "hm", "ek", "ti", "ks", "yan", "yab", "yen", "ykf"], ["none", "pq", "psl", "pwo", "pns", "phm", "pek", "pti", "pks", "pab", "pan", "pen", "pkf", "pol", "mgM", "feM", "CaM", "AlM", "sumM"]), ss_infos("ol_H18", "ol", 4, 3, 5, ["none", "mnt", "fa", "fo", "cfm"], ["none", "x", "c", "Q"], ["none", "xMgM1", "xFeM1", "xMgM2", "xFeM2", "xCaM2"]), ss_infos("opx_W24", "opx", 9, 8, 12, ["none", "en", "fs", "fm", "odi", "mgts", "cren", "obuf", "mess", "ojd"], ["none", "x", "y", "c", "Q", "f", "t", "cr", "j"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xCrM1", "xTiM1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT", "xAlT"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("lct_W24", "lct", 2, 1, 2, ["none", "nlc", "klc"], ["none", "n"], ["none", "xNaA", "xKA"]), ss_infos("mel_W24", "mel", 5, 4, 8, ["none", "geh", "ak", "fak", "nml", "fge"], ["none", "x", "n", "y", "f"], ["none", "xNaM1", "xCaM1", "xMgT1", "xFeT1", "xAlT1", "xFe3T1", "xAlT2", "xSiT2"]), ss_infos("nph_W24", "nph", 6, 5, 9, ["none", "neN", "neS", "neK", "neO", "neC", "neF"], ["none", "s", "k", "Q", "f", "c"], ["none", "xNaA1", "xKA1", "xCaA1", "xNaA2", "xKA2", "xvA2", "xAlT2", "xSiT2", "xFe3T2"]), ss_infos("kals_W24", "kals", 2, 1, 2, ["none", "nks", "kls"], ["none", "k"], ["none", "xKA", "xNaA"])], ["spl", "cpx", "g", "ilm", "liq", "ol", "opx", "fsp", "lct", "mel", "nph", "kals"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "ru", "sph", "O2", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("um", "Ultramafic (Evans & Frost., 2021)", 633, (62, 633, 634, 635, 636), ss_infos[ss_infos("fl_EF21", "fl", 2, 1, 2, ["none", "H2", "H2O"], ["none", "x"], ["none", "xH2", "xH2O"]), ss_infos("ol_H11", "ol", 2, 1, 2, ["none", "fo", "fa"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("br_E13", "br", 2, 1, 2, ["none", "br", "fbr"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("ch_EF21", "ch", 2, 1, 2, ["none", "chum", "chuf"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("atg_EF21", "atg", 5, 4, 8, ["none", "atgf", "fatg", "atgo", "aatg", "oatg"], ["none", "x", "y", "f", "t"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xSiT", "xAlT"]), ss_infos("g_H18", "g", 2, 1, 2, ["none", "py", "alm"], ["none", "x"], ["none", "xMgM1", "xFeM1"]), ss_infos("ta_EF21", "ta", 6, 5, 9, ["none", "ta", "fta", "tao", "tats", "ota", "tap"], ["none", "x", "y", "f", "v", "Q"], ["none", "xMgM1", "xFeM1", "xvM1", "xMgM23", "xFeM23", "xFe3M23", "xAlM23", "xSiT2", "xAlT2"]), ss_infos("chl_W14", "chl", 7, 6, 11, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin"], ["none", "x", "y", "f", "m", "t", "QA1"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("spi_W02", "spi", 3, 2, 4, ["none", "herc", "sp", "mt"], ["none", "x", "y"], ["none", "xAl", "xFe3", "xMg", "xFe2"]), ss_infos("opx_W14", "opx", 5, 4, 8, ["none", "en", "fs", "fm", "mgts", "fopx"], ["none", "x", "y", "f", "Q"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xAlT", "xSiT"]), ss_infos("po_E10", "po", 2, 1, 2, ["none", "trov", "trot"], ["none", "y"], ["none", "xfeM2", "xVM2"]), ss_infos("anth_D07", "anth", 5, 4, 9, ["none", "anth", "gedf", "fant", "a", "b"], ["none", "x", "y", "z", "a"], ["none", "xMgM4", "xFeM4", "xMgM13", "xFeM13", "xAlM2", "xMgM2", "xFeM2", "xAlT1", "xSiT1"])], ["fl", "ol", "br", "ch", "atg", "g", "ta", "chl", "spi", "opx", "po", "anth"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "pyr", "O2", "hem", "qfm", "qif", "nno", "hm", "mw", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("ume", "Ultramafic extended (Evans & Frost., 2021 with pl, amp and aug from Green et al., 2016)", 633, (62, 633, 634, 635, 636), ss_infos[ss_infos("fl_EF21", "fl", 2, 1, 2, ["none", "H2", "H2O"], ["none", "x"], ["none", "xH2", "xH2O"]), ss_infos("ol_H11", "ol", 2, 1, 2, ["none", "fo", "fa"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("br_E13", "br", 2, 1, 2, ["none", "br", "fbr"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("ch_EF21", "ch", 2, 1, 2, ["none", "chum", "chuf"], ["none", "x"], ["none", "xMg", "xFe"]), ss_infos("atg_EF21", "atg", 5, 4, 8, ["none", "atgf", "fatg", "atgo", "aatg", "oatg"], ["none", "x", "y", "f", "t"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xSiT", "xAlT"]), ss_infos("g_H18", "g", 2, 1, 2, ["none", "py", "alm"], ["none", "x"], ["none", "xMgM1", "xFeM1"]), ss_infos("ta_EF21", "ta", 6, 5, 9, ["none", "ta", "fta", "tao", "tats", "ota", "tap"], ["none", "x", "y", "f", "v", "Q"], ["none", "xMgM1", "xFeM1", "xvM1", "xMgM23", "xFeM23", "xFe3M23", "xAlM23", "xSiT2", "xAlT2"]), ss_infos("chl_W14", "chl", 7, 6, 11, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin"], ["none", "x", "y", "f", "m", "t", "QA1"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("spi_W02", "spi", 3, 2, 4, ["none", "herc", "sp", "mt"], ["none", "x", "y"], ["none", "xAl", "xFe3", "xMg", "xFe2"]), ss_infos("opx_W14", "opx", 5, 4, 8, ["none", "en", "fs", "fm", "mgts", "fopx"], ["none", "x", "y", "f", "Q"], ["none", "xMgM1", "xFeM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xAlT", "xSiT"]), ss_infos("po_E10", "po", 2, 1, 2, ["none", "trov", "trot"], ["none", "y"], ["none", "xfeM2", "xVM2"]), ss_infos("anth_D07", "anth", 5, 4, 9, ["none", "anth", "gedf", "fant", "a", "b"], ["none", "x", "y", "z", "a"], ["none", "xMgM4", "xFeM4", "xMgM13", "xFeM13", "xAlM2", "xMgM2", "xFeM2", "xAlT1", "xSiT1"]), ss_infos("fsp_H22", "pl4tr", 2, 1, 4, ["none", "ab", "an"], ["none", "ca"], ["none", "xNaA", "xCaA", "xAlTB", "xSiTB"]), ss_infos("amp_G16", "amp", 9, 8, 14, ["none", "tr", "tsm", "prgm", "glm", "cumm", "grnm", "a", "b", "mrb"], ["none", "x", "y", "z", "a", "c", "f", "Q1", "Q2"], ["none", "xvA", "xNaA", "xMgM13", "xFeM13", "xMgM2", "xFeM2", "xAlM2", "xFe3M2", "xCaM4", "xMgM4", "xFeM4", "xNaM4", "xSiT1", "xAlT1"]), ss_infos("aug_G16", "aug", 8, 7, 12, ["none", "di", "cenh", "cfs", "jdm", "acmm", "ocats", "dcats", "fmc"], ["none", "x", "y", "f", "z", "j", "Qfm", "Qa1"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT1", "xAlT1", "xSiT2", "xAlT2"])], ["fl", "ol", "br", "ch", "atg", "g", "ta", "chl", "spi", "opx", "po", "anth", "pl4tr", "amp", "aug"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "pyr", "O2", "hem", "qfm", "qif", "nno", "hm", "mw", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("mtl", "Mantle (Holland et al., 2013)", 633, (62, 633, 634, 635, 636), ss_infos[ss_infos("g_H13", "g", 6, 5, 8, ["none", "py", "alm", "gr", "maj", "gfm", "nagt"], ["none", "x", "c", "y", "Q", "n"], ["none", "xMgM1", "xFeM1", "xCaM1", "xNaM1", "xAlM2", "xMgM2", "xFeM2", "xSiM2"]), ss_infos("fp_H13", "fp", 2, 1, 2, ["none", "per", "fper"], ["none", "x"], ["none", "xMgM1", "xFeM1"]), ss_infos("mpv_H13", "mpv", 5, 4, 7, ["none", "mpv", "fpvm", "cpvm", "apv", "npvm"], ["none", "x", "y", "c", "n"], ["none", "xCaM1", "xMgM1", "xFeM1", "xNaM1", "xAlM1", "xAlM2", "xSiM2"]), ss_infos("cpv_H13", "cpv", 5, 4, 7, ["none", "mpv", "fpvm", "cpvm", "apv", "npvm"], ["none", "x", "y", "c", "n"], ["none", "xCaM1", "xMgM1", "xFeM1", "xNaM1", "xAlM1", "xAlM2", "xSiM2"]), ss_infos("crn_H13", "crn", 3, 2, 5, ["none", "cor", "mcor", "fcor"], ["none", "x", "y"], ["none", "xMgM1", "xFeM1", "xAlM1", "xAlM2", "xSiM2"]), ss_infos("cf_H13", "cf", 6, 5, 8, ["none", "macf", "cacf", "mscf", "fscf", "oscf", "nacfm"], ["none", "y", "x", "Q", "c", "n"], ["none", "xCaM1", "xMgM1", "xFeM1", "xNaM1", "xMgM2", "xFeM2", "xAlM2", "xSiM2"]), ss_infos("nal_H13", "nal", 7, 6, 10, ["none", "nanal", "canal", "manal", "msnal", "fsnal", "o1nal", "o2nal"], ["none", "y", "x", "Q1", "Q2", "c", "n"], ["none", "xCaM3", "xMgM3", "xFeM3", "xNaM3", "xMgM2", "xFeM2", "xMgM1", "xFeM1", "xAlM1", "xSiM1"]), ss_infos("aki_H13", "aki", 3, 2, 5, ["none", "aak", "mak", "fak"], ["none", "x", "y"], ["none", "xAlA", "xMgA", "xFeA", "xAlB", "xSiB"]), ss_infos("ol_H13", "ol", 2, 1, 2, ["none", "fo", "fa"], ["none", "x"], ["none", "pfo", "pfa"]), ss_infos("wad_H13", "wad", 2, 1, 2, ["none", "mwd", "fwd"], ["none", "x"], ["none", "pmwd", "pfwd"]), ss_infos("ring_H13", "ring", 2, 1, 2, ["none", "mrw", "frw"], ["none", "x"], ["none", "pmrw", "pfrw"]), ss_infos("cpx_H13", "cpx", 6, 5, 9, ["none", "di", "cfs", "cats", "jd", "cen", "cfm"], ["none", "x", "y", "o", "n", "Q"], ["none", "xMgM1", "xFeM1", "xAlM1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT", "xAlT"]), ss_infos("opx_H13", "opx", 5, 4, 8, ["none", "en", "fs", "fm", "odi", "mgts"], ["none", "x", "y", "c", "Q"], ["none", "xMgM1", "xFeM1", "xAlM1", "xCaM2", "xMgM2", "xFeM2", "xSiT", "xAlT"]), ss_infos("hpx_H13", "hpx", 5, 4, 8, ["none", "en", "fs", "fm", "odi", "hmts"], ["none", "x", "y", "c", "Q"], ["none", "xMgM1", "xFeM1", "xAlM1", "xCaM2", "xMgM2", "xFeM2", "xSiT", "xAlT"])], ["g", "fp", "mpv", "cpv", "crn", "cf", "nal", "aki", "ol", "wad", "ring", "cpx", "opx", "hpx"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and"]), db_infos("mpe", "Metapelite extended (White et al., 2014 with po from Evans & Frost., 2021, amp dio and aug from Green et al., 2016)", 62, (62, 633, 634, 635, 636), ss_infos[ss_infos("liq_W14", "liq", 8, 7, 10, ["none", "q4L", "abL", "kspL", "anL", "slL", "fo2L", "fa2L", "h2oL"], ["none", "q", "fsp", "na", "an", "ol", "x", "h2o"], ["none", "fac", "pq", "xab", "xksp", "pan", "psil", "pol", "xFe", "xMg", "ph2o"]), ss_infos("fsp_H22", "fsp", 3, 2, 5, ["none", "ab", "an", "san"], ["none", "ca", "k"], ["none", "xNaA", "xCaA", "xKA", "xAlTB", "xSiTB"]), ss_infos("bi_W14", "bi", 7, 6, 13, ["none", "phl", "annm", "obi", "east", "tbi", "fbi", "mmbi"], ["none", "x", "m", "y", "f", "t", "Q"], ["none", "xMgM3", "xMnM3", "xFeM3", "xFe3M3", "xTiM3", "xAlM3", "xMgM12", "xMnM12", "xFeM12", "xSiT", "xAlT", "xOHV", "xOV"]), ss_infos("g_W14", "g", 5, 4, 6, ["none", "py", "alm", "spss", "gr", "kho"], ["none", "x", "z", "m", "f"], ["none", "xMgX", "xFeX", "xMnX", "xCaX", "xAlY", "xFe3Y"]), ss_infos("ep_H11", "ep", 3, 2, 4, ["none", "cz", "ep", "fep"], ["none", "f", "Q"], ["none", "xFeM1", "xAlM1", "xFeM3", "xAlM3"]), ss_infos("ma_W14", "ma", 6, 5, 10, ["none", "mut", "celt", "fcelt", "pat", "ma", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("mu_W14", "mu", 6, 5, 10, ["none", "mut", "cel", "fcel", "pat", "ma", "fmu"], ["none", "x", "y", "f", "n", "c"], ["none", "xKA", "xNaA", "xCaA", "xMgM2A", "xFeM2A", "xAlM2A", "xAlM2B", "xFe3M2B", "xSiT1", "xAlT1"]), ss_infos("opx_W14", "opx", 7, 6, 11, ["none", "en", "fs", "fm", "mgts", "fopx", "mnopx", "odi"], ["none", "x", "m", "y", "f", "c", "Q"], ["none", "xMgM1", "xFeM1", "xMnM1", "xFe3M1", "xAlM1", "xMgM2", "xFeM2", "xMnM2", "xCaM2", "xSiT", "xAlT"]), ss_infos("sa_W14", "sa", 5, 4, 8, ["none", "spr4", "spr5", "fspm", "spro", "ospr"], ["none", "x", "y", "f", "Q"], ["none", "xMgM3", "xFeM3", "xFe3M3", "xAlM3", "xMgM456", "xFeM456", "xSiT", "xAlT"]), ss_infos("cd_W14", "cd", 4, 3, 5, ["none", "crd", "fcrd", "hcrd", "mncd"], ["none", "x", "m", "h"], ["none", "xFeX", "xMgX", "xMnX", "xH2OH", "xvH"]), ss_infos("st_W14", "st", 5, 4, 7, ["none", "mstm", "fst", "mnstm", "msto", "mstt"], ["none", "x", "m", "f", "t"], ["none", "xMgX", "xFeX", "xMnX", "xAlY", "xFe3Y", "xTiY", "xvY"]), ss_infos("chl_W14", "chl", 8, 7, 12, ["none", "clin", "afchl", "ames", "daph", "ochl1", "ochl4", "f3clin", "mmchl"], ["none", "x", "y", "f", "m", "QAl", "Q1", "Q4"], ["none", "xMgM1", "xMnM1", "xFeM1", "xAlM1", "xMgM23", "xFeM23", "xMgM4", "xFeM4", "xFe3M4", "xAlM4", "xSiT2", "xAlT2"]), ss_infos("ctd_W14", "ctd", 4, 3, 5, ["none", "mctd", "fctd", "mnct", "ctdo"], ["none", "x", "m", "f"], ["none", "xAlM1A", "xFe3M1A", "xFeM1B", "xMgM1B", "xMnM1B"]), ss_infos("sa_W02", "sp", 4, 3, 5, ["none", "herc", "sp", "mt", "usp"], ["none", "x", "y", "z"], ["none", "xAl", "xFe3", "xTi", "xMg", "xFe2"]), ss_infos("mt_W00", "mt", 3, 2, 5, ["none", "imt", "dmt", "usp"], ["none", "x", "Q"], ["none", "xTiM", "xFe3M", "xFeM", "xFe3T", "xFeT"]), ss_infos("ilm_W00", "ilm", 3, 2, 6, ["none", "oilm", "dilm", "dhem"], ["none", "x", "Q"], ["none", "xFe2A", "xTiA", "xFe3A", "xFe2B", "xTiB", "xFe3B"]), ss_infos("ilmm_W14", "ilmm", 5, 4, 7, ["none", "oilm", "dilm", "dhem", "geik", "pnt"], ["none", "i", "g", "m", "Q"], ["none", "xFeA", "xTiA", "xMgA", "xMnA", "xFe3A", "xFeB", "xTiB"]), ss_infos("occm_F11", "occm", 5, 4, 9, ["none", "cc", "odo", "mag", "sid", "oank"], ["none", "x", "j", "q", "v"], ["none", "xCaM1", "xMgM1", "xFeM1", "xCaM2a", "xMgM2a", "xFeM2a", "xCaM2b", "xMgM2b", "xFeM2b"]), ss_infos("fl_H03", "fl", 2, 1, 2, ["none", "H2O", "CO2"], ["none", "x"], ["none", "xH2O", "xCO2"]), ss_infos("po_E10", "po", 2, 1, 2, ["none", "trov", "trot"], ["none", "y"], ["none", "xfeM2", "xVM2"]), ss_infos("dio_G16", "dio", 7, 6, 12, ["none", "jd", "di", "hed", "acmm", "om", "cfm", "jac"], ["none", "x", "j", "t", "c", "qaf", "qfm"], ["none", "xMgM1m", "xFeM1m", "xFe3M1m", "xAlM1m", "xMgM1a", "xFeM1a", "xFe3M1a", "xAlM1a", "xNaM2c", "xCaM2c", "xNaM2n", "xCaM2n"]), ss_infos("aug_G16", "aug", 8, 7, 12, ["none", "di", "cenh", "cfs", "jdm", "acmm", "ocats", "dcats", "fmc"], ["none", "x", "y", "c", "z", "j", "qfm", "qal"], ["none", "xMgM1", "xFeM1", "xAlM1", "xFe3M1", "xMgM2", "xFeM2", "xCaM2", "xNaM2", "xSiT1", "xAlT1", "xSiT2", "xAlT2"]), ss_infos("amp_G16", "amp", 11, 10, 18, ["none", "tr", "tsm", "prgm", "glm", "cumm", "grnm", "a", "b", "mrb", "kprg", "tts"], ["none", "x", "y", "z", "a", "k", "c", "f", "t", "q1", "q2"], ["none", "xvA", "xNaA", "xKA", "xMgM13", "xFeM13", "xMgM2", "xFeM2", "xAlM2", "xFe3M2", "xTiM2", "xCaM4", "xMgM4", "xFeM4", "xNaM4", "xSiT1", "xAlT1", "xOHV", "xOV"]), ss_infos("", "oamp", 9, 8, 14, ["none", "anth", "ged", "ompa", "omgl", "otr", "fanth", "omrb", "amoa", "amob"], ["none", "x", "y", "z", "a", "c", "f", "q1", "q2"], ["none", "xvA", "xNaA", "xCaM4", "xNaM4", "xMgM4", "xFeM4", "xMgM13", "xFeM13", "xAlM2", "xFe3M2", "xMgM2", "xFeM2", "xAlT1", "xSiT1"])], ["liq", "fsp", "bi", "g", "ep", "ma", "mu", "opx", "sa", "cd", "st", "chl", "ctd", "sp", "mt", "ilm", "ilmm", "occm", "fl", "po", "dio", "aug", "amp", "oamp"], ["q", "crst", "trd", "coe", "stv", "ky", "sill", "and", "ru", "sph", "O2", "pyr", "gph", "law", "zo", "prl", "mpm", "pre", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]), db_infos("sb11", "Stixrude & Lithgow-Bertelloni (2011)", -1, (62, 633, 634, 635, 636), ss_infos[ss_infos("", "plg", 2, 2, 1, ["none", "an", "ab"], ["none", "", ""], ["none", ""]), ss_infos("", "sp", 2, 2, 2, ["none", "sp", "hc"], ["none", "", ""], ["none", "", ""]), ss_infos("", "ol", 2, 2, 1, ["none", "fo", "fa"], ["none", "", ""], ["none", ""]), ss_infos("", "wa", 2, 2, 1, ["none", "mgwa", "fewa"], ["none", "", ""], ["none", ""]), ss_infos("", "ri", 2, 2, 1, ["none", "mgri", "feri"], ["none", "", ""], ["none", ""]), ss_infos("", "opx", 4, 4, 2, ["none", "en", "fs", "mgts", "odi"], ["none", "", "", "", ""], ["none", "", ""]), ss_infos("", "cpx", 5, 5, 3, ["none", "di", "he", "cen", "cats", "jd"], ["none", "", "", "", "", ""], ["none", "", "", ""]), ss_infos("", "hpcpx", 2, 2, 1, ["none", "hpcen", "hpcfs"], ["none", "", ""], ["none", ""]), ss_infos("", "ak", 3, 3, 2, ["none", "mgak", "feak", "co"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "gtmj", 5, 5, 3, ["none", "py", "alm", "gr", "mgmj", "jdmj"], ["none", "", "", "", "", ""], ["none", "", "", ""]), ss_infos("", "pv", 3, 3, 2, ["none", "mgpv", "fepv", "alpv"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "ppv", 3, 3, 2, ["none", "mppv", "fppv", "appv"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "mw", 2, 2, 1, ["none", "pe", "wu"], ["none", "", ""], ["none", ""]), ss_infos("", "cf", 3, 3, 2, ["none", "mgcf", "fecf", "nacf"], ["none", "", "", ""], ["none", "", ""])], ["plg", "sp", "ol", "wa", "ri", "opx", "cpx", "hpcpx", "ak", "gtmj", "pv", "ppv", "mw", "cf"], ["neph", "ky", "st", "coe", "qtz", "capv", "co", "aMgO", "aFeO", "aAl2O3"]), db_infos("sb21", "Stixrude & Lithgow-Bertelloni (2021)", -1, (62, 633, 634, 635, 636), ss_infos[ss_infos("", "plg", 2, 2, 1, ["none", "an", "ab"], ["none", "", ""], ["none", ""]), ss_infos("", "sp", 2, 2, 2, ["none", "sp", "hc"], ["none", "", ""], ["none", "", ""]), ss_infos("", "ol", 2, 2, 1, ["none", "fo", "fa"], ["none", "", ""], ["none", ""]), ss_infos("", "wa", 2, 2, 1, ["none", "mgwa", "fewa"], ["none", "", ""], ["none", ""]), ss_infos("", "ri", 2, 2, 1, ["none", "mgri", "feri"], ["none", "", ""], ["none", ""]), ss_infos("", "opx", 4, 4, 2, ["none", "en", "fs", "mgts", "odi"], ["none", "", "", "", ""], ["none", "", ""]), ss_infos("", "cpx", 5, 5, 3, ["none", "di", "he", "cen", "cats", "jd"], ["none", "", "", "", "", ""], ["none", "", "", ""]), ss_infos("", "hpcpx", 2, 2, 1, ["none", "hpcen", "hpcfs"], ["none", "", ""], ["none", ""]), ss_infos("", "ak", 3, 3, 2, ["none", "mgak", "feak", "co"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "gtmj", 5, 5, 3, ["none", "py", "alm", "gr", "mgmj", "jdmj"], ["none", "", "", "", "", ""], ["none", "", "", ""]), ss_infos("", "pv", 3, 3, 2, ["none", "mgpv", "fepv", "alpv"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "ppv", 3, 3, 2, ["none", "mppv", "fppv", "appv"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "cf", 3, 3, 2, ["none", "mgcf", "fecf", "nacf"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "mw", 3, 3, 2, ["none", "pe", "wu", "anao"], ["none", "", "", ""], ["none", "", ""]), ss_infos("", "nal", 3, 3, 3, ["none", "mnal", "fnal", "nnal"], ["none", "", "", ""], ["none", "", "", ""])], ["plg", "sp", "ol", "wa", "ri", "opx", "cpx", "hpcpx", "ak", "gtmj", "pv", "ppv", "cf", "mw", "nal"], ["neph", "ky", "st", "coe", "qtz", "capv", "co", "aMgO", "aFeO", "aAl2O3"])]
    dbs     = ["mp","mb","mbe","ig","igad","um","ume","mtl","mpe","sb11","sb21"]
    id      = findall(dbs .== dtb)[1]

    return db_inf[id]
end

"""
    Function to retrieve the list of indexes of the solution phases to be removed from the minimization
"""
function remove_phases( list        :: Union{Nothing,Vector{String}},
                        dtb         :: String)

    if ~isnothing(list)
        db_inf = retrieve_solution_phase_information(dtb);
        rm_list = zeros(Int64,0);
        for i in list

            if ~isnothing(i)
                if i in db_inf.data_pp
                    idx = findfirst(db_inf.data_pp .== i)* -1       # negative index for pure phases
                    rm_list = vcat(rm_list,idx);
                elseif i in db_inf.ss_name
                    idx = findfirst(db_inf.ss_name .== i)
                    rm_list = vcat(rm_list,idx);
                else
                    println(" \"$phase_type\" is not a proper phase type, and thus cannot be deactivated")
                end
            end
        end

        if isempty(rm_list)
            rm_list = nothing
            print(" The list of phases to be removed appears to be empty...")
        end
    else
        rm_list = nothing
    end
    
    return rm_list;
end

"""
    Dat = Initialize_MAGEMin(db = "ig"; verbose::Union{Bool, Int64} = true)

Initializes MAGEMin on one or more threads, for the database `db`. You can suppress all output with `verbose=false`. `verbose=true` will give a brief summary of the result, whereas `verbose=1` will give more details about the computations.
"""
function Initialize_MAGEMin(db = "ig";  verbose     ::Union{Int64,Bool} = 0,
                                        dataset     ::Union{Nothing,Int64} = nothing,
                                        limitCaOpx  ::Int64             = 0,
                                        CaOpxLim    ::Float64           = 0.0,
                                        mbCpx       ::Int64             = 1,
                                        mbIlm       ::Int64             = 0,
                                        mpSp        ::Int64             = 0,
                                        mpIlm       ::Int64             = 0,
                                        ig_ed       ::Int64             = 0,
                                        buffer      ::String            = "NONE",
                                        solver      ::Int64             = 0         )

    gv, z_b, DB, splx_data = init_MAGEMin(db;   verbose     = verbose,
                                                dataset     = dataset,
                                                mbCpx       = mbCpx,
                                                mbIlm       = mbIlm,
                                                mpSp        = mpSp,
                                                mpIlm       = mpIlm,
                                                ig_ed       = ig_ed,
                                                limitCaOpx  = limitCaOpx,
                                                CaOpxLim    = CaOpxLim,
                                                buffer      = buffer,
                                                solver      = solver    );

    nt              = Threads.maxthreadid()
    list_gv         = Vector{typeof(gv)}(undef, nt)
    list_z_b        = Vector{typeof(z_b)}(undef, nt)
    list_DB         = Vector{typeof(DB)}(undef, nt)
    list_splx_data  = Vector{typeof(splx_data)}(undef, nt)

    if isa(verbose,Bool)
        if verbose
            verbose=0
        else
            verbose=-1
        end
    end

    for id in 1:nt
        gv, z_b, DB, splx_data = init_MAGEMin(db;   verbose     = verbose,
                                                    dataset     = dataset,
                                                    mbCpx       = mbCpx,
                                                    mbIlm       = mbIlm,
                                                    mpSp        = mpSp,
                                                    mpIlm       = mpIlm,
                                                    ig_ed       = ig_ed,
                                                    limitCaOpx  = limitCaOpx,
                                                    CaOpxLim    = CaOpxLim,
                                                    buffer      = buffer,
                                                    solver      = solver    );

        list_gv[id]         = gv
        list_z_b[id]        = z_b
        list_DB[id]         = DB
        list_splx_data[id]  = splx_data
    end

    return MAGEMin_Data(db, list_gv, list_z_b, list_DB, list_splx_data)
end


"""
    Finalize_MAGEMin(dat::MAGEMin_Data)
Finalizes MAGEMin and clears variables
"""
function Finalize_MAGEMin(dat::MAGEMin_Data)
    for id in 1:Threads.maxthreadid()
        LibMAGEMin.FreeDatabases(dat.gv[id], dat.DB[id], dat.z_b[id])
        # splx_data needs to be freed
     end
     return nothing
end


"""
    gv, DB = init_MAGEMin(;EM_database=0)

Initializes MAGEMin (including setting global options) and loads the Database.
"""
function  init_MAGEMin( db          :: String               =  "ig";
                        verbose     :: Union{Int64,Bool}    =   0,
                        dataset     :: Union{Nothing,Int}   =   nothing,
                        mbCpx       :: Int64                =   0,
                        mbIlm       :: Int64                =   0,
                        mpSp        :: Int64                =   0,
                        mpIlm       :: Int64                =   0,
                        ig_ed       :: Int64                =   0,
                        limitCaOpx  :: Int64                =   0,
                        CaOpxLim    :: Float64              =   1.0,
                        buffer      :: String               =  "NONE",
                        solver      :: Int64                =   0           )

    z_b         = LibMAGEMin.bulk_infos()
    gv          = LibMAGEMin.global_variables()
    splx_data   = LibMAGEMin.simplex_data();
    DB          = LibMAGEMin.Database()
    gv          = LibMAGEMin.global_variable_alloc( pointer_from_objref(z_b))

    sb = ["sb11","sb21"]
    if db in sb
        rg = "sb"
    else
        rg = "tc"
    end

    if rg == "tc"
        if db == "mp"
            gv.EM_database = 0
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "mb"
            gv.EM_database = 1
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "mbe"
            gv.EM_database = 11
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "ig"
            gv.EM_database = 2
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "igad"
            gv.EM_database = 3
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "um"
            gv.EM_database = 4
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "ume"
            gv.EM_database = 5
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "mtl"
            gv.EM_database = 6
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "mpe"
            gv.EM_database = 7
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        else
            print("Database not implemented... using default mp\n")
        end
    elseif rg == "sb"
        unsafe_copyto!(convert(Ptr{UInt8}, gv.research_group), pointer(rg), length(rg) + 1)
        if db == "sb11"
            gv.EM_database = 0
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        elseif db == "sb21"
            gv.EM_database = 1
            unsafe_copyto!(convert(Ptr{UInt8}, gv.db), pointer(db), length(db) + 1)
        else 
            print("Database not implemented... using default sb11\n")
            gv.EM_database = 0
        end
    end

    gv.verbose      = verbose
    gv.mbCpx        = mbCpx
    gv.mbIlm        = mbIlm
    gv.mpSp         = mpSp
    gv.mpIlm        = mpIlm
    gv.ig_ed        = ig_ed

    gv.limitCaOpx   = limitCaOpx
    gv.CaOpxLim     = CaOpxLim
    gv.solver       = solver
    gv.buffer       = pointer(buffer)

    if !isnothing(dataset) && rg == "tc" && dataset in available_TC_ds
        gv.EM_dataset = dataset
    end

    gv              = LibMAGEMin.SetupDatabase(gv, pointer_from_objref(z_b))
    gv              = LibMAGEMin.global_variable_init(gv, pointer_from_objref(z_b))
    DB              = LibMAGEMin.InitializeDatabases(gv, gv.EM_database)

    LibMAGEMin.init_simplex_A(      pointer_from_objref(splx_data), gv)
    LibMAGEMin.init_simplex_B_em(   pointer_from_objref(splx_data), gv)

    return gv, z_b, DB, splx_data
end

"""
    finalize_MAGEMin(gv,DB)
Cleans up the memory
"""
function  finalize_MAGEMin(gv,DB, z_b)
    LibMAGEMin.FreeDatabases(gv, DB, z_b)
    return nothing
end


# wrapper for single point minimization
function single_point_minimization(     P           ::  T1,
                                        T           ::  T1,
                                        MAGEMin_db  ::  MAGEMin_Data;
                                        light       ::  Bool                            = false,
                                        name_solvus ::  Bool                            = false,
                                        fixed_bulk  ::  Bool                            = false,
                                        test        ::  Int64                           = 0, # if using a build-in test case,
                                        X           ::  VecOrMat                        = nothing,      
                                        B           ::  Union{Nothing, T1 }             = nothing,
                                        G           ::  Union{Nothing, Vector{LibMAGEMin.mSS_data},Vector{Vector{LibMAGEMin.mSS_data}}}  = nothing,
                                        scp         ::  Int64                           = 0,
                                        iguess      ::  Bool                            = false,
                                        rm_list     ::  Union{Nothing, Vector{Int64}}   = nothing,
                                        W           ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}}          = nothing,
                                        Xoxides     = Vector{String},
                                        sys_in      = "mol",
                                        rg          = "tc",
                                        progressbar = true        # show a progress bar or not?
                                        ) where {T1 <: Float64}

    P   = [P];
    T   = [T];
    if X isa AbstractVector{Float64}
        X = [X]
    end
    if !isnothing(B)
        B = [B]
    end

    Out_PT     =   multi_point_minimization(    P,
                                                T,
                                                MAGEMin_db;
                                                light       =   light,
                                                name_solvus =   name_solvus,
                                                fixed_bulk  =   fixed_bulk,
                                                test        =   test,
                                                X           =   X,
                                                B           =   B,
                                                G           =   G,   
                                                scp         =   scp,
                                                iguess      =   iguess,
                                                rm_list     =   rm_list,
                                                W           =   W,
                                                Xoxides     =   Xoxides,
                                                sys_in      =   sys_in,
                                                rg          =   rg,
                                                progressbar =   progressbar);
    return Out_PT[1]

end


"""
Out_PT = multi_point_minimization(P           ::  T2,
                                  T           ::  T2,
                                  MAGEMin_db  ::  MAGEMin_Data;
                                  light       ::  Bool                            = false,
                                  name_solvus ::  Bool                            = false,
                                  fixed_bulk  ::  Bool                            = false,
                                  test        ::  Int64                           = 0, # if using a build-in test case,
                                  X           ::  VecOrMat                        = nothing,
                                  B           ::  Union{Nothing, Vector{T1}}  = nothing,
                                  G           ::  Union{Nothing, Vector{LibMAGEMin.mSS_data},Vector{Vector{LibMAGEMin.mSS_data}}}  = nothing,
                                  scp         ::  Int64                           = 0, 
                                  iguess      ::  Union{Vector{Bool},Bool}        = false,
                                  rm_list     ::  Union{Nothing, Vector{Int64}}   = nothing,
                                  W           ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}}  = nothing,
                                  Xoxides     = Vector{String},
                                  sys_in      :: String                           = "mol",
                                  rg          :: String                           = "tc",
                                  progressbar :: Bool                             = true,        # show a progress bar or not?
                                  callback_fn ::  Union{Nothing, Function} = nothing, 
                                  callback_int::  Int64 = 1
                                  )


Perform (parallel) MAGEMin calculations for a range of points as a function of pressure `P`, temperature `T` and/or composition `X`. The database `MAGEMin_db` must be initialised before calling the routine.
The bulk-rock composition can either be set to be one of the pre-defined build-in test cases, or can be specified specifically by passing `X`, `Xoxides` and `sys_in` (that specifies whether the input is in "mol" or "wt").

Below a few examples:

Example 1 - build-in test vs. pressure and temperature
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> n = 10
julia> P = rand(8:40.0,n)
julia> T = rand(800:1500.0,n)
julia> out = multi_point_minimization(P, T, data, test=0)
julia> Finalize_MAGEMin(data)
```

Example 2 - Specify constant bulk rock composition for all points:
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> n = 10
julia> P = fill(10.0,n)
julia> T = fill(1100.0,n)
julia> Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> X = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> sys_in = "wt"
julia> out = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
julia> Finalize_MAGEMin(data)
```

Example 3 - Different bulk rock composition for different points
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> P = [10.0, 20.0]
julia> T = [1100.0, 1200]
julia> Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> X1 = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> X2 = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> X = [X1,X2]
julia> sys_in = "wt"
julia> out = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
julia> Finalize_MAGEMin(data)
```

Activating multithreading on julia
===

To take advantage of multithreading, you need to start julia from the terminal with:
```bash
\$ julia -t auto
```
which will automatically use all threads on your machine. Alternatively, use `julia -t 4` to start it on 4 threads.
If you are interested to see what you can do on your machine type:
```
julia> versioninfo()
```

"""
function multi_point_minimization(P           ::  T2,
                                  T           ::  T2,
                                  MAGEMin_db  ::  MAGEMin_Data;
                                  light       ::  Bool                            = false,
                                  name_solvus ::  Bool                            = false,
                                  fixed_bulk  ::  Bool                            = false,
                                  test        ::  Int64                           = 0, # if using a build-in test case,
                                  X           ::  VecOrMat                        = nothing,
                                  B           ::  Union{Nothing, Vector{T1}}  = nothing,
                                  G           ::  Union{Nothing, Vector{LibMAGEMin.mSS_data},Vector{Vector{LibMAGEMin.mSS_data}}}  = nothing,
                                  scp         ::  Int64                           = 0, 
                                  iguess      ::  Union{Vector{Bool},Bool}        = false,
                                  rm_list     ::  Union{Nothing, Vector{Int64}}   = nothing,
                                  W           ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}}  = nothing,
                                  Xoxides     = Vector{String},
                                  sys_in      :: String                           = "mol",
                                  rg          :: String                           = "tc",
                                  progressbar :: Bool                             = true,        # show a progress bar or not?
                                  callback_fn ::  Union{Nothing, Function} = nothing, 
                                  callback_int::  Int64 = 1
                                  ) where {T1 <: Float64, T2 <: AbstractVector{Float64}}

    # Set the compositional info
    CompositionType::Int64 = 0;

    if isnothing(X)
        # Use one of the build-in tests
        # Create thread-local data
        for i in 1:Threads.maxthreadid()
            MAGEMin_db.gv[i] = use_predefined_bulk_rock(MAGEMin_db.gv[i], test, MAGEMin_db.db);
        end
        CompositionType = 0;    # build-in tests
    else
        if isa(X,Vector{Float64})
        # same bulk rock composition for the full diagram
        @assert length(X) == length(Xoxides)

            # Set the bulk rock composition for all points
            for i in 1:Threads.maxthreadid()
                MAGEMin_db.gv[i] = define_bulk_rock(MAGEMin_db.gv[i], X, Xoxides, sys_in, MAGEMin_db.db);
            end
            CompositionType = 1;    # specified bulk composition for all points
        else
            @assert length(X) == length(P)
            CompositionType = 2;    # different bulk rock composition for every point
        end
    end

    # initialize vectors
    Out_PT = light ? Vector{light_gmin_struct{Float32, Int8}}(undef, length(P)) : Vector{gmin_struct{Float64, Int64}}(undef, length(P))
    # main loop
    if progressbar
        progr = Progress(length(P), desc="Computing $(length(P)) points...") # progress meter
    end
    count  = 0;
    @threads :static for i in eachindex(P)

        # Get thread-local buffers. As of Julia v1.9, a dynamic scheduling of
        # the threads is the default setting. To avoid task migration and the
        # resulting concurrency issues, we restrict the loop to static scheduling.
        id          = Threads.threadid()
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        if CompositionType == 2
            gv = define_bulk_rock(gv, X[i], Xoxides, sys_in, MAGEMin_db.db);
        end

        Gi          = isnothing(G) ? nothing :  G[i]
        ig          =  isa(iguess, Vector{Bool}) ? iguess[i] : iguess

        buffer      = isnothing(B) ? 0.0 :      B[i] 
        out         = point_wise_minimization(  P[i], T[i], gv, z_b, DB, splx_data;
                                                light=light, buffer_n=buffer, name_solvus=name_solvus, fixed_bulk=fixed_bulk, Gi=Gi, W=W, scp=scp, iguess=ig, rm_list=rm_list)

        Out_PT[i]   = deepcopy(out)

        if progressbar
            next!(progr)
        end
        if mod(i,callback_int)==0 && !isnothing(callback_fn)
            count   += 1
            callback_fn(count, length(P), time())
        end
    end
    if progressbar
        finish!(progr)
    end
    if !isnothing(callback_fn)
        callback_fn(length(P), length(P), time())
    end

    return Out_PT
end

"""
    AMR_minimization(   init_sub    ::  Int64,
                        ref_lvl     ::  Int64,
                        Prange      ::  Union{T1, NTuple{2, T1}},
                        Trange      ::  Union{T1, NTuple{2, T1}},
                        MAGEMin_db  ::  MAGEMin_Data;
                        test        ::  Int64                           = 0, # if using a build-in test case,
                        X           ::  VecOrMat                        = nothing,
                        B           ::  Union{Nothing, T1, Vector{T1}}  = 0.0,
                        scp         ::  Int64                           = 0,  
                        iguess      ::  Union{Vector{Bool},Bool}        = false,
                        rm_list     ::  Union{Nothing, Vector{Int64}}   = nothing,
                        W           ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}}  = nothing,
                        Xoxides     = Vector{String},
                        sys_in      :: String                          = "mol",
                        rg          :: String                          = "tc",
                        progressbar :: Bool                             = true        # show a progress bar or not?
                    ) where {T1 <: Float64}

    Performs an Adaptive Mesh Refinement (AMR) minimization for a range of points as a function of pressure `Prange`, temperature `Trange` and/or composition `X`. The database `MAGEMin_db` must be initialised before calling the routine.

# Example
```julia
data        = Initialize_MAGEMin("mp", verbose=-1, solver=0);

init_sub    =  1
ref_lvl     =  2
Prange      = (1.0,10.0)
Trange      = (400.0,800.0)
Xoxides     = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
X           = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
sys_in      = "mol"    
out         = AMR_minimization(init_sub, ref_lvl, Prange, Trange, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
```

"""
function AMR_minimization(  init_sub    ::  Int64,
                            ref_lvl     ::  Int64,
                            Prange      ::  Union{T1, NTuple{2, T1}},
                            Trange      ::  Union{T1, NTuple{2, T1}},
                            MAGEMin_db  ::  MAGEMin_Data;
                            test        ::  Int64                           = 0, # if using a build-in test case,
                            X           ::  VecOrMat                        = nothing,
                            B           ::  Union{Nothing, T1, Vector{T1}}  = 0.0,
                            scp         ::  Int64                           = 0,  
                            iguess      ::  Union{Vector{Bool},Bool}        = false,
                            rm_list     ::  Union{Nothing, Vector{Int64}}   = nothing,
                            W           ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}}  = nothing,
                            Xoxides     =  Vector{String},
                            sys_in      ::  String                          = "mol",
                            rg          ::  String                          = "tc",
                            progressbar :: Bool                             = true        # show a progress bar or not?
                            ) where {T1 <: Float64}

    Out_XY          =  Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,0)
    data            =  initialize_AMR(Trange,Prange,init_sub)
    Hash_XY         =  Vector{UInt64}(undef,0)
    for irefine = 1:ref_lvl+1
        if irefine > 1
            data    = split_and_keep(data, Hash_XY)
            data    = AMR(data)
        end

        if isempty(data.split_cell_list)
            Out_XY_new      = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,length(data.points))
            n_new_points    = length(data.points)
            npoints         = data.points
        else
            Out_XY_new      = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,length(data.npoints))
            n_new_points    = length(data.npoints)
            npoints         = data.npoints
        end

        if n_new_points > 0
            Tvec = zeros(Float64,n_new_points);
            Pvec = zeros(Float64,n_new_points);
            Xvec = Vector{Vector{Float64}}(undef,n_new_points);
            Bvec = zeros(Float64,n_new_points);
            if !isempty(data.split_cell_list) && iguess == true
                Gvec = Vector{Vector{LibMAGEMin.mSS_data}}(undef,n_new_points);
            else
                Gvec = nothing;
            end
            for i = 1:n_new_points
                Tvec[i] = npoints[i][1];
                Pvec[i] = npoints[i][2];
                Bvec[i] = B;
                Xvec[i] = X;
                if !isempty(data.split_cell_list) && iguess == true
                    tmp = [Out_XY[data.npoints_ig[i][j]].mSS_vec for j=1:length(data.npoints_ig[i])]
                    Gvec[i] = vcat(tmp...)
                end
            end
            Out_XY_new  =   multi_point_minimization(   Pvec, Tvec, MAGEMin_db,
                                                        X=Xvec, B=Bvec, G=Gvec, Xoxides=Xoxides, sys_in=sys_in, scp=scp, iguess=iguess, rm_list=rm_list, rg=rg, test=test); 
        else
            println("There is no new point to compute...")
        end
        Out_XY      = vcat(Out_XY, Out_XY_new)

        # Compute hash for all points
        n_points    = length(Out_XY)

        Hash_XY     = Vector{UInt64}(undef,n_points)
        for i=1:n_points
            Hash_XY[i]      = hash(sort(Out_XY[i].ph))
        end

    end

    return Out_XY
end


"""
bulk_rock = use_predefined_bulk_rock(gv, test=-1, db="ig")

Returns the pre-defined bulk rock composition of a given test
"""
function use_predefined_bulk_rock(gv, test=0, db="ig")

    rg = unsafe_string(gv.research_group)

    if rg=="tc"
        if      db == "mp"
            gv.test = test
            gv = LibMAGEMin.get_bulk_metapelite(gv)
        elseif  db == "mb"
            gv.test = test
            gv = LibMAGEMin.get_bulk_metabasite(gv)
        elseif  db == "mbe"
            gv.test = test
            gv = LibMAGEMin.get_bulk_metabasite(gv)
        elseif  db == "ig"
            gv.test = test
            gv = LibMAGEMin.get_bulk_igneous(gv)
        elseif  db == "igad"
            gv.test = test
            gv = LibMAGEMin.get_bulk_igneous_igad(gv)
        elseif  db == "um"
            gv.test = test
            gv = LibMAGEMin.get_bulk_ultramafic(gv)
        elseif  db == "ume"
            gv.test = test
            gv = LibMAGEMin.get_bulk_ultramafic_ext(gv)
        elseif  db == "mtl"
            gv.test = test
            gv = LibMAGEMin.get_bulk_mantle(gv)
        elseif  db == "mpe"
            gv.test = test
            gv = LibMAGEMin.get_bulk_metapelite_ext(gv)
        else
            print("Database not implemented...\n")
        end
    elseif rg=="sb"
        if  db == "sb11"
            gv.test = test
            gv = LibMAGEMin.get_bulk_stx11(gv)
        elseif  db == "sb21"
            gv.test = test
            gv = LibMAGEMin.get_bulk_stx21(gv)
        else
            print("Database not implemented...\n")
        end
    else
        print("Research group not implemented...\n")
    end

    LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)

    return gv
end

"""
    data = use_predefined_bulk_rock(data::MAGEMin_Data, test=0)

Returns the pre-defined bulk rock composition of a given test
"""
function use_predefined_bulk_rock(data::MAGEMin_Data, test=0)  
    nt = Threads.maxthreadid()
    for id in 1:nt
        data.gv[id] =  use_predefined_bulk_rock(data.gv[id], test, data.db);
    end
    return data
end

"""
    gv = define_bulk_rock(gv, bulk_in, bulk_in_ox, sys_in, db)

Defines the bulk-rock composition in `gv` structure using the input `bulk_in` and `bulk_in_ox`, converting it to the appropriate format for MAGEMin.
"""
function define_bulk_rock(gv, bulk_in, bulk_in_ox, sys_in,db)

    bulk_rock, ox   = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in,db)     # conversion changes the system unit to mol
    unsafe_copyto!(gv.bulk_rock, pointer(bulk_rock), gv.len_ox)            # copy the bulk-rock

    LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)

    return gv
end


function normalize(vector::AbstractVector{Float64})
    return vector ./ sum(vector)
end

"""
    bulk_mol = wt2mol(bulk_wt, bulk_ox)
Converts bulk-rock composition from wt to mol fraction
"""
function wt2mol(    bulk_wt     :: AbstractVector{Float64},
                    bulk_ox     :: AbstractVector{String}) 

    ref_ox          = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "MnO"; "H2O"; "CO2"; "S"];
    ref_MolarMass   = [60.08; 101.96; 56.08; 40.30; 71.85; 159.69; 94.2; 61.98; 79.88; 16.0; 151.99; 70.937; 18.015; 44.01; 32.06];      #Molar mass of oxides

    bulk_mol = zeros(length(bulk_ox));
    bulk_wt  = normalize(bulk_wt)

    for i = axes(bulk_ox,1)
        id = findfirst(ref_ox .== bulk_ox[i]);
        bulk_mol[i] = bulk_wt[i]/ref_MolarMass[id];
    end

    bulk_mol .= bulk_mol ./sum(bulk_mol) .* 100.0

    return bulk_mol
end


"""
    mol2wt(bulk_mol::Vector{Float64}, bulk_ox::Vector{String})

    Converts bulk-rock composition from mol to wt fraction
"""
function mol2wt(    bulk_mol     :: AbstractVector{Float64},
                    bulk_ox      :: AbstractVector{String}) 

    ref_ox          = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "MnO"; "H2O"; "CO2"; "S"];
    ref_MolarMass   = [60.08; 101.96; 56.08; 40.30; 71.85; 159.69; 94.2; 61.98; 79.88; 16.0; 151.99; 70.937; 18.015; 44.01; 32.06];      #Molar mass of oxides

    bulk_wt = zeros(length(bulk_ox));
    bulk_mol = normalize(bulk_mol)

    for i = axes(bulk_ox,1)
        id = findfirst(ref_ox .== bulk_ox[i]);
        bulk_wt[i] = bulk_mol[i]*ref_MolarMass[id];
    end

    bulk_wt .= bulk_wt ./sum(bulk_wt) .* 100.0

    return bulk_wt
end

"""
    MAGEMin_bulk, MAGEMin_ox; = convertBulk4MAGEMin(bulk_in::T1,bulk_in_ox::Vector{String},sys_in::String,db::String) where {T1 <: AbstractVector{Float64}}

receives bulk-rock composition in [mol,wt] fraction and associated oxide list and sends back bulk-rock composition converted for MAGEMin use

"""
function convertBulk4MAGEMin(   bulk_in     :: T1,
                                bulk_in_ox  :: Vector{String},
                                sys_in      :: String,
                                db          :: String ) where {T1 <: AbstractVector{Float64}}

    bulk_in = normalize(bulk_in);                            

	ref_ox          = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "MnO"; "H2O"; "CO2"; "S"];
	ref_MolarMass   = [60.08; 101.96; 56.08; 40.30; 71.85; 159.69; 94.2; 61.98; 79.88; 16.0; 151.99; 70.937; 18.015; 44.01; 32.06];      #Molar mass of oxides

    if db       == "mp"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"];
    elseif db   == "mb"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "H2O"];
    elseif db   == "mbe"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "H2O"];
    elseif db   == "ig"
	    MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    elseif db   == "igad"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"];
    elseif db   == "um"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "MgO"; "FeO"; "O"; "H2O"; "S"];
    elseif db   == "ume"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "MgO"; "FeO"; "O"; "H2O"; "S"; "CaO";"Na2O"];
    elseif db   == "mtl"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO";"Na2O"];
    elseif db   == "mpe"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"; "CO2"; "S"];
    elseif db   == "sb11"
        MAGEMin_ox      = ["SiO2"; "CaO"; "Al2O3"; "FeO"; "MgO";"Na2O"];
    elseif db   == "sb21"
        MAGEMin_ox      = ["SiO2"; "CaO"; "Al2O3"; "FeO"; "MgO";"Na2O"];
    else
        print("Database not implemented...\n")
    end

    # here we kick out oxides that are not part of the database
    filter = setdiff(bulk_in_ox,MAGEMin_ox)
    keep_ids = []
    if !isempty(filter)
        tmp_id    = findall(filter .!= "Fe2O3")
        if !isempty(tmp_id)
            filter = filter[tmp_id]
            
            for i=1:length(bulk_in_ox)
                if bulk_in_ox[i] in filter
                else
                    push!(keep_ids,i)
                end
            end

            bulk_in_ox  = bulk_in_ox[keep_ids]
            bulk_in     = bulk_in[keep_ids]
        end
    end

	MAGEMin_bulk    = zeros(length(MAGEMin_ox));
    bulk            = zeros(length(MAGEMin_ox));
    
	# convert to mol, if system unit = wt
	if sys_in == "wt"
		for i=1:length(bulk_in_ox)
            id = findall(ref_ox .== bulk_in_ox[i])[1];
			bulk[i] = bulk_in[i]/ref_MolarMass[id];
		end
    elseif sys_in == "mol"
		for i=1:length(bulk_in_ox)
			bulk[i] = bulk_in[i];
		end
    else
        println("System unit not implemented -> use 'mol' or 'wt' -> falling back to 'mol'")
	end

	bulk = normalize(bulk);

	for i=1:length(MAGEMin_ox)
        id = findall(bulk_in_ox .== MAGEMin_ox[i]);
		if ~isempty(id)
			MAGEMin_bulk[i] = bulk[id[1]];
		end
	end
    idFe2O3 = findall(bulk_in_ox .== "Fe2O3");

    if ~isempty(idFe2O3)
        idFeO = findall(MAGEMin_ox .== "FeO");
        MAGEMin_bulk[idFeO[1]] += bulk[idFe2O3[1]]*2.0;

        idO = findall(MAGEMin_ox .== "O");
        MAGEMin_bulk[idO[1]] += bulk[idFe2O3[1]];
    end

    MAGEMin_bulk .= normalize(MAGEMin_bulk);

    # check which component can safely be put to 0.0
    d = []
    c = collect(1:length(MAGEMin_ox))
    if db       == "mp"
        c = findall(MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O" .&& MAGEMin_ox .!= "MnO" .&& MAGEMin_ox .!= "H2O");
        d = findall(MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O" .|| MAGEMin_ox .!= "MnO");
    elseif db   == "mb"               #for the metabasite database it is better to set a low value for H2O as dry system haw not been validated by Eleanor
        c = findall(MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O");
        d = findall(MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O");
    elseif db   == "mbe"               #for the metabasite database it is better to set a low value for H2O as dry system haw not been validated by Eleanor
        c = findall(MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O");
        d = findall(MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O");
    elseif db   == "ig"
        c = findall(MAGEMin_ox .!= "K2O" .&& MAGEMin_ox .!= "Cr2O3" .&& MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O" .&& MAGEMin_ox .!= "H2O");
        d = findall(MAGEMin_ox .== "K2O" .|| MAGEMin_ox .== "Cr2O3" .|| MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O");# .|| MAGEMin_ox .== "H2O");
    elseif db   == "igad"
        c = findall(MAGEMin_ox .!= "Cr2O3" .&& MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O");
        d = findall(MAGEMin_ox .== "Cr2O3" .|| MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O");
    elseif db == "um"
        c = findall( MAGEMin_ox .!= "S" .&& MAGEMin_ox .!= "O" .&& MAGEMin_ox .!= "H2O");
        d = findall( MAGEMin_ox .== "S" .|| MAGEMin_ox .== "O");
    elseif db == "ume"
        c = findall( MAGEMin_ox .!= "S" .&& MAGEMin_ox .!= "O" .&& MAGEMin_ox .!= "H2O");
        d = findall( MAGEMin_ox .== "S" .|| MAGEMin_ox .== "O");
    elseif db == "mpe"
        c = findall(MAGEMin_ox .!= "CO2" .&& MAGEMin_ox .!= "S" .&& MAGEMin_ox .!= "TiO2" .&& MAGEMin_ox .!= "O" .&& MAGEMin_ox .!= "MnO" .&& MAGEMin_ox .!= "H2O");
        d = findall(MAGEMin_ox .== "CO2" .||MAGEMin_ox .== "S" .|| MAGEMin_ox .== "TiO2" .|| MAGEMin_ox .== "O" .|| MAGEMin_ox .!= "MnO");
    end

    id0 = findall(MAGEMin_bulk[c] .< 1e-4)
    if ~isempty(id0)
        MAGEMin_bulk[c[id0]] .= 1e-4;
    end

    id1 = findall(MAGEMin_bulk[d] .< 2e-5 .&& MAGEMin_bulk[d] .> -2e-5)
    if ~isempty(id1)
        MAGEMin_bulk[d[id1]] .= 0.0;
    end
    MAGEMin_bulk .= normalize(MAGEMin_bulk).*100.0

    return MAGEMin_bulk, MAGEMin_ox;
end


"""
    oint_wise_minimization(         P       ::Float64,
                                    T       ::Float64,
                                    gv,
                                    z_b,
                                    DB,
                                    splx_data;
                                    light       = false,
                                    name_solvus = false,
                                    fixed_bulk  = false,
                                    buffer_n    = 0.0,
                                    Gi          = nothing,
                                    scp         = 0,
                                    iguess      = false,
                                    rm_list     = nothing,
                                    W           = nothing   )

Computes the stable assemblage at `P` [kbar], `T` [C] and for a given bulk rock composition


# Example 1

This is an example of how to use it for a predefined bulk rock composition:
```julia
julia> db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
julia> gv, z_b, DB, splx_data      = init_MAGEMin(db);
julia> test        = 0;
julia> sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
julia> gv          = use_predefined_bulk_rock(gv, test, db)
julia> P           = 8.0;
julia> T           = 800.0;
julia> gv.verbose  = -1;        # switch off any verbose
julia> out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
Pressure          : 8.0      [kbar]
Temperature       : 800.0    [Celsius]
     Stable phase | Fraction (mol fraction)
              opx   0.24229
               ol   0.58808
              cpx   0.14165
              spl   0.02798
     Stable phase | Fraction (wt fraction)
              opx   0.23908
               ol   0.58673
              cpx   0.14583
              spl   0.02846
Gibbs free energy : -797.749183  (26 iterations; 94.95 ms)
Oxygen fugacity          : 9.645393319147175e-12
```

# Example 2
And here a case in which you specify your own bulk rock composition.
We convert that in the correct format, using the `convertBulk4MAGEMin` function.
```julia
julia> db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
julia> gv, z_b, DB, splx_data      = init_MAGEMin(db);
julia> bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> sys_in     = "wt"
julia> gv         = define_bulk_rock(gv, bulk_in, bulk_in_ox, sys_in, db);
julia> P,T         = 10.0, 1100.0;
julia> gv.verbose  = -1;        # switch off any verbose
julia> out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
Pressure          : 10.0      [kbar]
Temperature       : 1100.0    [Celsius]
     Stable phase | Fraction (mol fraction)
             pl4T   0.01114
              liq   0.74789
              cpx   0.21862
              opx   0.02154
     Stable phase | Fraction (wt fraction)
             pl4T   0.01168
              liq   0.72576
              cpx   0.23872
              opx   0.02277
Gibbs free energy : -907.27887  (47 iterations; 187.11 ms)
Oxygen fugacity          : 0.02411835177808492
julia> finalize_MAGEMin(gv,DB)
```

"""
function point_wise_minimization(   P       ::Float64,
                                    T       ::Float64,
                                    gv,
                                    z_b,
                                    DB,
                                    splx_data;
                                    light       = false,
                                    name_solvus = false,
                                    fixed_bulk  = false,
                                    buffer_n    = 0.0,
                                    Gi          = nothing,
                                    scp         = 0,
                                    iguess      = false,
                                    rm_list     = nothing,
                                    W           = nothing   )

    gv.buffer_n     =   buffer_n;
    input_data      =   LibMAGEMin.io_data();           # zero (not used actually)
    z_b.T           =   T + 273.15;                     # in K

    if P < 0.001
        P = 0.001
    end
    
    z_b.P           =   P
    gv.numPoint     =   1; 							    # the number of the current point */

    # Perform the point-wise minimization after resetting variables
    gv      = LibMAGEMin.reset_gv(gv,z_b, DB.PP_ref_db, DB.SS_ref_db)
    z_b     = LibMAGEMin.reset_z_b_bulk(	gv,	z_b	   )

    LibMAGEMin.reset_simplex_A(pointer_from_objref(splx_data), z_b, gv)
    LibMAGEMin.reset_simplex_B_em(pointer_from_objref(splx_data), gv)

    LibMAGEMin.reset_cp(gv,z_b, DB.cp)
    LibMAGEMin.reset_SS(gv,z_b, DB.SS_ref_db)
    LibMAGEMin.reset_sp(gv, DB.sp)

    if ~isnothing(rm_list)
        # if the list of phase to be removed is not empty then we first activate all combination
        # the following entries are set to 2, as only values of 0 or 1 are considered in the C code
        # i.e. that all phases are first set active before removing the ones in the list (with exception of the restricted ones)
        gv.mbCpx        = 2
        gv.mbIlm        = 2
        gv.mpSp         = 2
        gv.mpIlm        = 2
    end
    

    gv      = LibMAGEMin.ComputeG0_point(gv.EM_database, z_b, gv, DB.PP_ref_db,DB.SS_ref_db);

    # here we can over-ride default W's
    if ~isnothing(W)
        n_over  = length(W)
        Pw      = z_b.P
        Tw      = z_b.T
        for n=1:n_over

            if gv.EM_database == W[n].dtb
                ss          = W[n].ss_ids
                SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
                n_W         = SS_ref_db[ss].n_w;
                if W[n].n_Ws == n_W
                    # override = 1
                    # unsafe_copyto!(pointer(SS_ref_db[ss].override), pointer(override), 1) 
                    # SS_ref_db[ss].override = 1;                                                     # set the override flag to 1
                    # Wdef    = unsafe_wrap(Vector{Cdouble},SS_ref_db[ss].W, SS_ref_db[ss].n_w);      # retrieve default Ws
                    new_Ws  = zeros(n_W)
                    for i=1:n_W
                        new_Ws[i] = W[n].Ws[i,1] + W[n].Ws[i,2]*Tw + W[n].Ws[i,3]*Pw 
                    end
                    unsafe_copyto!(SS_ref_db[ss].W, pointer(new_Ws), SS_ref_db[ss].n_w) 
                else
                    print(" Wrong number of W's, please make sure the custom Ws are linked to the right solution model\n Ws override will be ignored\n")
                    println(" n_W target= $(n_W), n_W provided = $(W[n].n_Ws)")
                end
            end
            
        end
    end

    # gv      = LibMAGEMin.ComputeG0_point(gv.EM_database, z_b, gv, DB.PP_ref_db,DB.SS_ref_db);

    #= THIS IS WHERE pwm_init ends =#
    if ~isnothing(rm_list)

        SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
        pp_flags    = unsafe_wrap(Vector{Ptr{Int32}},gv.pp_flags, gv.len_pp);

        # here we manage a special case for igneous database
        if gv.EM_database == 2
            ss_names    = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss))
            pp_names    = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.PP_list, gv.len_pp))
            id_fl       = findfirst(ss_names .== "fl")
            if id_fl in rm_list
                id_H2O      = findfirst(pp_names .== "H2O")
                flags_on    = zeros(Int32,5); flags_on[1] = 1;
                unsafe_copyto!(pp_flags[id_H2O], pointer(flags_on), 5)
            end
        end

        flags_off = zeros(Int32,5);
        for i in rm_list
            if i > 0    # solution phase
                id = i
                unsafe_copyto!(SS_ref_db[id].ss_flags, pointer(flags_off), 5)
            else        # pure phase
                id = abs(i)
                unsafe_copyto!(pp_flags[id], pointer(flags_off), 5)
            end
        end
    end

    if iguess == true && Gi !== nothing
        SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
        PP_ref_db   = unsafe_wrap(Vector{LibMAGEMin.PP_ref},DB.PP_ref_db,gv.len_pp);
        n_SS_PC     = unsafe_wrap(Vector{Cint},gv.n_SS_PC,gv.len_ss);
        rg          = unsafe_string(gv.research_group)
        
        # retrieve dimensions
        np          = z_b.nzEl_val
        nzEl_array  = unsafe_wrap(Vector{Cint},z_b.nzEl_array, gv.len_ox) .+ 1
        nzEl_array  = nzEl_array[1:np]
    
        PC_read = Vector{LibMAGEMin.PC_type}(undef,gv.len_ss)
        if rg == "tc"
            LibMAGEMin.TC_PC_init(PC_read,gv)
        elseif rg == "sb"
            LibMAGEMin.SB_PC_init(PC_read,gv)
        end
    
        if fixed_bulk == true
            gv.fixed_bulk = 1
            # Declare array to be copied in splx_data
            A_jll       = zeros(np,np)
            g0_A_jll    = zeros(np)
            ph_id_A_jll = zeros(Int32,np,4)
            n_pc_ss     = zeros(gv.len_ss)

            # fill the arrays to be copied in splx_data
            for i = 1:np
                if Gi[i].ph_type == "pp"
                    ph_id = Gi[i].ph_id+1
                    g0_A_jll[i] = PP_ref_db[ph_id].gbase*PP_ref_db[ph_id].factor
                    A_jll[i,:]  = Gi[i].comp_Ppc[nzEl_array]

                    ph_id_A_jll[i,1] = 1
                    ph_id_A_jll[i,2] = ph_id-1
                    ph_id_A_jll[i,3] = 0
                    ph_id_A_jll[i,4] = 0
                elseif Gi[i].ph_type == "ss"
                    ph_id   = Gi[i].ph_id+1
                    ph      = Gi[i].ph_name

                    unsafe_copyto!(SS_ref_db[ph_id].gb_lvl,SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
                    unsafe_copyto!(SS_ref_db[ph_id].iguess,pointer(Gi[i].xeos_Ppc), SS_ref_db[ph_id].n_xeos)
                    if rg == "tc"
                        SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                    elseif rg == "sb"
                        SS_ref_db[ph_id] = LibMAGEMin.SB_PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                    end
                    g0_A_jll[i] = SS_ref_db[ph_id].df
                    A_jll[i,:]  = Gi[i].comp_Ppc[nzEl_array]
                    ph_id_A_jll[i,1] = 3
                    ph_id_A_jll[i,2] = ph_id-1
                    ph_id_A_jll[i,3] = 0
                    ph_id_A_jll[i,4] = n_pc_ss[ph_id]
                    n_pc_ss[ph_id]  += 1
                elseif Gi[i].ph_type == "ss_em"
                    ph_id   = Gi[i].ph_id+1
                    em_id   = Gi[i].em_id+1
                    ape     = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].ape, SS_ref_db[ph_id].n_em)
                    gbase   = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
                    comp_ptr= unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].Comp, SS_ref_db[ph_id].n_em)
                    Comp    = unsafe_wrap(Vector{Cdouble},comp_ptr[em_id], gv.len_ox)
                    factor 	= z_b.fbc/ape[em_id]
                    ph      = Gi[i].ph_name

                    g0_A_jll[i] = gbase[em_id]*factor;
                    A_jll[i,:]  = Comp[nzEl_array]*factor
                    ph_id_A_jll[i,1] = 2
                    ph_id_A_jll[i,2] = ph_id-1
                    ph_id_A_jll[i,3] = 0
                    ph_id_A_jll[i,4] = em_id-1
                end
            end
            # copy to the appropriate places
            ph_id_A = unsafe_wrap(Vector{Ptr{Int32}},splx_data.ph_id_A, np)

            for i=1:np
                unsafe_copyto!(ph_id_A[i],pointer(ph_id_A_jll[i,:]),4)
            end

            unsafe_copyto!(splx_data.A,pointer(vec(A_jll)),np*np)
            unsafe_copyto!(splx_data.A1,pointer(vec(A_jll)),np*np)
            unsafe_copyto!(splx_data.g0_A,pointer(g0_A_jll),np)
        else
            gv.fixed_bulk = 0
        end

        # add pseudocompounds
        n_mSS = length(Gi)
        for i = 1:n_mSS
    
            if Gi[i].ph_type == "ss" || Gi[i].ph_type == "ss_em"
                ph          = Gi[i].ph_name
                ph_id       = Gi[i].ph_id+1
                n_xeos      = SS_ref_db[ph_id].n_xeos
                n_em        = SS_ref_db[ph_id].n_em

                tot_pc      = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].tot_pc, 1)
                id_pc       = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].id_pc, 1)

                if tot_pc[1] < n_SS_PC[ph_id]   # here we make sure we have the space to store the pseudocompound

                    m_pc        = id_pc[1]+1;
                    info        = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].info,           SS_ref_db[ph_id].n_pc)
                    factor_pc   = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].factor_pc,   SS_ref_db[ph_id].n_pc)
                    DF_pc       = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].DF_pc,       SS_ref_db[ph_id].n_pc)
                    G_pc        = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].G_pc,        SS_ref_db[ph_id].n_pc)
        
                    ptr_comp_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].comp_pc,SS_ref_db[ph_id].n_pc)
                    ptr_p_pc    = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].p_pc,   SS_ref_db[ph_id].n_pc)
                    ptr_xeos_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].xeos_pc,SS_ref_db[ph_id].n_pc)
        
                    unsafe_copyto!(SS_ref_db[ph_id].gb_lvl,SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
                    xeos        = Gi[i].xeos_Ppc
        
                    # get solution phase information for given compositional variables
                    unsafe_copyto!(SS_ref_db[ph_id].iguess,pointer(xeos), n_xeos)
                    if rg == "tc"
                        SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                    elseif rg == "sb"
                        SS_ref_db[ph_id] = LibMAGEMin.SB_PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                    end
                    # copy solution phase composition
                    ss_comp     = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].ss_comp, gv.len_ox)
                    # println("ph: $ph comp: $ss_comp")
                    comp_pc     = unsafe_wrap(Vector{Cdouble}, ptr_comp_pc[m_pc], gv.len_ox)
                    comp_pc    .= ss_comp .* SS_ref_db[ph_id].factor;
        
                    # copy endmember fraction
                    p           = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].p, n_em)
                    p_pc        = unsafe_wrap(Vector{Cdouble}, ptr_p_pc[m_pc], n_em)
                    p_pc       .= p
        
                    # copy compositional variables
                    xeos_pc     = unsafe_wrap(Vector{Cdouble}, ptr_xeos_pc[m_pc], n_xeos)
                    xeos_pc    .= xeos
        
                    info[m_pc]      = 1;
                    factor_pc[m_pc] = SS_ref_db[ph_id].factor;
                    DF_pc[m_pc]     = SS_ref_db[ph_id].df;
                    G_pc[m_pc]      = SS_ref_db[ph_id].df;
        
                    tot_pc    .+= 1;
                    id_pc     .+= 1;
                end
            end
        end
    
        gv.leveling_mode = 1
    end

    # computing minimization
    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(gv.EM_database, input_data, z_b, gv, pointer_from_objref(splx_data),	DB.PP_ref_db, DB.SS_ref_db, DB.cp);

    # Postprocessing (NOTE: we should switch off printing if gv.verbose=0)
    gv = LibMAGEMin.ComputePostProcessing(z_b, gv, DB.PP_ref_db, DB.SS_ref_db, DB.cp)

    # Fill structure
    LibMAGEMin.fill_output_struct(gv, pointer_from_objref(splx_data),	z_b, DB.PP_ref_db, DB.SS_ref_db, DB.cp, DB.sp );

    # Print output to screen
    LibMAGEMin.PrintOutput(gv, 0, 1, DB, time, z_b);

    # Transform results to a more convenient julia struct
    if light == true
        out = deepcopy(create_light_gmin_struct(DB,gv));
    else    
        out = deepcopy(create_gmin_struct(DB, gv, time; name_solvus = name_solvus));
    end
    # here we compute specific heat capacity using reactions
    if (scp == 1)
        mSS_vec     = deepcopy(out.mSS_vec)
        dT          = 2.0;
        dP          = 0.002;
        out_W       = point_wise_minimization_with_guess(mSS_vec, P, T-dT, gv, z_b, DB, splx_data)
        out_E       = point_wise_minimization_with_guess(mSS_vec, P, T+dT, gv, z_b, DB, splx_data)

        hcp         = -(T+273.15)*(out_E.G_system + out_W.G_system - 2.0*out.G_system)/(dT*dT);
        # hcp         = (T+273.15)*(out_E.entropy - out.entropy)/(dT); # entropy way

        out_N       = point_wise_minimization_with_guess(mSS_vec, P+dP, T, gv, z_b, DB, splx_data)
        out_NE      = point_wise_minimization_with_guess(mSS_vec, P+dP, T+dT, gv, z_b, DB, splx_data)
        dGdT_N 		= (out_NE.G_system - out_N.G_system)	/(dT);
        dGdT_P 		= (out_E.G_system - out.G_system)	    /(dT);

        # out.entropy     .= -(out_E.entropy - out.entropy)/(dT);
        out.entropy     .= -(out_E.G_system - out.G_system)/(dT);
        out.enthalpy    .= out.entropy*(T+273.15) .+ out.G_system;
        out.s_cp   .= hcp/out.M_sys*1e6;
        out.alpha  .= 1.0/( (out_N.G_system - out.G_system)/dP * 10.0)*((dGdT_N-dGdT_P)/(dP))
    end

    return out
end



"""
    out = point_wise_minimization(P::Number,T::Number, data::MAGEMin_Data)

Performs a point-wise optimization for a given pressure `P` and temperature `T` for the data specified in the MAGEMin database `MAGEMin_Data` (where also compoition is specified)
"""
point_wise_minimization(P       ::  Number,
                        T       ::  Number,
                        gv,
                        z_b,
                        DB,
                        splx_data;
                        buffer_n::  Float64     = 0.0,
                        Gi      ::  Union{Nothing, Vector{LibMAGEMin.mSS_data}}  = nothing,
                        scp     ::  Int64       = 0,
                        iguess  ::  Bool        = false,
                        rm_list ::  Union{Nothing, Vector{Int64}}   = nothing,
                        name_solvus::Bool       = false,
                        fixed_bulk::Bool        = false,
                        W       ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}} = nothing) = 
                        point_wise_minimization(Float64(P),Float64(T), gv, z_b, DB, splx_data; buffer_n, Gi, scp, iguess, rm_list, name_solvus, fixed_bulk, W)

point_wise_minimization(P       ::  Number,
                        T       ::  Number,
                        gv      ::  LibMAGEMin.global_variables,
                        z_b     ::  LibMAGEMin.bulk_infos,
                        DB      ::  LibMAGEMin.Database,
                        splx_data:: LibMAGEMin.simplex_datas,
                        sys_in  ::  String;
                        buffer_n::  Float64     = 0.0,
                        Gi      ::  Union{Nothing, Vector{LibMAGEMin.mSS_data}}  = nothing,
                        scp     ::  Int64       = 0,
                        iguess  ::  Bool        = false,
                        rm_list ::  Union{Nothing, Vector{Int64}}   = nothing,
                        name_solvus::Bool       = false,
                        fixed_bulk::Bool        = false,
                        W       ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}} = nothing) = 
                        point_wise_minimization(Float64(P),Float64(T), gv, z_b, DB, splx_data; buffer_n, Gi, scp, iguess, rm_list, name_solvus, fixed_bulk, W)

point_wise_minimization(P       ::  Number,
                        T       ::  Number,
                        data    ::  MAGEMin_Data;
                        buffer_n::  Float64     = 0.0,
                        Gi      ::  Union{Nothing, Vector{LibMAGEMin.mSS_data}}  = nothing,
                        scp     ::  Int64       = 0,
                        iguess  ::  Bool        = false,
                        rm_list ::  Union{Nothing, Vector{Int64}}   = nothing,
                        name_solvus::Bool       = false,
                        fixed_bulk::Bool        = false,
                        W       ::  Union{Nothing, Vector{MAGEMin_C.W_data{Float64, Int64}}} = nothing) = 
                        point_wise_minimization(Float64(P),Float64(T), data.gv[1], data.z_b[1], data.DB[1], data.splx_data[1]; buffer_n, Gi, scp, iguess, rm_list, name_solvus, fixed_bulk, W)


"""
    pwm_init(P::Float64,T::Float64, gv, z_b, DB, splx_data; G0 = true)

    Function to provide single point minimization and give access to G0 and W's (Margules) parameters
    The objective here is to be able to use MAGEMin for thermodynamic database inversion/calibration

# Example

```julia
dtb     = "mp"
gv, z_b, DB, splx_data  = init_MAGEMin(dtb);

Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in  = "wt"
gv      = define_bulk_rock(gv, X, Xoxides, sys_in, dtb);
P, T    = 6.0, 500.0
gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);

# Modify gbase or W's here

out     = pwm_run(gv, z_b, DB, splx_data);
```
which gives:

```
julia> out       = pwm_run(gv, z_b, DB, splx_data);
Status             :            0 
Mass residual      : +3.23681e-15
Rank               :            0 
Point              :            1 
Temperature        :   +500.00000       [C] 
Pressure           :     +6.00000       [kbar]

SOL = [G: -791.477] (12 iterations, 11.03 ms)
GAM = [-943.275969,-1739.945965,-770.603948,-668.434121,-339.621452,-886.160861,-804.045787,-992.801306,-480.273909,-480.273909]

Phase :       bi      chl        g      opx      fsp        q      sph      H2O       zo 
Mode  :  0.18022  0.00161  0.00174  0.04239  0.33178  0.15155  0.01325  0.24445  0.03301 
```

"""
function pwm_init(P::Float64,T::Float64, gv, z_b, DB, splx_data; G0 = true)

    # input_data      =   LibMAGEMin.io_data();                           # zero (not used actually)

    z_b.T           =   T + 273.15    # in K
    z_b.P           =   P

    gv.numPoint     = 1; 							    # the number of the current point */

    # Perform the point-wise minimization after resetting variables
    gv      = LibMAGEMin.reset_gv(gv,z_b, DB.PP_ref_db, DB.SS_ref_db)
    z_b     = LibMAGEMin.reset_z_b_bulk(	gv,	z_b	   )

    LibMAGEMin.reset_simplex_A(pointer_from_objref(splx_data), z_b, gv)
    LibMAGEMin.reset_simplex_B_em(pointer_from_objref(splx_data), gv)

    LibMAGEMin.reset_cp(gv,z_b, DB.cp)
    LibMAGEMin.reset_SS(gv,z_b, DB.SS_ref_db)
    LibMAGEMin.reset_sp(gv, DB.sp)

    if G0 == true
        gv      = LibMAGEMin.ComputeG0_point(gv.EM_database, z_b, gv, DB.PP_ref_db,DB.SS_ref_db);
    end

    return (gv, z_b, DB, splx_data)
end

pwm_init(P::Number,T::Number, gv, z_b, DB, splx_data) = pwm_init(Float64(P),Float64(T), gv, z_b, DB, splx_data; G0 = true)


function pwm_run(gv, z_b, DB, splx_data; name_solvus = false)
    input_data      =   LibMAGEMin.io_data();                           # zero (not used actually)

    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(gv.EM_database, input_data, z_b, gv, pointer_from_objref(splx_data),	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

    # Postprocessing (NOTE: we should switch off printing if gv.verbose=0)
    gv = LibMAGEMin.ComputePostProcessing(z_b, gv, DB.PP_ref_db, DB.SS_ref_db, DB.cp)

    # Fill structure
    LibMAGEMin.fill_output_struct(gv, pointer_from_objref(splx_data), z_b, DB.PP_ref_db,DB.SS_ref_db, DB.cp, DB.sp );

    # Print output to screen
    LibMAGEMin.PrintOutput(gv, 0, 1, DB, time, z_b);

    # Transform results to a more convenient julia struct
    out = create_gmin_struct(DB, gv, time; name_solvus = name_solvus);

    # LibMAGEMin.FreeDatabases(gv, DB, z_b);

    return out
end


"""
    out = create_gmin_struct(gv, z_b)

This extracts the output of a pointwise MAGEMin optimization and adds it into a julia structure
"""
function create_gmin_struct(DB, gv, time; name_solvus = false)

    stb      = unsafe_load(DB.sp)

    MAGEMin_ver = unsafe_string(stb.MAGEMin_ver)
    dataset  = unsafe_string(stb.dataset)
    database = unsafe_string(stb.database)
    buffer   = unsafe_string(stb.buffer)
    buffer_n = stb.buffer_n
    G_system = stb.G
    Gamma    = unsafe_wrap(Vector{Cdouble},stb.gamma,gv.len_ox)
    P_kbar   = stb.P
    T_C      = stb.T-273.15
    X        = [stb.X]
    M_sys    = stb.M_sys

    # Bulk rock info (total, melt, solid, fluid)
    bulk     = unsafe_wrap(Vector{Cdouble},stb.bulk,   gv.len_ox)
    bulk_M   = unsafe_wrap(Vector{Cdouble},stb.bulk_M, gv.len_ox)
    bulk_S   = unsafe_wrap(Vector{Cdouble},stb.bulk_S, gv.len_ox)
    bulk_F   = unsafe_wrap(Vector{Cdouble},stb.bulk_F, gv.len_ox)

    # Bulk rock info (total, melt, solid, fluid)
    bulk_wt     = unsafe_wrap(Vector{Cdouble},stb.bulk_wt,   gv.len_ox)
    bulk_M_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_M_wt, gv.len_ox)
    bulk_S_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_S_wt, gv.len_ox)
    bulk_F_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_F_wt, gv.len_ox)

    # Solid, melt, fluid fractions
    frac_M   = stb.frac_M
    frac_S   = stb.frac_S
    frac_F   = stb.frac_F

    # Solid, melt, fluid fractions
    frac_M_wt   = stb.frac_M_wt
    frac_S_wt   = stb.frac_S_wt
    frac_F_wt   = stb.frac_F_wt

    frac_M_vol   = stb.frac_M_vol
    frac_S_vol   = stb.frac_S_vol
    frac_F_vol   = stb.frac_F_vol

    # Solid, melt, fluid densities
    alpha   = [stb.alpha]
    V       = stb.V
    # cp      = stb.cp
    s_cp    = [stb.s_cp]
    rho     = stb.rho
    rho_M   = stb.rho_M
    rho_S   = stb.rho_S
    rho_F   = stb.rho_F

    # Oxygen fugacity
    fO2     = stb.fO2
    dQFM    = stb.dQFM

    # Activities
    aH2O    = stb.aH2O
    aSiO2   = stb.aSiO2
    aTiO2   = stb.aTiO2
    aAl2O3  = stb.aAl2O3
    aMgO    = stb.aMgO
    aFeO    = stb.aFeO

    # thermodynamic properties
    entropy = [stb.entropy]
    enthalpy= [stb.enthalpy]

    # Stable assemblage info
    n_ph     =  stb.n_ph        # total # of stable phases
    n_PP     =  stb.n_PP        # number of pure phases
    n_SS     =  stb.n_SS        # number of solid solutions
    n_mSS    =  stb.n_mSS        # number of solid solutions

    ph_frac     =  unsafe_wrap(Vector{Cdouble},stb.ph_frac,         n_ph)
    ph_frac_wt  =  unsafe_wrap(Vector{Cdouble},stb.ph_frac_wt,      n_ph)
    ph_frac_1at =  unsafe_wrap(Vector{Cdouble},stb.ph_frac_1at,     n_ph)
    ph_frac_vol =  unsafe_wrap(Vector{Cdouble},stb.ph_frac_vol,     n_ph)
    ph_type     =  unsafe_wrap(Vector{Cint},   stb.ph_type,         n_ph)
    ph_id       =  unsafe_wrap(Vector{Cint},   stb.ph_id  ,         n_ph)
    ph_id_db    =  unsafe_wrap(Vector{Cint},   stb.ph_id_db ,       n_ph)

    # println("ph_frac $ph_frac ph_frac_wt $ph_frac_wt")
    ph          =  unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.ph, n_ph)) # stable phases
    sol_name    =  unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.sol_name, n_ph)) # stable phases

    # extract info about compositional variables of the solution models:
    SS_vec  = convert.(LibMAGEMin.SS_data, unsafe_wrap(Vector{LibMAGEMin.stb_SS_phase},stb.SS,n_SS))

    if name_solvus == true
        for i=1:n_SS
            ph[i] = get_mineral_name(database, ph[i], SS_vec[i])
        end
    end

    # create dictionaries for easy access to the phase symbols; has to be after the name_solvus is set!
    SS_syms = Dict( Symbol("$(ph[i])") => i for i=1:n_SS )
    PP_syms = Dict( Symbol("$(ph[i])") => i-n_SS for i=n_SS+1:n_SS+n_PP )


    # extract information about metastable solution phases
    mSS_vec = convert.(LibMAGEMin.mSS_data, unsafe_wrap(Vector{LibMAGEMin.mstb_SS_phase},stb.mSS,n_mSS))

    # Info about the endmembers:
    PP_vec  = convert.(LibMAGEMin.PP_data, unsafe_wrap(Vector{LibMAGEMin.stb_PP_phase},stb.PP,n_PP))

    # Names of oxides:
    oxides   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.oxides, gv.len_ox))
    elements = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.elements, gv.len_ox))

    # Numerics
    bulk_res_norm   =  gv.BR_norm
    iter            =  gv.global_ite
    time_ms         =  time*1000.0

    Vs_S = stb.Vs_S
    Vp_S = stb.Vp_S
    if isinf(Vp_S)
        Vp_S = NaN
    end
    if isinf(Vs_S)
        Vs_S = NaN
    end

    if frac_M > 0.0
        eta_M = compute_melt_viscosity_G08(oxides, bulk_M, T_C; A = -4.55)
    else
        eta_M = NaN
    end


    # Store all in output struct
    out = gmin_struct{Float64,Int64}( MAGEMin_ver, dataset, database, buffer, buffer_n, G_system, Gamma, P_kbar, T_C, X, M_sys,
                bulk, bulk_M, bulk_S, bulk_F,
                bulk_wt, bulk_M_wt, bulk_S_wt, bulk_F_wt,
                frac_M, frac_S, frac_F,
                frac_M_wt, frac_S_wt, frac_F_wt,
                frac_M_vol, frac_S_vol, frac_F_vol,
                alpha, V, #=cp,=# s_cp,
                rho, rho_M, rho_S, rho_F, eta_M,
                fO2, dQFM, aH2O, aSiO2, aTiO2, aAl2O3, aMgO, aFeO,
                n_PP, n_SS, n_mSS,
                ph_frac, ph_frac_wt, ph_frac_1at, ph_frac_vol, ph_type, ph_id, ph_id_db, ph, sol_name,
                SS_syms, PP_syms,
                SS_vec,  mSS_vec, PP_vec,
                oxides,  elements,
                stb.Vp, stb.Vs, Vp_S, Vs_S, stb.bulkMod, stb.shearMod, stb.bulkModulus_M,  stb.bulkModulus_S, stb.shearModulus_S,
                entropy, enthalpy,
                iter, bulk_res_norm, time_ms, stb.status)

   return out
end


"""
    out = create_light_gmin_struct(gv, z_b)

This extracts the output of a pointwise MAGEMin optimization and adds it into a julia structure
"""
function create_light_gmin_struct(DB,gv)

    stb         = unsafe_load(DB.sp)
    n_ph        =  stb.n_ph        # total # of stable phases
    n_PP        =  stb.n_PP        # number of pure phases
    n_SS        =  stb.n_SS        # number of solid solutions

    P_kbar      = Float32(stb.P)
    T_C         = Float32(stb.T-273.15)

    ph_frac_wt  =  Float32.(unsafe_wrap(Vector{Cdouble},  stb.ph_frac_wt,         n_ph))
    ph_type     =  Int8.(unsafe_wrap(Vector{Cint},        stb.ph_type,            n_ph))
    ph_id_db    =  Int8.(unsafe_wrap(Vector{Cint},        stb.ph_id_db,           n_ph))

    frac_F_wt   = Float32.(stb.frac_F_wt)
    frac_S_wt   = Float32.(stb.frac_S_wt)
    frac_M_wt   = Float32.(stb.frac_M_wt)

    bulk_S_wt   = Float32.(unsafe_wrap(Vector{Cdouble},stb.bulk_S_wt, gv.len_ox))
    bulk_F_wt   = Float32.(unsafe_wrap(Vector{Cdouble},stb.bulk_F_wt, gv.len_ox))
    bulk_M_wt   = Float32.(unsafe_wrap(Vector{Cdouble},stb.bulk_M_wt, gv.len_ox))

    rho_S       = Float32.(stb.rho_S)
    rho_F       = Float32.(stb.rho_F)
    rho_M       = Float32.(stb.rho_M)

    s_cp        = Float32.([stb.s_cp])
    alpha       = Float32.([stb.alpha])
    out = light_gmin_struct{Float32,Int8}(  P_kbar, T_C, ph_frac_wt, ph_type, ph_id_db,
                                            frac_S_wt, frac_F_wt, frac_M_wt,
                                            bulk_S_wt, bulk_F_wt, bulk_M_wt,
                                            rho_S, rho_F, rho_M,
                                            s_cp,alpha)

   return out
end

# Print brief info about pointwise calculation result
function show(io::IO, g::gmin_struct)
    println(io, "Pressure          : $(g.P_kbar)      [kbar]")
    println(io, "Temperature       : $(round(g.T_C,digits=4))    [Celsius]")

    println(io, "     Stable phase | Fraction (mol fraction) ")
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac[i], digits=5)) ")
    end
    println(io, "     Stable phase | Fraction (wt fraction) ")
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac_wt[i], digits=5)) ")
    end
    println(io, "     Stable phase | Fraction (vol fraction) ")
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac_vol[i], digits=5)) ")
    end
    println(io, "Gibbs free energy : $(round(g.G_system,digits=6))  ($(g.iter) iterations; $(round(g.time_ms,digits=2)) ms)")
    if g.status == 5
        println(io, "WARNING: calculation did not converge ----------------------------")
    end
    println(io, "Oxygen fugacity          : $(g.fO2)")
    println(io, "Delta QFM                : $(g.dQFM)")
end

"""
    print_info(g::gmin_struct)

Prints a more extensive overview of the simulation results
"""
function print_info(g::gmin_struct)

    println("Stable phases @ {$(round(g.P_kbar,digits=4)), $(round(g.T_C,digits=4))}  kbar/°C :")
    for i=1:length(g.ph)
        print("$(lpad(g.ph[i],6," ")) ")
    end
    print("  \n \n")

    # ==
    println("Compositional variables (solution phase):")
    for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")
        for j=1:length(g.SS_vec[i].compVariables)
            print("$(lpad(round(g.SS_vec[i].compVariables[j],digits=5),8," ")) ")
        end
        print("\n")
    end
    print("\n")
    # ==

    # ==
    println("End-members fraction [mol% 1 atom basis](solution phase):")
    for i=1:g.n_SS

        print("                ")
        for j=1:length(g.SS_vec[i].emNames)
            print("$(lpad(g.SS_vec[i].emNames[j],8," ")) ")
        end
        print("\n")
        print("$(lpad(g.ph[i],15," ")) ")
        for j=1:length(g.SS_vec[i].emFrac)
            print("$(lpad(round(g.SS_vec[i].emFrac[j],digits=5),8," ")) ")
        end
        print("\n")
    end
    print("\n")
    # ==
    println("End-members fraction [wt%](solution phase):")
    for i=1:g.n_SS

        print("                ")
        for j=1:length(g.SS_vec[i].emNames)
            print("$(lpad(g.SS_vec[i].emNames[j],8," ")) ")
        end
        print("\n")
        print("$(lpad(g.ph[i],15," ")) ")
        for j=1:length(g.SS_vec[i].emFrac_wt)
            print("$(lpad(round(g.SS_vec[i].emFrac_wt[j],digits=5),8," ")) ")
        end
        print("\n")
    end
    print("\n")
    # ==
      # ==
      println("Oxide compositions [mol% 1 atom basis] (normalized):")
      print("                ")
      for i=1:length(g.oxides)
          print("$(lpad(g.oxides[i],8," ")) ")
      end
      print("\n")

      print("$(lpad("SYS",15)) ")
      for i=1:length(g.oxides)
          print("$(lpad(round(g.bulk[i],digits=5),8," ")) ")
      end
      print("\n")

      for i=1:g.n_SS
          print("$(lpad(g.ph[i],15," ")) ")
          for j=1:length(g.oxides)
              print("$(lpad(round(g.SS_vec[i].Comp[j],digits=5),8," ")) ")
          end
          print("\n")
      end

      for i=1:g.n_PP
        print("$(lpad(g.ph[i+g.n_SS],15," ")) ")
        for j=1:length(g.oxides)
            print("$(lpad(round(g.PP_vec[i].Comp[j],digits=5),8," ")) ")
        end
        print("\n")
      end
      print("\n")

      # ==
    # ==

    # ==
      # ==
      println("Oxide compositions [wt%] (normalized):")
      print("                ")
      for i=1:length(g.oxides)
          print("$(lpad(g.oxides[i],8," ")) ")
      end
      print("\n")

      print("$(lpad("SYS",15)) ")
      for i=1:length(g.oxides)
          print("$(lpad(round(g.bulk_wt[i],digits=5),8," ")) ")
      end
      print("\n")

      for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")
        for j=1:length(g.oxides)
            print("$(lpad(round(g.SS_vec[i].Comp_wt[j],digits=5),8," ")) ")
        end
        print("\n")
    end

      print("\n")
      for i=1:g.n_PP
        print("$(lpad(g.ph[i+g.n_SS],15," ")) ")
        for j=1:length(g.oxides)
            print("$(lpad(round(g.PP_vec[i].Comp_wt[j],digits=5),8," ")) ")
        end
        print("\n")
      end
      print("\n")

      # ==
    # ==



    println("Stable mineral assemblage:")
    println("          phase  mode[mol1at] mode[wt]        f           G        V       Cp  rho[kg/m3]  Thermal_Exp Entropy[J/K] Enthalpy[J] BulkMod[GPa] ShearMod[GPa]   Vp[km/s]   Vs[km/s]")
    for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")
        print("$(lpad(round(g.ph_frac[i],digits=5),13," ")) ")
        print("$(lpad(round(g.ph_frac_wt[i],digits=5),8," ")) ")
        print("$(lpad(round(g.ph_frac_vol[i],digits=5),8," ")) ")
        print("$(lpad(round(g.SS_vec[i].f,digits=5),8," ")) ")
        print("$(lpad(round(g.SS_vec[i].G,digits=5),8," ")) ")
        print("$(lpad(round(g.SS_vec[i].V,digits=5),8," ")) ")
        print("$(lpad(round(g.SS_vec[i].cp,digits=5),8," ")) ")
        print("$(lpad(round(g.SS_vec[i].rho,digits=5),11," "))   ")
        print("$(lpad(round(g.SS_vec[i].alpha,digits=5),8," "))   ")
        print("$(lpad(round(g.SS_vec[i].entropy,digits=5),8," "))   ")
        print("$(lpad(round(g.SS_vec[i].enthalpy,digits=5),8," "))   ")
        print("$(lpad(round(g.SS_vec[i].bulkMod,digits=5),12," ")) ")
        print("$(lpad(round(g.SS_vec[i].shearMod,digits=5),13," ")) ")
        print("$(lpad(round(g.SS_vec[i].Vp,digits=5),10," ")) ")
        print("$(lpad(round(g.SS_vec[i].Vs,digits=5),10," ")) ")

        print("\n")
    end
    for i=1:g.n_PP
        print("$(lpad(g.ph[i+g.n_SS],15," ")) ")
        print("$(lpad(round(g.ph_frac[i],digits=5),13," ")) ")
        print("$(lpad(round(g.ph_frac_wt[i],digits=5),8," ")) ")
        print("$(lpad(round(g.ph_frac_vol[i],digits=5),8," ")) ")
        print("$(lpad(round(g.PP_vec[i].f,digits=5),8," ")) ")
        print("$(lpad(round(g.PP_vec[i].G,digits=5),8," ")) ")
        print("$(lpad(round(g.PP_vec[i].V,digits=5),8," ")) ")
        print("$(lpad(round(g.PP_vec[i].cp,digits=5),8," ")) ")
        print("$(lpad(round(g.PP_vec[i].rho,digits=5),11," "))   ")
        print("$(lpad(round(g.PP_vec[i].alpha,digits=5),8," "))   ")
        print("$(lpad(round(g.PP_vec[i].entropy,digits=5),8," "))   ")
        print("$(lpad(round(g.PP_vec[i].enthalpy,digits=5),8," "))   ")
        print("$(lpad(round(g.PP_vec[i].bulkMod,digits=5),12," ")) ")
        print("$(lpad(round(g.PP_vec[i].shearMod,digits=5),13," ")) ")
        print("$(lpad(round(g.PP_vec[i].Vp,digits=5),10," ")) ")
        print("$(lpad(round(g.PP_vec[i].Vs,digits=5),10," ")) ")

        print("\n")
    end

    print("$(lpad("SYS",15," ")) ")
    print("$(lpad(round(sum(g.ph_frac),digits=5),13," ")) ")
    print("$(lpad(round(sum(g.ph_frac_wt),digits=5),8," ")) ")
    print("$(lpad(round(sum(g.ph_frac_vol),digits=5),8," ")) ")
    print("$(lpad(round(g.G_system,digits=5),20," ")) ")
    print("$(lpad(round(g.alpha[1],digits=5),29," ")) ")
    print("$(lpad(round(g.cp,digits=5),29," ")) ")
    print("$(lpad(round(g.s_cp[1],digits=5),29," ")) ")
    print("$(lpad(round(g.V,digits=5),29," ")) ")
    print("$(lpad(round(g.rho,digits=5),29," ")) ")
    print("$(lpad(round(g.entropy,digits=5),21," ")) ")
    print("$(lpad(round(g.enthalpy,digits=5),13," ")) ")
    print("$(lpad(round(g.bulkMod,digits=5),13," ")) ")
    print("$(lpad(round(g.shearMod,digits=5),13," ")) ")
    print("$(lpad(round(g.Vp,digits=5),10," ")) ")
    print("$(lpad(round(g.Vs,digits=5),10," ")) ")
    print("\n")
    print("\n")
    # ==

    println("Gamma (chemical potential of oxides):")
    for i=1:length(g.oxides)
        println("  $(lpad(g.oxides[i],6," "))  $(rpad(round(g.Gamma[i],digits=5),15," ")) ")
    end
    print("\n")

    println("delta Gibbs energy (G-hyperplane distance):")
    for i=1:g.n_SS
        println("  $(lpad(g.ph[i],6," "))  $(rpad(g.SS_vec[i].deltaG,15," ")) ")
    end
    for i=1:g.n_PP
        println("  $(lpad(g.ph[i+g.n_SS],6," "))  $(rpad(g.PP_vec[i].deltaG,15," ")) ")
    end

    print("\n")

    println("mass residual :  $(g.bulk_res_norm)")
    println("# iterations  :  $(g.iter)")
    println("status        :  $(g.status)")

    print("\n")

end


"""
    out = function point_wise_minimization_with_guess(mSS_vec :: Vector{mSS_data}, P, T, gv, z_b, DB, splx_data)
"""
function point_wise_minimization_with_guess(    mSS_vec ::  Vector{LibMAGEMin.mSS_data},
                                                P       ::  Float64,
                                                T       ::  Float64,
                                                gv      ::  LibMAGEMin.global_variables,
                                                z_b     ::  LibMAGEMin.bulk_infos,
                                                DB      ::  LibMAGEMin.Database,
                                                splx_data:: LibMAGEMin.simplex_datas)

    # initialize MAGEMin up to G0 computation included
    gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);
    gv.verbose = -1


    # println("reasearch group: $rg")
    ############################################################################
    # retrieve Solution Phases information
    rg          = unsafe_string(gv.research_group)
    SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    # retrieve dimensions
    np          = z_b.nzEl_val
    nzEl_array  = unsafe_wrap(Vector{Cint},z_b.nzEl_array, gv.len_ox) .+ 1
    nzEl_array  = nzEl_array[1:np]

    PC_read = Vector{LibMAGEMin.PC_type}(undef,gv.len_ss)
    if rg == "tc"
        LibMAGEMin.TC_PC_init(PC_read,gv)
    elseif rg == "sb"
        LibMAGEMin.SB_PC_init(PC_read,gv)
    end

    # add pseudocompounds
    n_mSS = length(mSS_vec)
    for i = 1:n_mSS

        if mSS_vec[i].ph_type == "ss"
            ph          = mSS_vec[i].ph_name
            ph_id       = mSS_vec[i].ph_id+1
            n_xeos      = SS_ref_db[ph_id].n_xeos
            n_em        = SS_ref_db[ph_id].n_em

            tot_pc      = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].tot_pc, 1)
            id_pc       = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].id_pc, 1)
            info        = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].info, SS_ref_db[ph_id].n_pc)
            factor_pc   = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].factor_pc, SS_ref_db[ph_id].n_pc)
            DF_pc       = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].DF_pc, SS_ref_db[ph_id].n_pc)
            G_pc        = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].G_pc, SS_ref_db[ph_id].n_pc)

            # here we make sure the pseudocompound can be added
            if id_pc[1]+1 < SS_ref_db[ph_id].n_pc 

                m_pc        = id_pc[1]+1;
                ptr_comp_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].comp_pc,SS_ref_db[ph_id].n_pc)
                ptr_p_pc    = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].p_pc,SS_ref_db[ph_id].n_pc)
                ptr_xeos_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].xeos_pc,SS_ref_db[ph_id].n_pc)

                unsafe_copyto!(SS_ref_db[ph_id].gb_lvl,SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
                xeos        = mSS_vec[i].xeos_Ppc

                # retrieve bounds
                bounds_ref      = zeros( n_xeos,2)
                ptr_bounds_ref  = unsafe_wrap(Vector{Ptr{Cdouble}}, SS_ref_db[ph_id].bounds_ref, n_xeos)

                for k=1:n_xeos
                    bounds_ref[k,:] = unsafe_wrap(Vector{Cdouble}, ptr_bounds_ref[k], 2)
                    if xeos[k] < bounds_ref[k,1]
                        xeos[k] = bounds_ref[k,1]
                    elseif xeos[k] > bounds_ref[k,2]
                        xeos[k] = bounds_ref[k,2]
                    end
                end

                # get solution phase information for given compositional variables
                unsafe_copyto!(SS_ref_db[ph_id].iguess,pointer(xeos), n_xeos)
                if rg == "tc"
                    SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                elseif rg == "sb"
                    SS_ref_db[ph_id] = LibMAGEMin.SB_PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
                end

                # copy solution phase composition
                ss_comp     = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].ss_comp, gv.len_ox)
  
                # println("ph $ph n_pc $(SS_ref_db[ph_id].n_pc) m_pc $m_pc")
                comp_pc     = unsafe_wrap(Vector{Cdouble}, ptr_comp_pc[m_pc], gv.len_ox)
                comp_pc    .= ss_comp .* SS_ref_db[ph_id].factor;

                # copy endmember fraction
                p           = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].p, n_em)
                p_pc        = unsafe_wrap(Vector{Cdouble}, ptr_p_pc[m_pc], n_em)
                p_pc       .= p

                # copy compositional variables
                xeos_pc     = unsafe_wrap(Vector{Cdouble}, ptr_xeos_pc[m_pc], n_xeos)
                xeos_pc    .= xeos

                info[m_pc]      = 1;
                factor_pc[m_pc] = SS_ref_db[ph_id].factor;
                DF_pc[m_pc]     = SS_ref_db[ph_id].df;
                G_pc[m_pc]      = SS_ref_db[ph_id].df;

                tot_pc .+= 1;
                id_pc  .+= 1;
            end
        end
    end

    gv.leveling_mode = 1
    out = deepcopy(pwm_run(gv, z_b, DB, splx_data))

    return out
end



"""
    out = point_wise_metastability(out, P, T, data)

    where out is a MAGEMin_C output structure of type MAGEMin_C.gmin_struct{Float64, Int64}, P is pressure, T is temperature and data is the MAGEMin_C working data structure. 

    This function computes the metastability of the solution phases in `out` at pressure `P` and temperature `T`.

# Example:

```julia
using MAGEMin_C

data    = Initialize_MAGEMin("mp", verbose=-1; solver=0);
P,T     = 6.0, 630.0
Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in  = "wt"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

Pmeta, Tmeta       = 6.0, 500.0
out2    = point_wise_metastability(out, Pmeta, Tmeta, data)
```

which gives:

````
julia> out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
Pressure          : 6.0      [kbar]
Temperature       : 630.0    [Celsius]
     Stable phase | Fraction (mol fraction) 
               bi   0.18975 
              fsp   0.33769 
              opx   0.06615 
                q   0.14587 
              sph   0.01103 
              H2O   0.24952 
     Stable phase | Fraction (wt fraction) 
               bi   0.21347 
              fsp   0.44152 
              opx   0.07259 
                q   0.17078 
              sph   0.01404 
              H2O   0.0876 
     Stable phase | Fraction (vol fraction) 
               bi   0.16535 
              fsp   0.38301 
              opx   0.05008 
                q   0.1507 
              sph   0.00921 
              H2O   0.24165 
Gibbs free energy : -806.707117  (23 iterations; 15.5 ms)
Oxygen fugacity          : 11.422305637291458
Delta QFM                : 29.789994195381436

julia> Pmeta, Tmeta       = 6.0, 500.0
(6.0, 500.0)

julia> out2    = point_wise_metastability(out, Pmeta, Tmeta, data)
Pressure          : 6.0      [kbar]
Temperature       : 500.0    [Celsius]
     Stable phase | Fraction (mol fraction) 
               bi   0.18975 
              fsp   0.33769 
              opx   0.06615 
                q   0.14587 
              sph   0.01103 
              H2O   0.24952 
     Stable phase | Fraction (wt fraction) 
               bi   0.21347 
              fsp   0.44152 
              opx   0.07259 
                q   0.17078 
              sph   0.01404 
              H2O   0.0876 
     Stable phase | Fraction (vol fraction) 
               bi   0.16819 
              fsp   0.39062 
              opx   0.05102 
                q   0.1528 
              sph   0.0094 
              H2O   0.22795 
Gibbs free energy : -791.460287  (0 iterations; 0.04 ms)
Oxygen fugacity          : 11.256905434522107
Delta QFM                : 34.26950061699225
````

"""
function point_wise_metastability(  out     :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    P       :: Float64,
                                    T       :: Float64,
                                    gv, z_b, DB, splx_data)

    mSS_vec = deepcopy(out.mSS_vec)                                
    gv      = define_bulk_rock(gv, out.bulk, out.oxides, "mol", out.database);
    # initialize MAGEMin up to G0 computation included
    gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);
    gv.verbose = -1

    ############################################################################
    # retrieve Pure Phases information
    PP_ref_db   = unsafe_wrap(Vector{LibMAGEMin.PP_ref},DB.PP_ref_db,gv.len_pp);

    # retrieve Solution Phases information
    SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    # retrieve dimensions
    np          = z_b.nzEl_val
    nzEl_array  = unsafe_wrap(Vector{Cint},z_b.nzEl_array, gv.len_ox) .+ 1
    nzEl_array  = nzEl_array[1:np]

    # Declare array to be copied in splx_data
    A_jll       = zeros(np,np)
    g0_A_jll    = zeros(np)
    ph_id_A_jll = zeros(Int32,np,4)

    n_pc_ss     = zeros(gv.len_ss)

    rg          = unsafe_string(gv.research_group)

    PC_read = Vector{LibMAGEMin.PC_type}(undef,gv.len_ss)
    if rg == "tc"
        LibMAGEMin.TC_PC_init(PC_read,gv)
    elseif rg == "sb"
        LibMAGEMin.SB_PC_init(PC_read,gv)
    end
    
    # fill the arrays to be copied in splx_data
    for i = 1:np
        if mSS_vec[i].ph_type == "pp"
            ph_id = mSS_vec[i].ph_id+1
            g0_A_jll[i] = PP_ref_db[ph_id].gbase*PP_ref_db[ph_id].factor
            A_jll[i,:]  = mSS_vec[i].comp_Ppc[nzEl_array]

            ph_id_A_jll[i,1] = 1
            ph_id_A_jll[i,2] = ph_id-1
            ph_id_A_jll[i,3] = 0
            ph_id_A_jll[i,4] = 0
        elseif mSS_vec[i].ph_type == "ss"
            ph_id   = mSS_vec[i].ph_id+1
            ph      = mSS_vec[i].ph_name

            unsafe_copyto!(SS_ref_db[ph_id].gb_lvl,SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
            unsafe_copyto!(SS_ref_db[ph_id].iguess,pointer(mSS_vec[i].xeos_Ppc), SS_ref_db[ph_id].n_xeos)
            SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)

            g0_A_jll[i] = SS_ref_db[ph_id].df
            A_jll[i,:]  = mSS_vec[i].comp_Ppc[nzEl_array]
            ph_id_A_jll[i,1] = 3
            ph_id_A_jll[i,2] = ph_id-1
            ph_id_A_jll[i,3] = 0
            ph_id_A_jll[i,4] = n_pc_ss[ph_id]
            n_pc_ss[ph_id]  += 1
        elseif mSS_vec[i].ph_type == "ss_em"
            ph_id   = mSS_vec[i].ph_id+1
            em_id   = mSS_vec[i].em_id+1
            ape     = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].ape, SS_ref_db[ph_id].n_em)
            gbase   = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
            comp_ptr= unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].Comp, SS_ref_db[ph_id].n_em)
            Comp    = unsafe_wrap(Vector{Cdouble},comp_ptr[em_id], gv.len_ox)
            factor 	= z_b.fbc/ape[em_id]
            ph      = mSS_vec[i].ph_name

            g0_A_jll[i] = gbase[em_id]*factor;
            A_jll[i,:]  = Comp[nzEl_array]*factor
            ph_id_A_jll[i,1] = 2
            ph_id_A_jll[i,2] = ph_id-1
            ph_id_A_jll[i,3] = 0
            ph_id_A_jll[i,4] = em_id-1
        end
    end

    # copy to the appropriate places
    ph_id_A = unsafe_wrap(Vector{Ptr{Int32}},splx_data.ph_id_A, np)

    for i=1:np
        unsafe_copyto!(ph_id_A[i],pointer(ph_id_A_jll[i,:]),4)
    end

    unsafe_copyto!(splx_data.A,pointer(vec(A_jll)),np*np)
    unsafe_copyto!(splx_data.A1,pointer(vec(A_jll)),np*np)
    unsafe_copyto!(splx_data.g0_A,pointer(g0_A_jll),np)

    # add pseudocompounds
    n_mSS = length(mSS_vec)
    for i = 1:n_mSS

        if mSS_vec[i].ph_type == "ss"
            ph          = mSS_vec[i].ph_name
            ph_id       = mSS_vec[i].ph_id+1
            n_xeos      = SS_ref_db[ph_id].n_xeos
            n_em        = SS_ref_db[ph_id].n_em

            tot_pc      = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].tot_pc, 1)
            id_pc       = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].id_pc, 1)
            info        = unsafe_wrap(Vector{Cint},SS_ref_db[ph_id].info, gv.max_n_mSS)
            factor_pc   = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].factor_pc, gv.max_n_mSS)
            DF_pc       = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].DF_pc, gv.max_n_mSS)
            G_pc        = unsafe_wrap(Vector{Cdouble},SS_ref_db[ph_id].G_pc, gv.max_n_mSS)

            m_pc        = id_pc[1]+1;
            ptr_comp_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].comp_pc,gv.max_n_mSS)
            ptr_p_pc    = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].p_pc,gv.max_n_mSS)
            ptr_xeos_pc = unsafe_wrap(Vector{Ptr{Cdouble}},SS_ref_db[ph_id].xeos_pc,gv.max_n_mSS)

            unsafe_copyto!(SS_ref_db[ph_id].gb_lvl,SS_ref_db[ph_id].gbase, SS_ref_db[ph_id].n_em)
            xeos        = mSS_vec[i].xeos_Ppc

            # retrieve bounds
            bounds_ref      = zeros( n_xeos,2)
            ptr_bounds_ref  = unsafe_wrap(Vector{Ptr{Cdouble}}, SS_ref_db[ph_id].bounds_ref, n_xeos)

            for k=1:n_xeos
                bounds_ref[k,:] = unsafe_wrap(Vector{Cdouble}, ptr_bounds_ref[k], 2)
                if xeos[k] < bounds_ref[k,1]
                    xeos[k] = bounds_ref[k,1]
                elseif xeos[k] > bounds_ref[k,2]
                    xeos[k] = bounds_ref[k,2]
                end
            end

            # get solution phase information for given compositional variables
            unsafe_copyto!(SS_ref_db[ph_id].iguess,pointer(xeos), n_xeos)
            SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)

            # copy solution phase composition
            ss_comp     = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].ss_comp, gv.len_ox)
            comp_pc     = unsafe_wrap(Vector{Cdouble}, ptr_comp_pc[m_pc], gv.len_ox)
            comp_pc    .= ss_comp .* SS_ref_db[ph_id].factor;

            # copy endmember fraction
            p           = unsafe_wrap(Vector{Cdouble}, SS_ref_db[ph_id].p, n_em)
            p_pc        = unsafe_wrap(Vector{Cdouble}, ptr_p_pc[m_pc], n_em)
            p_pc       .= p

            # copy compositional variables
            xeos_pc     = unsafe_wrap(Vector{Cdouble}, ptr_xeos_pc[m_pc], n_xeos)
            xeos_pc    .= xeos

            info[m_pc]      = 1;
            factor_pc[m_pc] = SS_ref_db[ph_id].factor;
            DF_pc[m_pc]     = SS_ref_db[ph_id].df;
            G_pc[m_pc]      = SS_ref_db[ph_id].df;

            tot_pc .+= 1;
            id_pc  .+= 1;
        end
    end

    gv.leveling_mode    = 2
    gv.solver           = 3    #this deactivates minimization 
    out                 = deepcopy(pwm_run(gv, z_b, DB, splx_data))

    return out
end

point_wise_metastability(           out     :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    P       :: Float64,
                                    T       :: Float64,
                                    data    ::  MAGEMin_Data) =
                                    point_wise_metastability(out, P, T, data.gv[1], data.z_b[1], data.DB[1], data.splx_data[1])



# The following section add post-processing routines
include("TE_ph_models.jl")
include("TE_partitioning.jl")
include("Zircon_saturation.jl")
include("Sulfur_saturation.jl")
include("export2CSV.jl")
include("External_routines.jl")

# Loading Adaptive mesh refinement functions
include("AMR.jl")



const out_TE_struct = out_tepm
const out_struct    = gmin_struct{Float64, Int64}

