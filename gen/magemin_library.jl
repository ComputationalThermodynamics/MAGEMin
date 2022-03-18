module LibMAGEMin

using MAGEMin_jll
export MAGEMin_jll

using CEnum

@cenum var"##Ctag#316"::UInt32 begin
    _tc_ds634_ = 0
end

struct EM_db
    Name::NTuple{20, Cchar}
    Comp::NTuple{12, Cdouble}
    input_1::NTuple{3, Cdouble}
    input_2::NTuple{4, Cdouble}
    input_3::NTuple{11, Cdouble}
end

function Access_EM_DB(id, EM_database)
    ccall((:Access_EM_DB, libMAGEMin), EM_db, (Cint, Cint), id, EM_database)
end

struct EM2id
    EM_tag::NTuple{20, Cchar}
    id::Cint
    hh::Cint
end

function find_EM_id(EM_tag)
    ccall((:find_EM_id, libMAGEMin), Cint, (Ptr{Cchar},), EM_tag)
end

struct PP2id
    PP_tag::NTuple{20, Cchar}
    id::Cint
    hh::Cint
end

function find_PP_id(PP_tag)
    ccall((:find_PP_id, libMAGEMin), Cint, (Ptr{Cchar},), PP_tag)
end

function get_EM_DB_names(EM_database)
    ccall((:get_EM_DB_names, libMAGEMin), Ptr{Ptr{Cchar}}, (Cint,), EM_database)
end

struct bulk_info
    P::Cdouble
    T::Cdouble
    R::Cdouble
    bulk_rock::Ptr{Cdouble}
    nzEl_val::Cint
    zEl_val::Cint
    nzEl_array::Ptr{Cint}
    zEl_array::Ptr{Cint}
    apo::Ptr{Cdouble}
    fbc::Cdouble
    masspo::Ptr{Cdouble}
end

function zeros_in_bulk(bulk, P, T)
    ccall((:zeros_in_bulk, libMAGEMin), bulk_info, (Ptr{Cdouble}, Cdouble, Cdouble), bulk, P, T)
end

# no prototype is found for this function at Initialize.h:143:17, please use with caution
function global_variable_init()
    ccall((:global_variable_init, libMAGEMin), Cint, ())
end

function norm_array(array, size)
    ccall((:norm_array, libMAGEMin), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), array, size)
end

function get_bulk(bulk_rock, test, n_El)
    ccall((:get_bulk, libMAGEMin), Cvoid, (Ptr{Cdouble}, Cint, Cint), bulk_rock, test, n_El)
end

function runMAGEMin(argc, argv)
    ccall((:runMAGEMin, libMAGEMin), Cint, (Cint, Ptr{Ptr{Cchar}}), argc, argv)
end

struct SS_refs
    P::Cdouble
    T::Cdouble
    R::Cdouble
    EM_list::Ptr{Ptr{Cchar}}
    ss_n::Cdouble
    delta_ss_n::Cdouble
    ss_flags::Ptr{Cint}
    CstFactor::Cint
    n_pc::Cint
    tot_pc::Cint
    id_pc::Cint
    n_swap::Ptr{Cint}
    info::Ptr{Cint}
    G_pc::Ptr{Cdouble}
    DF_pc::Ptr{Cdouble}
    comp_pc::Ptr{Ptr{Cdouble}}
    p_pc::Ptr{Ptr{Cdouble}}
    mu_pc::Ptr{Ptr{Cdouble}}
    xeos_pc::Ptr{Ptr{Cdouble}}
    factor_pc::Ptr{Cdouble}
    ub_pc::Ptr{Cdouble}
    lb_pc::Ptr{Cdouble}
    solvus_id::Ptr{Cint}
    min_mode::Cint
    is_liq::Cint
    symmetry::Cint
    n_em::Cint
    n_xeos::Cint
    n_sf::Cint
    n_w::Cint
    eye::Ptr{Ptr{Cdouble}}
    W::Ptr{Cdouble}
    v::Ptr{Cdouble}
    sum_v::Cdouble
    n_v::Cint
    sf_ok::Cint
    Comp::Ptr{Ptr{Cdouble}}
    gbase::Ptr{Cdouble}
    mu_array::Ptr{Ptr{Cdouble}}
    gb_lvl::Ptr{Cdouble}
    factor::Cdouble
    box_bounds::Ptr{Ptr{Cdouble}}
    box_bounds_default::Ptr{Ptr{Cdouble}}
    z_em::Ptr{Cdouble}
    n_guess::Cint
    iguess::Ptr{Cdouble}
    dguess::Ptr{Cdouble}
    nz_em::Cint
    check_df::Cdouble
    forced_stop::Cint
    xeos_sf_ok_saved::Cint
    status::Cint
    nlopt_verb::Cint
    tol_sf::Ptr{Cdouble}
    lb::Ptr{Cdouble}
    ub::Ptr{Cdouble}
    opt::Cint
    fbc::Cdouble
    sum_apep::Cdouble
    p::Ptr{Cdouble}
    ape::Ptr{Cdouble}
    mat_phi::Ptr{Cdouble}
    mu_Gex::Ptr{Cdouble}
    sf::Ptr{Cdouble}
    dsf::Ptr{Cdouble}
    mu::Ptr{Cdouble}
    dfx::Ptr{Cdouble}
    dp_dx::Ptr{Ptr{Cdouble}}
    df::Cdouble
    df_raw::Cdouble
    LM_time::Cdouble
    ss_comp::Ptr{Cdouble}
    xi_em::Ptr{Cdouble}
    sum_xi::Cdouble
    xeos::Ptr{Cdouble}
    xeos_sf_ok::Ptr{Cdouble}
    density::Ptr{Cdouble}
    phase_density::Cdouble
    volume::Cdouble
    mass::Cdouble
end

const SS_ref = SS_refs

struct IODATA
    n_phase::Cint
    P::Cdouble
    T::Cdouble
    bulk::Ptr{Cdouble}
    in_gam::Ptr{Cdouble}
    phase_names::Ptr{Ptr{Cchar}}
    phase_xeos::Ptr{Ptr{Cdouble}}
    phase_emp::Ptr{Ptr{Cdouble}}
end

const io_data = IODATA

struct OUTDATA
    n_phase::Cint
    P::Cdouble
    T::Cdouble
    G_system::Cdouble
    BR_norm::Cdouble
    dG_norm::Cdouble
    number::Cint
    status::Cint
    Gamma::Ptr{Cdouble}
    n_SS::Cint
    n_PP::Cint
    StableSolutions::Ptr{Ptr{Cchar}}
    StableFractions::Ptr{Cdouble}
    Phasedensity::Ptr{Cdouble}
    max_num_EM::Cint
    n_em::Ptr{Cint}
    xEOS::Ptr{Ptr{Cdouble}}
    p_EM::Ptr{Ptr{Cdouble}}
end

const out_data = OUTDATA

struct csd_phase_sets
    name::Ptr{Cchar}
    split::Cint
    in_iter::Cint
    id::Cint
    n_xeos::Cint
    n_em::Cint
    n_sf::Cint
    sf_ok::Cint
    ss_flags::Ptr{Cint}
    ss_n::Cdouble
    delta_ss_n::Cdouble
    df::Cdouble
    factor::Cdouble
    min_time::Cdouble
    sum_xi::Cdouble
    sum_dxi::Cdouble
    z_em::Ptr{Cdouble}
    p_em::Ptr{Cdouble}
    xi_em::Ptr{Cdouble}
    lvlxeos::Ptr{Cdouble}
    dguess::Ptr{Cdouble}
    xeos::Ptr{Cdouble}
    dpdx::Ptr{Ptr{Cdouble}}
    dfx::Ptr{Cdouble}
    mu::Ptr{Cdouble}
    mu0::Ptr{Cdouble}
    delta_mu::Ptr{Cdouble}
    sf::Ptr{Cdouble}
    ss_comp::Ptr{Cdouble}
    gbase::Ptr{Cdouble}
    mass::Cdouble
    volume::Cdouble
    phase_density::Cdouble
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_shearModulus::Cdouble
end

const csd_phase_set = csd_phase_sets

struct global_variables
    version::Ptr{Cchar}
    verbose::Cint
    outpath::Ptr{Cchar}
    Mode::Cint
    numDiff::Ptr{Ptr{Cdouble}}
    n_Diff::Cint
    relax_PGE::Cdouble
    relax_PGE_val::Cdouble
    PC_df_add::Cdouble
    PC_min_dist::Cdouble
    PC_check_val::Cdouble
    check_PC::Cint
    check_PC_ite::Cint
    act_varFac_stab::Cdouble
    len_pp::Cint
    len_ss::Cint
    len_ox::Cint
    max_n_cp::Cint
    len_cp::Cint
    ox::Ptr{Ptr{Cchar}}
    gam_tot::Ptr{Cdouble}
    del_gam_tot::Ptr{Cdouble}
    delta_gam_tot::Ptr{Cdouble}
    n_flags::Cint
    PP_list::Ptr{Ptr{Cchar}}
    SS_list::Ptr{Ptr{Cchar}}
    pp_n::Ptr{Cdouble}
    pp_xi::Ptr{Cdouble}
    delta_pp_n::Ptr{Cdouble}
    delta_pp_xi::Ptr{Cdouble}
    pp_flags::Ptr{Ptr{Cint}}
    numPoint::Cint
    global_ite::Cint
    LVL_time::Cdouble
    em2ss_shift::Cdouble
    bnd_filter_pc::Cdouble
    n_pc::Cint
    max_G_pc::Cdouble
    n_SS_PC::Ptr{Cint}
    SS_PC_stp::Ptr{Cdouble}
    eps_sf_pc::Cdouble
    verifyPC::Ptr{Cint}
    n_solvi::Ptr{Cint}
    id_solvi::Ptr{Ptr{Cint}}
    ineq_res::Cdouble
    obj_tol::Cdouble
    newly_added::Ptr{Cint}
    box_size_mode_1::Cdouble
    maxeval::Cint
    maxeval_mode_1::Cint
    bnd_val::Cdouble
    init_prop::Ptr{Cdouble}
    A_PGE::Ptr{Cdouble}
    b_PGE::Ptr{Cdouble}
    dn_cp::Ptr{Cdouble}
    dn_pp::Ptr{Cdouble}
    cp_id::Ptr{Cint}
    pp_id::Ptr{Cint}
    fc_norm_t1::Cdouble
    outter_PGE_ite::Cint
    inner_PGE_ite::Cint
    inner_PGE_ite_time::Cdouble
    n_phase::Cint
    n_pp_phase::Cint
    n_cp_phase::Cint
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    br_liq_x::Cdouble
    max_fac::Cdouble
    max_br::Cdouble
    br_max_rlx::Cdouble
    ur_1::Cint
    ur_2::Cint
    ur_3::Cint
    ur_f::Cint
    div::Cint
    dGamma::Ptr{Cdouble}
    PGE_mass_norm::Ptr{Cdouble}
    PGE_total_norm::Ptr{Cdouble}
    gamma_norm::Ptr{Cdouble}
    ite_time::Ptr{Cdouble}
    G_system::Cdouble
    G_system_mu::Cdouble
    br_norm::Cdouble
    br_max_tol::Cdouble
    alpha::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    remove_dG_val::Cdouble
    remove_sum_xi::Cdouble
    ph_change::Cint
    A::Ptr{Ptr{Cdouble}}
    b::Ptr{Cdouble}
    mass_residual::Ptr{Cdouble}
    BR_norm::Cdouble
    gb_P_eps::Cdouble
    gb_T_eps::Cdouble
end

const global_variable = global_variables

struct PP_refs
    Name::NTuple{20, Cchar}
    Comp::NTuple{11, Cdouble}
    gbase::Cdouble
    gb_lvl::Cdouble
    factor::Cdouble
    phase_density::Cdouble
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_shearModulus::Cdouble
    volume::Cdouble
    mass::Cdouble
end

const PP_ref = PP_refs

struct Database
    PP_ref_db::Ptr{PP_ref}
    SS_ref_db::Ptr{SS_ref}
    cp::Ptr{csd_phase_set}
    EM_names::Ptr{Ptr{Cchar}}
end

const Databases = Database

function InitializeDatabases(gv, EM_database)
    ccall((:InitializeDatabases, libMAGEMin), Databases, (global_variable, Cint), gv, EM_database)
end

function FreeDatabases(gv, DB)
    ccall((:FreeDatabases, libMAGEMin), Cvoid, (global_variable, Databases), gv, DB)
end

function ComputeEquilibrium_Point(EM_database, input_data, Mode, z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:ComputeEquilibrium_Point, libMAGEMin), global_variable, (Cint, io_data, Cint, bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, input_data, Mode, z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function ComputePostProcessing(EM_database, z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:ComputePostProcessing, libMAGEMin), Cvoid, (Cint, bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function ReadCommandLineOptions(gv, argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, n_pc_out, maxeval_out, get_version_out)
    ccall((:ReadCommandLineOptions, libMAGEMin), global_variable, (global_variable, Cint, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), gv, argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, n_pc_out, maxeval_out, get_version_out)
end

function PrintOutput(gv, rank, l, DB, time_taken, z_b)
    ccall((:PrintOutput, libMAGEMin), Cvoid, (global_variable, Cint, Cint, Databases, Cdouble, bulk_info), gv, rank, l, DB, time_taken, z_b)
end

# typedef void ( * sf_type ) ( unsigned m , double * result , unsigned n , const double * x , double * grad , void * data )
const sf_type = Ptr{Cvoid}

function NLopt_global_opt_function()
    ccall((:NLopt_global_opt_function, libMAGEMin), Cint, ())
end

function NLopt_opt_function()
    ccall((:NLopt_opt_function, libMAGEMin), Cint, ())
end

function PGE()
    ccall((:PGE, libMAGEMin), Cint, ())
end

function norm_vector(array, n)
    ccall((:norm_vector, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, n)
end

struct ss_pc
    xeos_pc::NTuple{11, Cdouble}
end

struct PC_refs
    ss_pc_xeos::Ptr{ss_pc}
end

const PC_ref = PC_refs

function SS_PC_init_function(SS_PC_xeos, iss, name)
    ccall((:SS_PC_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_PC_xeos, iss, name)
end

function dump_init(gv)
    ccall((:dump_init, libMAGEMin), Cvoid, (Cint,), gv)
end

function dump_results_function(gv, z_b, PP_ref_db, SS_ref_db, cp)
    ccall((:dump_results_function, libMAGEMin), Cvoid, (Cint, bulk_info, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), gv, z_b, PP_ref_db, SS_ref_db, cp)
end

function mergeParallelFiles(gv)
    ccall((:mergeParallelFiles, libMAGEMin), Cvoid, (Cint,), gv)
end

function mergeParallel_LocalMinima_Files(gv)
    ccall((:mergeParallel_LocalMinima_Files, libMAGEMin), Cvoid, (Cint,), gv)
end

function mergeParallel_LevellingGamma_Files(gv)
    ccall((:mergeParallel_LevellingGamma_Files, libMAGEMin), Cvoid, (Cint,), gv)
end

function sum_array(array, size)
    ccall((:sum_array, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, size)
end

function check_sign(v1, v2)
    ccall((:check_sign, libMAGEMin), Cint, (Cdouble, Cdouble), v1, v2)
end

function G_EM_function(EM_database, bulk_rock, P, T, name, state)
    ccall((:G_EM_function, libMAGEMin), PP_ref, (Cint, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), EM_database, bulk_rock, P, T, name, state)
end

function G_SS_EM_function()
    ccall((:G_SS_EM_function, libMAGEMin), Cint, ())
end

function G_SS_INIT_EM_function()
    ccall((:G_SS_INIT_EM_function, libMAGEMin), Cint, ())
end

function CP_INIT_function()
    ccall((:CP_INIT_function, libMAGEMin), Cint, ())
end

function read_in_data(gv, input_data, file_name, n_points)
    ccall((:read_in_data, libMAGEMin), Cvoid, (Cint, Ptr{Cint}, Ptr{Cchar}, Cint), gv, input_data, file_name, n_points)
end

function AddResults_output_struct(gv, z_b, P, T, DB, output)
    ccall((:AddResults_output_struct, libMAGEMin), Cvoid, (Cint, bulk_info, Cdouble, Cdouble, Cint, Cint), gv, z_b, P, T, DB, output)
end

function InitializeOutput()
    ccall((:InitializeOutput, libMAGEMin), Cint, ())
end

function FreeOutput(output)
    ccall((:FreeOutput, libMAGEMin), Cvoid, (Cint,), output)
end

struct ketopt_t
    ind::Cint
    opt::Cint
    arg::Ptr{Cchar}
    longidx::Cint
    i::Cint
    pos::Cint
    n_args::Cint
end

struct ko_longopt_t
    name::Ptr{Cchar}
    has_arg::Cint
    val::Cint
end

function ketopt_permute(argv, j, n)
    ccall((:ketopt_permute, libMAGEMin), Cvoid, (Ptr{Ptr{Cchar}}, Cint, Cint), argv, j, n)
end

function ketopt(s, argc, argv, permute, ostr, longopts)
    ccall((:ketopt, libMAGEMin), Cint, (Ptr{ketopt_t}, Cint, Ptr{Ptr{Cchar}}, Cint, Ptr{Cchar}, Ptr{ko_longopt_t}), s, argc, argv, permute, ostr, longopts)
end

# typedef double ( * obj_type ) ( unsigned n , const double * x , double * grad , void * SS_ref_db )
const obj_type = Ptr{Cvoid}

function SS_objective_init_function(SS_objective, gv)
    ccall((:SS_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, Cint), SS_objective, gv)
end

function p2x_bi(SS_ref_db, eps)
    ccall((:p2x_bi, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_cd(SS_ref_db, eps)
    ccall((:p2x_cd, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_cpx(SS_ref_db, eps)
    ccall((:p2x_cpx, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_ep(SS_ref_db, eps)
    ccall((:p2x_ep, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_fl(SS_ref_db, eps)
    ccall((:p2x_fl, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_g(SS_ref_db, eps)
    ccall((:p2x_g, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_hb(SS_ref_db, eps)
    ccall((:p2x_hb, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_ilm(SS_ref_db, eps)
    ccall((:p2x_ilm, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_liq(SS_ref_db, eps)
    ccall((:p2x_liq, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_mu(SS_ref_db, eps)
    ccall((:p2x_mu, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_ol(SS_ref_db, eps)
    ccall((:p2x_ol, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_opx(SS_ref_db, eps)
    ccall((:p2x_opx, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_pl4T(SS_ref_db, eps)
    ccall((:p2x_pl4T, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function p2x_spn(SS_ref_db, eps)
    ccall((:p2x_spn, libMAGEMin), Cvoid, (Cint, Cdouble), SS_ref_db, eps)
end

function obj_bi(n, x, grad, SS_ref_db)
    ccall((:obj_bi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_cd(n, x, grad, SS_ref_db)
    ccall((:obj_cd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_cpx(n, x, grad, SS_ref_db)
    ccall((:obj_cpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ep(n, x, grad, SS_ref_db)
    ccall((:obj_ep, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_fl(n, x, grad, SS_ref_db)
    ccall((:obj_fl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_g(n, x, grad, SS_ref_db)
    ccall((:obj_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_hb(n, x, grad, SS_ref_db)
    ccall((:obj_hb, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_liq(n, x, grad, SS_ref_db)
    ccall((:obj_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mu(n, x, grad, SS_ref_db)
    ccall((:obj_mu, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ol(n, x, grad, SS_ref_db)
    ccall((:obj_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_opx(n, x, grad, SS_ref_db)
    ccall((:obj_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_pl4T(n, x, grad, SS_ref_db)
    ccall((:obj_pl4T, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_spn(n, x, grad, SS_ref_db)
    ccall((:obj_spn, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function PC_PX_function()
    ccall((:PC_PX_function, libMAGEMin), Cint, ())
end

function PC_function()
    ccall((:PC_function, libMAGEMin), Cint, ())
end

function P2X()
    ccall((:P2X, libMAGEMin), Cint, ())
end

function get_phase_id(gv, name)
    ccall((:get_phase_id, libMAGEMin), Cint, (Cint, Ptr{Cchar}), gv, name)
end

function phase_update_function()
    ccall((:phase_update_function, libMAGEMin), Cint, ())
end

function phase_merge_function()
    ccall((:phase_merge_function, libMAGEMin), Cint, ())
end

struct str
    value::Cdouble
    index::Cint
end

function cmp_dbl(a, b)
    ccall((:cmp_dbl, libMAGEMin), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), a, b)
end

function cmp_int(a, b)
    ccall((:cmp_int, libMAGEMin), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), a, b)
end

function pp_min_function(gv, z_b, PP_ref_db)
    ccall((:pp_min_function, libMAGEMin), Cvoid, (Cint, bulk_info, Ptr{Cint}), gv, z_b, PP_ref_db)
end

function init_em_db()
    ccall((:init_em_db, libMAGEMin), Cint, ())
end

function inverseMatrix(A1, n)
    ccall((:inverseMatrix, libMAGEMin), Cvoid, (Ptr{Cdouble}, Cint), A1, n)
end

function VecMatMul(B1, A1, B, n)
    ccall((:VecMatMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), B1, A1, B, n)
end

function MatVecMul(A1, br, n_vec, n)
    ccall((:MatVecMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), A1, br, n_vec, n)
end

function Levelling()
    ccall((:Levelling, libMAGEMin), Cint, ())
end

struct simplex_datas
    gamma_ps::Ptr{Cdouble}
    gamma_ss::Ptr{Cdouble}
    gamma_tot::Ptr{Cdouble}
    gamma_delta::Ptr{Cdouble}
    min_F::Cdouble
    ph2swp::Cint
    n_swp::Cint
    swp::Cint
    pivot::Ptr{Cint}
    A::Ptr{Cdouble}
    A1::Ptr{Cdouble}
    ph_id_A::Ptr{Ptr{Cint}}
    g0_A::Ptr{Cdouble}
    dG_A::Ptr{Cdouble}
    n_vec::Ptr{Cdouble}
    n_Ox::Cint
    len_ox::Cint
    n_pp::Cint
    n_em_ss::Cint
    B::Ptr{Cdouble}
    B1::Ptr{Cdouble}
    ph_id_B::Ptr{Cint}
    g0_B::Cdouble
    dG_B::Cdouble
    n_B::Cint
    n_local_min::Cint
    n_filter::Cint
end

const simplex_data = simplex_datas

function print_levelling(z_b, gv, PP_ref_db, SS_ref_db)
    ccall((:print_levelling, libMAGEMin), Cvoid, (bulk_info, Cint, Ptr{Cint}, Ptr{Cint}), z_b, gv, PP_ref_db, SS_ref_db)
end

function SS_UPDATE_function()
    ccall((:SS_UPDATE_function, libMAGEMin), Cint, ())
end

function CP_UPDATE_function()
    ccall((:CP_UPDATE_function, libMAGEMin), Cint, ())
end

function split_cp()
    ccall((:split_cp, libMAGEMin), Cint, ())
end

function ss_min_PGE(mode, i, gv, z_b, SS_ref_db, cp)
    ccall((:ss_min_PGE, libMAGEMin), Cvoid, (Cint, Cint, Cint, bulk_info, Ptr{Cint}, Ptr{Cint}), mode, i, gv, z_b, SS_ref_db, cp)
end

function reset_global_variables()
    ccall((:reset_global_variables, libMAGEMin), Cint, ())
end

function init_ss_db()
    ccall((:init_ss_db, libMAGEMin), Cint, ())
end

function reset_phases()
    ccall((:reset_phases, libMAGEMin), Cint, ())
end

function SS_ref_destroy(gv, SS_ref_db)
    ccall((:SS_ref_destroy, libMAGEMin), Cvoid, (Cint, Ptr{Cint}), gv, SS_ref_db)
end

function CP_destroy(gv, cp)
    ccall((:CP_destroy, libMAGEMin), Cvoid, (Cint, Ptr{Cint}), gv, cp)
end

function _DCDCT_fct(id, result, A, n_act_sf, n_xeos)
    ccall((:_DCDCT_fct, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Cint, Cint), id, result, A, n_act_sf, n_xeos)
end

function _DC_Null_fct(id, result, A, B, n_xeos, n_act_sf)
    ccall((:_DC_Null_fct, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Ptr{Cdouble}, Cint, Cint), id, result, A, B, n_xeos, n_act_sf)
end

function _Epsilon_C_fct(id, result, A, b, n_xeos, n_sf)
    ccall((:_Epsilon_C_fct, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint), id, result, A, b, n_xeos, n_sf)
end

function _Epsilon_J_fct(result, A, b, n_xeos)
    ccall((:_Epsilon_J_fct, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), result, A, b, n_xeos)
end

function _I_DC_Null_fct(id, result, A, B, eye, n_act_sf, n_xeos)
    ccall((:_I_DC_Null_fct, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Ptr{Cdouble}, Cint, Cint), id, result, A, B, eye, n_act_sf, n_xeos)
end

function _FillEyeMatrix(A, n)
    ccall((:_FillEyeMatrix, libMAGEMin), Cvoid, (Ptr{Cdouble}, Cint), A, n)
end

function get_act_sf_id(result, A, n)
    ccall((:get_act_sf_id, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Cint), result, A, n)
end

function MatMatMul(A, nrowA, B, ncolB, common, C)
    ccall((:MatMatMul, libMAGEMin), Cvoid, (Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Cint, Ptr{Ptr{Cdouble}}), A, nrowA, B, ncolB, common, C)
end

function pseudo_inverse(matrix, B, m, n)
    ccall((:pseudo_inverse, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint), matrix, B, m, n)
end

function get_act_sf(A, n)
    ccall((:get_act_sf, libMAGEMin), Cint, (Ptr{Cdouble}, Cint), A, n)
end

function get_active_em(array, n)
    ccall((:get_active_em, libMAGEMin), Cint, (Ptr{Cdouble}, Cint), array, n)
end

function EndsWithTail(name, tail)
    ccall((:EndsWithTail, libMAGEMin), Cint, (Ptr{Cchar}, Ptr{Cchar}), name, tail)
end

function RootBracketed(x1, x2)
    ccall((:RootBracketed, libMAGEMin), Cint, (Cdouble, Cdouble), x1, x2)
end

function sign(x)
    ccall((:sign, libMAGEMin), Cdouble, (Cdouble,), x)
end

function AFunction(x, data)
    ccall((:AFunction, libMAGEMin), Cdouble, (Cdouble, Ptr{Cdouble}), x, data)
end

function Minimum(x1, x2)
    ccall((:Minimum, libMAGEMin), Cdouble, (Cdouble, Cdouble), x1, x2)
end

function Maximum(x1, x2)
    ccall((:Maximum, libMAGEMin), Cdouble, (Cdouble, Cdouble), x1, x2)
end

function euclidean_distance(array1, array2, n)
    ccall((:euclidean_distance, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), array1, array2, n)
end

function partial_euclidean_distance(array1, array2, n)
    ccall((:partial_euclidean_distance, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), array1, array2, n)
end

function VecVecMul(B0, B1, n)
    ccall((:VecVecMul, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), B0, B1, n)
end

function BrentRoots(x1, x2, data, Tolerance, mode, maxIterations, valueAtRoot, niter, error)
    ccall((:BrentRoots, libMAGEMin), Cdouble, (Cdouble, Cdouble, Ptr{Cdouble}, Cdouble, Cint, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), x1, x2, data, Tolerance, mode, maxIterations, valueAtRoot, niter, error)
end

function print_cp(gv, cp)
    ccall((:print_cp, libMAGEMin), Cvoid, (Cint, Ptr{Cint}), gv, cp)
end

function print_SS_informations(gv, SS_ref_db, iss)
    ccall((:print_SS_informations, libMAGEMin), Cvoid, (Cint, Cint, Cint), gv, SS_ref_db, iss)
end

function rotate_hyperplane()
    ccall((:rotate_hyperplane, libMAGEMin), Cint, ())
end

function raw_hyperplane()
    ccall((:raw_hyperplane, libMAGEMin), Cint, ())
end

function restrict_SS_HyperVolume()
    ccall((:restrict_SS_HyperVolume, libMAGEMin), Cint, ())
end

function check_SS_bounds()
    ccall((:check_SS_bounds, libMAGEMin), Cint, ())
end

struct TMatrix
    m::Ptr{Ptr{Cdouble}}
    nRows::Cint
    nCols::Cint
end

const TMATRIX = TMatrix

function createMatrix(nRows, nCols)
    ccall((:createMatrix, libMAGEMin), TMATRIX, (Cint, Cint), nRows, nCols)
end

function rref(stoeMat, pivot, tolerance)
    ccall((:rref, libMAGEMin), TMATRIX, (TMATRIX, Ptr{Cint}, Cdouble), stoeMat, pivot, tolerance)
end

function freeMatrix(oMatrix)
    ccall((:freeMatrix, libMAGEMin), Cvoid, (TMATRIX,), oMatrix)
end

function cleanUpMatrix(stoeMat, tolerance)
    ccall((:cleanUpMatrix, libMAGEMin), Cvoid, (TMATRIX, Cdouble), stoeMat, tolerance)
end

function get_pp_id()
    ccall((:get_pp_id, libMAGEMin), Cint, ())
end

function get_ss_id()
    ccall((:get_ss_id, libMAGEMin), Cint, ())
end

struct UT_hash_handle
    tbl::Ptr{Cvoid} # tbl::Ptr{UT_hash_table}
    prev::Ptr{Cvoid}
    next::Ptr{Cvoid}
    hh_prev::Ptr{UT_hash_handle}
    hh_next::Ptr{UT_hash_handle}
    key::Ptr{Cvoid}
    keylen::Cuint
    hashv::Cuint
end

function Base.getproperty(x::UT_hash_handle, f::Symbol)
    f === :tbl && return Ptr{UT_hash_table}(getfield(x, f))
    return getfield(x, f)
end

struct UT_hash_table
    buckets::Ptr{Cvoid} # buckets::Ptr{UT_hash_bucket}
    num_buckets::Cuint
    log2_num_buckets::Cuint
    num_items::Cuint
    tail::Ptr{UT_hash_handle}
    hho::Cptrdiff_t
    ideal_chain_maxlen::Cuint
    nonideal_items::Cuint
    ineff_expands::Cuint
    noexpand::Cuint
    signature::UInt32
end

function Base.getproperty(x::UT_hash_table, f::Symbol)
    f === :buckets && return Ptr{UT_hash_bucket}(getfield(x, f))
    return getfield(x, f)
end

struct UT_hash_bucket
    hh_head::Ptr{UT_hash_handle}
    count::Cuint
    expand_mult::Cuint
end

const n_em_db = 291

const nEl = 11

const ko_no_argument = 0

const ko_required_argument = 1

const ko_optional_argument = 2

# Skipping MacroDefinition: UTHASH_VERSION 2.1.0

const HASH_NONFATAL_OOM = 0

const HASH_INITIAL_NUM_BUCKETS = Cuint(32)

const HASH_INITIAL_NUM_BUCKETS_LOG2 = Cuint(5)

const HASH_BKT_CAPACITY_THRESH = Cuint(10)

const HASH_FCN = HASH_JEN

const HASH_BLOOM_BYTELEN = Cuint(0)

const HASH_SIGNATURE = Cuint(0xa0111fe1)

const HASH_BLOOM_SIGNATURE = Cuint(0xb12220f2)

# exports
const PREFIXES = ["CX", "clang_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
