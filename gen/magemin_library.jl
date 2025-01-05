module LibMAGEMin

using MAGEMin_jll
export MAGEMin_jll

using CEnum

#
# START OF PROLOGUE
#
using MAGEMin_C, MAGEMin_jll
const HASH_JEN = 0;


function __init__()
    if isfile("libMAGEMin.dylib")
        global libMAGEMin = joinpath(pwd(),"libMAGEMin.dylib")
        println("Using locally compiled version of libMAGEMin.dylib")
    else
        global libMAGEMin = MAGEMin_jll.libMAGEMin
        println("Using libMAGEMin.dylib from MAGEMin_jll")
    end
end

#
# END OF PROLOGUE
#

mutable struct EM_db_
    Name::NTuple{20, Cchar}
    Comp::NTuple{16, Cdouble}
    input_1::NTuple{3, Cdouble}
    input_2::NTuple{4, Cdouble}
    input_3::NTuple{11, Cdouble}
    input_4::NTuple{3, Cdouble}
    EM_db_() = new()
end

const EM_db = EM_db_

mutable struct FS_db_
    Name::NTuple{20, Cchar}
    Comp::NTuple{17, Cdouble}
    input_1::NTuple{4, Cdouble}
    input_2::NTuple{7, Cdouble}
    input_3::NTuple{1, Cdouble}
    FS_db_() = new()
end

const FS_db = FS_db_

function Access_EM_DB(id, EM_dataset)
    ccall((:Access_EM_DB, libMAGEMin), EM_db, (Cint, Cint), id, EM_dataset)
end

function Access_FS_DB(id)
    ccall((:Access_FS_DB, libMAGEMin), FS_db, (Cint,), id)
end

struct PP_refs
    Name::NTuple{20, Cchar}
    Comp::NTuple{15, Cdouble}
    Comp_mol::NTuple{15, Cdouble}
    Comp_wt::NTuple{15, Cdouble}
    gbase::Cdouble
    gb_lvl::Cdouble
    factor::Cdouble
    factor_norm::Cdouble
    phase_density::Cdouble
    phase_shearModulus::Cdouble
    phase_bulkModulus::Cdouble
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_isoTbulkModulus::Cdouble
    thetaExp::Cdouble
    phase_entropy::Cdouble
    phase_enthalpy::Cdouble
    volume::Cdouble
    mass::Cdouble
    charge::Cdouble
end

const PP_ref = PP_refs

function G_EM_function(research_group, EM_dataset, len_ox, id, bulk_rock, apo, P, T, name, state)
    ccall((:G_EM_function, libMAGEMin), PP_ref, (Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), research_group, EM_dataset, len_ox, id, bulk_rock, apo, P, T, name, state)
end

function sum_array(array, size)
    ccall((:sum_array, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, size)
end

function check_sign(v1, v2)
    ccall((:check_sign, libMAGEMin), Cint, (Cdouble, Cdouble), v1, v2)
end

function TC_G_EM_function(EM_database, len_ox, id, bulk_rock, apo, P, T, name, state)
    ccall((:TC_G_EM_function, libMAGEMin), PP_ref, (Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), EM_database, len_ox, id, bulk_rock, apo, P, T, name, state)
end

mutable struct solvent_properties
    g::Cdouble
    density::Cdouble
    epsilon::Cdouble
    Z::Cdouble
    solvent_properties() = new()
end

const solvent_prop = solvent_properties

function G_FS_function(len_ox, wat, id, bulk_rock, ElH, apo, P, T, name, state)
    ccall((:G_FS_function, libMAGEMin), PP_ref, (Cint, Ptr{solvent_prop}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), len_ox, wat, id, bulk_rock, ElH, apo, P, T, name, state)
end

function propSolvent_JN91_calc(wat, TK)
    ccall((:propSolvent_JN91_calc, libMAGEMin), Cvoid, (Ptr{solvent_prop}, Cdouble), wat, TK)
end

function propSolvent_FE97_calc(wat, Pbar, TK)
    ccall((:propSolvent_FE97_calc, libMAGEMin), Cvoid, (Ptr{solvent_prop}, Cdouble, Cdouble), wat, Pbar, TK)
end

function propSolvent_SV14_calc(wat, Pbar, TK)
    ccall((:propSolvent_SV14_calc, libMAGEMin), Cvoid, (Ptr{solvent_prop}, Cdouble, Cdouble), wat, Pbar, TK)
end

function rho_wat_calc(wat, Pbar, TK, opt)
    ccall((:rho_wat_calc, libMAGEMin), Cvoid, (Ptr{solvent_prop}, Cdouble, Cdouble, Ptr{Cchar}), wat, Pbar, TK, opt)
end

mutable struct EM_db_sb_
    Name::NTuple{50, Cchar}
    FullName::NTuple{80, Cchar}
    Equation::NTuple{90, Cchar}
    Comp::NTuple{6, Cdouble}
    input_1::NTuple{10, Cdouble}
    input_2::NTuple{3, Cdouble}
    EM_db_sb_() = new()
end

const EM_db_sb = EM_db_sb_

function Access_SB_EM_DB(id, EM_dataset)
    ccall((:Access_SB_EM_DB, libMAGEMin), EM_db_sb, (Cint, Cint), id, EM_dataset)
end

function SB_G_EM_function(EM_database, len_ox, id, bulk_rock, apo, P, T, name, state)
    ccall((:SB_G_EM_function, libMAGEMin), PP_ref, (Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), EM_database, len_ox, id, bulk_rock, apo, P, T, name, state)
end

# typedef double ( * nlopt_func ) ( unsigned n , const double * x , double * gradient , /* NULL if not needed */ void * func_data )
const nlopt_func = Ptr{Cvoid}

# typedef void ( * nlopt_mfunc ) ( unsigned m , double * result , unsigned n , const double * x , double * gradient , /* NULL if not needed */ void * func_data )
const nlopt_mfunc = Ptr{Cvoid}

# typedef void ( * nlopt_precond ) ( unsigned n , const double * x , const double * v , double * vpre , void * data )
const nlopt_precond = Ptr{Cvoid}

@cenum nlopt_algorithm::UInt32 begin
    NLOPT_GN_DIRECT = 0
    NLOPT_GN_DIRECT_L = 1
    NLOPT_GN_DIRECT_L_RAND = 2
    NLOPT_GN_DIRECT_NOSCAL = 3
    NLOPT_GN_DIRECT_L_NOSCAL = 4
    NLOPT_GN_DIRECT_L_RAND_NOSCAL = 5
    NLOPT_GN_ORIG_DIRECT = 6
    NLOPT_GN_ORIG_DIRECT_L = 7
    NLOPT_GD_STOGO = 8
    NLOPT_GD_STOGO_RAND = 9
    NLOPT_LD_LBFGS_NOCEDAL = 10
    NLOPT_LD_LBFGS = 11
    NLOPT_LN_PRAXIS = 12
    NLOPT_LD_VAR1 = 13
    NLOPT_LD_VAR2 = 14
    NLOPT_LD_TNEWTON = 15
    NLOPT_LD_TNEWTON_RESTART = 16
    NLOPT_LD_TNEWTON_PRECOND = 17
    NLOPT_LD_TNEWTON_PRECOND_RESTART = 18
    NLOPT_GN_CRS2_LM = 19
    NLOPT_GN_MLSL = 20
    NLOPT_GD_MLSL = 21
    NLOPT_GN_MLSL_LDS = 22
    NLOPT_GD_MLSL_LDS = 23
    NLOPT_LD_MMA = 24
    NLOPT_LN_COBYLA = 25
    NLOPT_LN_NEWUOA = 26
    NLOPT_LN_NEWUOA_BOUND = 27
    NLOPT_LN_NELDERMEAD = 28
    NLOPT_LN_SBPLX = 29
    NLOPT_LN_AUGLAG = 30
    NLOPT_LD_AUGLAG = 31
    NLOPT_LN_AUGLAG_EQ = 32
    NLOPT_LD_AUGLAG_EQ = 33
    NLOPT_LN_BOBYQA = 34
    NLOPT_GN_ISRES = 35
    NLOPT_AUGLAG = 36
    NLOPT_AUGLAG_EQ = 37
    NLOPT_G_MLSL = 38
    NLOPT_G_MLSL_LDS = 39
    NLOPT_LD_SLSQP = 40
    NLOPT_LD_CCSAQ = 41
    NLOPT_GN_ESCH = 42
    NLOPT_GN_AGS = 43
    NLOPT_NUM_ALGORITHMS = 44
end

function nlopt_algorithm_name(a)
    ccall((:nlopt_algorithm_name, libMAGEMin), Ptr{Cchar}, (nlopt_algorithm,), a)
end

function nlopt_algorithm_to_string(algorithm)
    ccall((:nlopt_algorithm_to_string, libMAGEMin), Ptr{Cchar}, (nlopt_algorithm,), algorithm)
end

function nlopt_algorithm_from_string(name)
    ccall((:nlopt_algorithm_from_string, libMAGEMin), nlopt_algorithm, (Ptr{Cchar},), name)
end

@cenum nlopt_result::Int32 begin
    NLOPT_FAILURE = -1
    NLOPT_INVALID_ARGS = -2
    NLOPT_OUT_OF_MEMORY = -3
    NLOPT_ROUNDOFF_LIMITED = -4
    NLOPT_FORCED_STOP = -5
    NLOPT_NUM_FAILURES = -6
    NLOPT_SUCCESS = 1
    NLOPT_STOPVAL_REACHED = 2
    NLOPT_FTOL_REACHED = 3
    NLOPT_XTOL_REACHED = 4
    NLOPT_MAXEVAL_REACHED = 5
    NLOPT_MAXTIME_REACHED = 6
    NLOPT_NUM_RESULTS = 7
end

function nlopt_result_to_string(algorithm)
    ccall((:nlopt_result_to_string, libMAGEMin), Ptr{Cchar}, (nlopt_result,), algorithm)
end

function nlopt_result_from_string(name)
    ccall((:nlopt_result_from_string, libMAGEMin), nlopt_result, (Ptr{Cchar},), name)
end

function nlopt_srand(seed)
    ccall((:nlopt_srand, libMAGEMin), Cvoid, (Culong,), seed)
end

function nlopt_srand_time()
    ccall((:nlopt_srand_time, libMAGEMin), Cvoid, ())
end

function nlopt_version(major, minor, bugfix)
    ccall((:nlopt_version, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), major, minor, bugfix)
end

mutable struct nlopt_opt_s end

const nlopt_opt = Ptr{nlopt_opt_s}

function nlopt_create(algorithm, n)
    ccall((:nlopt_create, libMAGEMin), nlopt_opt, (nlopt_algorithm, Cuint), algorithm, n)
end

function nlopt_destroy(opt)
    ccall((:nlopt_destroy, libMAGEMin), Cvoid, (nlopt_opt,), opt)
end

function nlopt_copy(opt)
    ccall((:nlopt_copy, libMAGEMin), nlopt_opt, (nlopt_opt,), opt)
end

function nlopt_optimize(opt, x, opt_f)
    ccall((:nlopt_optimize, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}, Ptr{Cdouble}), opt, x, opt_f)
end

function nlopt_set_min_objective(opt, f, f_data)
    ccall((:nlopt_set_min_objective, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, Ptr{Cvoid}), opt, f, f_data)
end

function nlopt_set_max_objective(opt, f, f_data)
    ccall((:nlopt_set_max_objective, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, Ptr{Cvoid}), opt, f, f_data)
end

function nlopt_set_precond_min_objective(opt, f, pre, f_data)
    ccall((:nlopt_set_precond_min_objective, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, nlopt_precond, Ptr{Cvoid}), opt, f, pre, f_data)
end

function nlopt_set_precond_max_objective(opt, f, pre, f_data)
    ccall((:nlopt_set_precond_max_objective, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, nlopt_precond, Ptr{Cvoid}), opt, f, pre, f_data)
end

function nlopt_get_algorithm(opt)
    ccall((:nlopt_get_algorithm, libMAGEMin), nlopt_algorithm, (nlopt_opt,), opt)
end

function nlopt_get_dimension(opt)
    ccall((:nlopt_get_dimension, libMAGEMin), Cuint, (nlopt_opt,), opt)
end

function nlopt_get_errmsg(opt)
    ccall((:nlopt_get_errmsg, libMAGEMin), Ptr{Cchar}, (nlopt_opt,), opt)
end

function nlopt_set_param(opt, name, val)
    ccall((:nlopt_set_param, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cchar}, Cdouble), opt, name, val)
end

function nlopt_get_param(opt, name, defaultval)
    ccall((:nlopt_get_param, libMAGEMin), Cdouble, (nlopt_opt, Ptr{Cchar}, Cdouble), opt, name, defaultval)
end

function nlopt_has_param(opt, name)
    ccall((:nlopt_has_param, libMAGEMin), Cint, (nlopt_opt, Ptr{Cchar}), opt, name)
end

function nlopt_num_params(opt)
    ccall((:nlopt_num_params, libMAGEMin), Cuint, (nlopt_opt,), opt)
end

function nlopt_nth_param(opt, n)
    ccall((:nlopt_nth_param, libMAGEMin), Ptr{Cchar}, (nlopt_opt, Cuint), opt, n)
end

function nlopt_set_lower_bounds(opt, lb)
    ccall((:nlopt_set_lower_bounds, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, lb)
end

function nlopt_set_lower_bounds1(opt, lb)
    ccall((:nlopt_set_lower_bounds1, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, lb)
end

function nlopt_set_lower_bound(opt, i, lb)
    ccall((:nlopt_set_lower_bound, libMAGEMin), nlopt_result, (nlopt_opt, Cint, Cdouble), opt, i, lb)
end

function nlopt_get_lower_bounds(opt, lb)
    ccall((:nlopt_get_lower_bounds, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, lb)
end

function nlopt_set_upper_bounds(opt, ub)
    ccall((:nlopt_set_upper_bounds, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, ub)
end

function nlopt_set_upper_bounds1(opt, ub)
    ccall((:nlopt_set_upper_bounds1, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, ub)
end

function nlopt_set_upper_bound(opt, i, ub)
    ccall((:nlopt_set_upper_bound, libMAGEMin), nlopt_result, (nlopt_opt, Cint, Cdouble), opt, i, ub)
end

function nlopt_get_upper_bounds(opt, ub)
    ccall((:nlopt_get_upper_bounds, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, ub)
end

function nlopt_remove_inequality_constraints(opt)
    ccall((:nlopt_remove_inequality_constraints, libMAGEMin), nlopt_result, (nlopt_opt,), opt)
end

function nlopt_add_inequality_constraint(opt, fc, fc_data, tol)
    ccall((:nlopt_add_inequality_constraint, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, Ptr{Cvoid}, Cdouble), opt, fc, fc_data, tol)
end

function nlopt_add_precond_inequality_constraint(opt, fc, pre, fc_data, tol)
    ccall((:nlopt_add_precond_inequality_constraint, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, nlopt_precond, Ptr{Cvoid}, Cdouble), opt, fc, pre, fc_data, tol)
end

function nlopt_add_inequality_mconstraint(opt, m, fc, fc_data, tol)
    ccall((:nlopt_add_inequality_mconstraint, libMAGEMin), nlopt_result, (nlopt_opt, Cuint, nlopt_mfunc, Ptr{Cvoid}, Ptr{Cdouble}), opt, m, fc, fc_data, tol)
end

function nlopt_remove_equality_constraints(opt)
    ccall((:nlopt_remove_equality_constraints, libMAGEMin), nlopt_result, (nlopt_opt,), opt)
end

function nlopt_add_equality_constraint(opt, h, h_data, tol)
    ccall((:nlopt_add_equality_constraint, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, Ptr{Cvoid}, Cdouble), opt, h, h_data, tol)
end

function nlopt_add_precond_equality_constraint(opt, h, pre, h_data, tol)
    ccall((:nlopt_add_precond_equality_constraint, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_func, nlopt_precond, Ptr{Cvoid}, Cdouble), opt, h, pre, h_data, tol)
end

function nlopt_add_equality_mconstraint(opt, m, h, h_data, tol)
    ccall((:nlopt_add_equality_mconstraint, libMAGEMin), nlopt_result, (nlopt_opt, Cuint, nlopt_mfunc, Ptr{Cvoid}, Ptr{Cdouble}), opt, m, h, h_data, tol)
end

function nlopt_set_stopval(opt, stopval)
    ccall((:nlopt_set_stopval, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, stopval)
end

function nlopt_get_stopval(opt)
    ccall((:nlopt_get_stopval, libMAGEMin), Cdouble, (nlopt_opt,), opt)
end

function nlopt_set_ftol_rel(opt, tol)
    ccall((:nlopt_set_ftol_rel, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, tol)
end

function nlopt_get_ftol_rel(opt)
    ccall((:nlopt_get_ftol_rel, libMAGEMin), Cdouble, (nlopt_opt,), opt)
end

function nlopt_set_ftol_abs(opt, tol)
    ccall((:nlopt_set_ftol_abs, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, tol)
end

function nlopt_get_ftol_abs(opt)
    ccall((:nlopt_get_ftol_abs, libMAGEMin), Cdouble, (nlopt_opt,), opt)
end

function nlopt_set_xtol_rel(opt, tol)
    ccall((:nlopt_set_xtol_rel, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, tol)
end

function nlopt_get_xtol_rel(opt)
    ccall((:nlopt_get_xtol_rel, libMAGEMin), Cdouble, (nlopt_opt,), opt)
end

function nlopt_set_xtol_abs1(opt, tol)
    ccall((:nlopt_set_xtol_abs1, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, tol)
end

function nlopt_set_xtol_abs(opt, tol)
    ccall((:nlopt_set_xtol_abs, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, tol)
end

function nlopt_get_xtol_abs(opt, tol)
    ccall((:nlopt_get_xtol_abs, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, tol)
end

function nlopt_set_x_weights1(opt, w)
    ccall((:nlopt_set_x_weights1, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, w)
end

function nlopt_set_x_weights(opt, w)
    ccall((:nlopt_set_x_weights, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, w)
end

function nlopt_get_x_weights(opt, w)
    ccall((:nlopt_get_x_weights, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, w)
end

function nlopt_set_maxeval(opt, maxeval)
    ccall((:nlopt_set_maxeval, libMAGEMin), nlopt_result, (nlopt_opt, Cint), opt, maxeval)
end

function nlopt_get_maxeval(opt)
    ccall((:nlopt_get_maxeval, libMAGEMin), Cint, (nlopt_opt,), opt)
end

function nlopt_get_numevals(opt)
    ccall((:nlopt_get_numevals, libMAGEMin), Cint, (nlopt_opt,), opt)
end

function nlopt_set_maxtime(opt, maxtime)
    ccall((:nlopt_set_maxtime, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, maxtime)
end

function nlopt_get_maxtime(opt)
    ccall((:nlopt_get_maxtime, libMAGEMin), Cdouble, (nlopt_opt,), opt)
end

function nlopt_force_stop(opt)
    ccall((:nlopt_force_stop, libMAGEMin), nlopt_result, (nlopt_opt,), opt)
end

function nlopt_set_force_stop(opt, val)
    ccall((:nlopt_set_force_stop, libMAGEMin), nlopt_result, (nlopt_opt, Cint), opt, val)
end

function nlopt_get_force_stop(opt)
    ccall((:nlopt_get_force_stop, libMAGEMin), Cint, (nlopt_opt,), opt)
end

function nlopt_set_local_optimizer(opt, local_opt)
    ccall((:nlopt_set_local_optimizer, libMAGEMin), nlopt_result, (nlopt_opt, nlopt_opt), opt, local_opt)
end

function nlopt_set_population(opt, pop)
    ccall((:nlopt_set_population, libMAGEMin), nlopt_result, (nlopt_opt, Cuint), opt, pop)
end

function nlopt_get_population(opt)
    ccall((:nlopt_get_population, libMAGEMin), Cuint, (nlopt_opt,), opt)
end

function nlopt_set_vector_storage(opt, dim)
    ccall((:nlopt_set_vector_storage, libMAGEMin), nlopt_result, (nlopt_opt, Cuint), opt, dim)
end

function nlopt_get_vector_storage(opt)
    ccall((:nlopt_get_vector_storage, libMAGEMin), Cuint, (nlopt_opt,), opt)
end

function nlopt_set_default_initial_step(opt, x)
    ccall((:nlopt_set_default_initial_step, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, x)
end

function nlopt_set_initial_step(opt, dx)
    ccall((:nlopt_set_initial_step, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}), opt, dx)
end

function nlopt_set_initial_step1(opt, dx)
    ccall((:nlopt_set_initial_step1, libMAGEMin), nlopt_result, (nlopt_opt, Cdouble), opt, dx)
end

function nlopt_get_initial_step(opt, x, dx)
    ccall((:nlopt_get_initial_step, libMAGEMin), nlopt_result, (nlopt_opt, Ptr{Cdouble}, Ptr{Cdouble}), opt, x, dx)
end

# typedef void * ( * nlopt_munge ) ( void * p )
const nlopt_munge = Ptr{Cvoid}

function nlopt_set_munge(opt, munge_on_destroy, munge_on_copy)
    ccall((:nlopt_set_munge, libMAGEMin), Cvoid, (nlopt_opt, nlopt_munge, nlopt_munge), opt, munge_on_destroy, munge_on_copy)
end

# typedef void * ( * nlopt_munge2 ) ( void * p , void * data )
const nlopt_munge2 = Ptr{Cvoid}

function nlopt_munge_data(opt, munge, data)
    ccall((:nlopt_munge_data, libMAGEMin), Cvoid, (nlopt_opt, nlopt_munge2, Ptr{Cvoid}), opt, munge, data)
end

# typedef double ( * nlopt_func_old ) ( int n , const double * x , double * gradient , /* NULL if not needed */ void * func_data )
const nlopt_func_old = Ptr{Cvoid}

function nlopt_minimize(algorithm, n, f, f_data, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, maxeval, maxtime)
    ccall((:nlopt_minimize, libMAGEMin), nlopt_result, (nlopt_algorithm, Cint, nlopt_func_old, Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Cint, Cdouble), algorithm, n, f, f_data, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, maxeval, maxtime)
end

function nlopt_minimize_constrained(algorithm, n, f, f_data, m, fc, fc_data, fc_datum_size, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, maxeval, maxtime)
    ccall((:nlopt_minimize_constrained, libMAGEMin), nlopt_result, (nlopt_algorithm, Cint, nlopt_func_old, Ptr{Cvoid}, Cint, nlopt_func_old, Ptr{Cvoid}, Cptrdiff_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Cint, Cdouble), algorithm, n, f, f_data, m, fc, fc_data, fc_datum_size, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, maxeval, maxtime)
end

function nlopt_minimize_econstrained(algorithm, n, f, f_data, m, fc, fc_data, fc_datum_size, p, h, h_data, h_datum_size, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, htol_rel, htol_abs, maxeval, maxtime)
    ccall((:nlopt_minimize_econstrained, libMAGEMin), nlopt_result, (nlopt_algorithm, Cint, nlopt_func_old, Ptr{Cvoid}, Cint, nlopt_func_old, Ptr{Cvoid}, Cptrdiff_t, Cint, nlopt_func_old, Ptr{Cvoid}, Cptrdiff_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Cdouble, Cdouble, Cint, Cdouble), algorithm, n, f, f_data, m, fc, fc_data, fc_datum_size, p, h, h_data, h_datum_size, lb, ub, x, minf, minf_max, ftol_rel, ftol_abs, xtol_rel, xtol_abs, htol_rel, htol_abs, maxeval, maxtime)
end

function nlopt_get_local_search_algorithm(deriv, nonderiv, maxeval)
    ccall((:nlopt_get_local_search_algorithm, libMAGEMin), Cvoid, (Ptr{nlopt_algorithm}, Ptr{nlopt_algorithm}, Ptr{Cint}), deriv, nonderiv, maxeval)
end

function nlopt_set_local_search_algorithm(deriv, nonderiv, maxeval)
    ccall((:nlopt_set_local_search_algorithm, libMAGEMin), Cvoid, (nlopt_algorithm, nlopt_algorithm, Cint), deriv, nonderiv, maxeval)
end

function nlopt_get_stochastic_population()
    ccall((:nlopt_get_stochastic_population, libMAGEMin), Cint, ())
end

function nlopt_set_stochastic_population(pop)
    ccall((:nlopt_set_stochastic_population, libMAGEMin), Cvoid, (Cint,), pop)
end

struct ss_pc
    xeos_pc::NTuple{25, Cdouble}
end

mutable struct PC_refs
    ss_pc_xeos::Ptr{ss_pc}
    PC_refs() = new()
end

const PC_ref = PC_refs

# typedef double ( * PC_type ) ( unsigned n , const double * x , double * grad , void * SS_ref_db )
const PC_type = Ptr{Cvoid}

mutable struct oxide_datas
    n_ox::Cint
    oxName::NTuple{16, NTuple{20, Cchar}}
    elName::NTuple{16, NTuple{20, Cchar}}
    oxMass::NTuple{16, Cdouble}
    atPerOx::NTuple{16, Cdouble}
    ElEntropy::NTuple{16, Cdouble}
    OPerOx::NTuple{16, Cdouble}
    catPerOx::NTuple{16, Cdouble}
    oxide_datas() = new()
end

const oxide_data = oxide_datas

mutable struct global_variables
    version::Ptr{Cchar}
    verbose::Cint
    outpath::Ptr{Cchar}
    Mode::Cint
    pdev::Ptr{Ptr{Cdouble}}
    n_em_db::Cint
    research_group::Ptr{Cchar}
    EM_database::Cint
    EM_dataset::Cint
    n_Diff::Cint
    leveling_mode::Cint
    status::Cint
    solver::Cint
    solver_switch_T::Cdouble
    output_matlab::Cint
    tot_min_time::Cdouble
    tot_time::Cdouble
    buffer::Ptr{Cchar}
    buffer_n::Cdouble
    limitCaOpx::Cint
    CaOpxLim::Cdouble
    mbCpx::Cint
    mbIlm::Cint
    mpSp::Cint
    mpIlm::Cint
    fluidSpec::Cint
    n_fs_db::Cint
    test::Cint
    bulk_rock::Ptr{Cdouble}
    arg_bulk::Ptr{Cdouble}
    arg_gamma::Ptr{Cdouble}
    n_points::Cint
    File::Ptr{Cchar}
    db::Ptr{Cchar}
    Phase::Ptr{Cchar}
    sys_in::Ptr{Cchar}
    ipiv::Ptr{Cint}
    lwork::Cint
    work::Ptr{Cdouble}
    n_min::Ptr{Cint}
    LP::Cint
    PGE::Cint
    mean_sum_xi::Cdouble
    sigma_sum_xi::Cdouble
    min_melt_T::Cdouble
    relax_PGE_val::Cdouble
    PC_df_add::Cdouble
    PC_min_dist::Cdouble
    PC_check_val1::Cdouble
    PC_check_val2::Cdouble
    PC_checked::Cint
    check_PC1::Cint
    check_PC2::Cint
    len_pp::Cint
    len_ss::Cint
    len_ox::Cint
    maxlen_ox::Cint
    max_n_cp::Cint
    max_n_mSS::Cint
    max_ss_size_cp::Cint
    len_cp::Cint
    ox::Ptr{Ptr{Cchar}}
    gam_tot::Ptr{Cdouble}
    gam_tot_0::Ptr{Cdouble}
    delta_gam_tot::Ptr{Cdouble}
    n_flags::Cint
    PP_list::Ptr{Ptr{Cchar}}
    SS_list::Ptr{Ptr{Cchar}}
    pp_n::Ptr{Cdouble}
    pp_n_mol::Ptr{Cdouble}
    pp_n_wt::Ptr{Cdouble}
    pp_xi::Ptr{Cdouble}
    delta_pp_n::Ptr{Cdouble}
    delta_pp_xi::Ptr{Cdouble}
    pp_flags::Ptr{Ptr{Cint}}
    numPoint::Cint
    global_ite::Cint
    H2O_id::Cint
    CO2_id::Cint
    S_id::Cint
    Al2O3_id::Cint
    CaO_id::Cint
    Na2O_id::Cint
    FeO_id::Cint
    MgO_id::Cint
    SiO2_id::Cint
    K2O_id::Cint
    O_id::Cint
    TiO2_id::Cint
    Cr2O3_id::Cint
    MnO_id::Cint
    LVL_time::Cdouble
    em2ss_shift::Cdouble
    bnd_filter_pc::Cdouble
    bnd_filter_pge::Cdouble
    max_G_pc::Cdouble
    n_SS_PC::Ptr{Cint}
    SS_PC_stp::Ptr{Cdouble}
    eps_sf_pc::Cdouble
    n_pc::Cint
    n_Ppc::Cint
    max_LP_ite::Cint
    save_Ppc_val::Cdouble
    launch_PGE::Cint
    n_ss_ph::Ptr{Cint}
    verifyPC::Ptr{Cint}
    n_solvi::Ptr{Cint}
    maxgmTime::Cdouble
    obj_tol::Cdouble
    box_size_mode_PGE::Cdouble
    box_size_mode_LP::Cdouble
    maxeval::Cint
    bnd_val::Cdouble
    obj_refine_fac::Cdouble
    A_PGE::Ptr{Cdouble}
    A0_PGE::Ptr{Cdouble}
    b_PGE::Ptr{Cdouble}
    dn_cp::Ptr{Cdouble}
    dn_pp::Ptr{Cdouble}
    cp_id::Ptr{Cint}
    pp_id::Ptr{Cint}
    fc_norm_t1::Cdouble
    inner_PGE_ite::Cint
    inner_PGE_ite_time::Cdouble
    n_phase::Cint
    n_pp_phase::Cint
    n_cp_phase::Cint
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    it_1::Cint
    ur_1::Cdouble
    it_2::Cint
    ur_2::Cdouble
    it_3::Cint
    ur_3::Cdouble
    it_f::Cint
    div::Cint
    dGamma::Ptr{Cdouble}
    gibbs_ev::Ptr{Cdouble}
    PGE_mass_norm::Ptr{Cdouble}
    Alg::Ptr{Cint}
    gamma_norm::Ptr{Cdouble}
    ite_time::Ptr{Cdouble}
    G_system::Cdouble
    G_system_mu::Cdouble
    br_max_tol::Cdouble
    alpha::Cdouble
    n_ss_array::Ptr{Cint}
    ph_change::Cint
    merge_value::Cdouble
    re_in_n::Cdouble
    re_in_df::Cdouble
    min_df::Cdouble
    pc_composite_dist::Cdouble
    A::Ptr{Ptr{Cdouble}}
    A2::Ptr{Ptr{Cdouble}}
    b::Ptr{Cdouble}
    b1::Ptr{Cdouble}
    tmp1::Ptr{Cdouble}
    tmp2::Ptr{Cdouble}
    tmp3::Ptr{Cdouble}
    pc_id::Ptr{Cint}
    mass_residual::Ptr{Cdouble}
    BR_norm::Cdouble
    poisson_ratio::Cdouble
    gb_P_eps::Cdouble
    gb_T_eps::Cdouble
    system_density::Cdouble
    system_entropy::Cdouble
    system_enthalpy::Cdouble
    system_cp::Cdouble
    system_expansivity::Cdouble
    system_bulkModulus::Cdouble
    system_shearModulus::Cdouble
    system_Vp::Cdouble
    system_Vs::Cdouble
    system_volume::Cdouble
    system_fO2::Cdouble
    system_deltaQFM::Cdouble
    system_aH2O::Cdouble
    system_aSiO2::Cdouble
    system_aTiO2::Cdouble
    system_aAl2O3::Cdouble
    system_aMgO::Cdouble
    system_aFeO::Cdouble
    melt_density::Cdouble
    melt_bulkModulus::Cdouble
    melt_fraction::Cdouble
    solid_fraction::Cdouble
    solid_density::Cdouble
    solid_bulkModulus::Cdouble
    solid_shearModulus::Cdouble
    solid_Vp::Cdouble
    solid_Vs::Cdouble
    V_cor::Ptr{Cdouble}
    global_variables() = new()
end

const global_variable = global_variables

function runMAGEMin(argc, argv)
    ccall((:runMAGEMin, libMAGEMin), Cint, (Cint, Ptr{Ptr{Cchar}}), argc, argv)
end

function find_EM_id(em_tag)
    ccall((:find_EM_id, libMAGEMin), Cint, (Ptr{Cchar},), em_tag)
end

function find_FS_id(em_tag)
    ccall((:find_FS_id, libMAGEMin), Cint, (Ptr{Cchar},), em_tag)
end

# typedef double ( * obj_type ) ( unsigned n , const double * x , double * grad , void * SS_ref_db )
const obj_type = Ptr{Cvoid}

mutable struct simplex_datas
    gamma_ps::Ptr{Cdouble}
    gamma_ss::Ptr{Cdouble}
    gamma_tot::Ptr{Cdouble}
    gamma_delta::Ptr{Cdouble}
    dG_B_tol::Cdouble
    min_F_tol::Cdouble
    min_F::Cdouble
    ph2swp::Cint
    n_swp::Cint
    swp::Cint
    pivot::Ptr{Cint}
    A::Ptr{Cdouble}
    Alu::Ptr{Cdouble}
    A1::Ptr{Cdouble}
    ph_id_A::Ptr{Ptr{Cint}}
    g0_A::Ptr{Cdouble}
    dG_A::Ptr{Cdouble}
    n_vec::Ptr{Cdouble}
    stage::Ptr{Cint}
    n_Ox::Cint
    n_pp::Cint
    n_em_ss::Cint
    B::Ptr{Cdouble}
    B1::Ptr{Cdouble}
    ph_id_B::Ptr{Cint}
    g0_B::Cdouble
    dG_B::Cdouble
    simplex_datas() = new()
end

const simplex_data = simplex_datas

struct SS_refs
    P::Cdouble
    T::Cdouble
    R::Cdouble
    len_ox::Cint
    ElEntropy::Ptr{Cdouble}
    g::Cdouble
    Z::Cdouble
    densityW::Cdouble
    epsilon::Cdouble
    fName::Ptr{Cchar}
    EM_list::Ptr{Ptr{Cchar}}
    CV_list::Ptr{Ptr{Cchar}}
    SF_list::Ptr{Ptr{Cchar}}
    ss_flags::Ptr{Cint}
    n_pc::Cint
    tot_pc::Ptr{Cint}
    id_pc::Ptr{Cint}
    info::Ptr{Cint}
    G_pc::Ptr{Cdouble}
    DF_pc::Ptr{Cdouble}
    comp_pc::Ptr{Ptr{Cdouble}}
    p_pc::Ptr{Ptr{Cdouble}}
    xeos_pc::Ptr{Ptr{Cdouble}}
    factor_pc::Ptr{Cdouble}
    n_Ppc::Cint
    tot_Ppc::Cint
    id_Ppc::Cint
    info_Ppc::Ptr{Cint}
    G_Ppc::Ptr{Cdouble}
    DF_Ppc::Ptr{Cdouble}
    comp_Ppc::Ptr{Ptr{Cdouble}}
    p_Ppc::Ptr{Ptr{Cdouble}}
    mu_Ppc::Ptr{Ptr{Cdouble}}
    xeos_Ppc::Ptr{Ptr{Cdouble}}
    solvus_id::Ptr{Cint}
    n_cat::Cint
    C::Ptr{Ptr{Cdouble}}
    N::Ptr{Ptr{Cdouble}}
    Vec1::Ptr{Cdouble}
    Vec2::Ptr{Cdouble}
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
    sf_id::Cint
    Comp::Ptr{Ptr{Cdouble}}
    gbase::Ptr{Cdouble}
    mu_array::Ptr{Ptr{Cdouble}}
    gb_lvl::Ptr{Cdouble}
    factor::Cdouble
    bounds::Ptr{Ptr{Cdouble}}
    bounds_ref::Ptr{Ptr{Cdouble}}
    d_em::Ptr{Cdouble}
    z_em::Ptr{Cdouble}
    n_guess::Cint
    iguess::Ptr{Cdouble}
    dguess::Ptr{Cdouble}
    mguess::Ptr{Cdouble}
    orderVar::Cint
    idOrderVar::Ptr{Cdouble}
    status::Cint
    lb::Ptr{Cdouble}
    ub::Ptr{Cdouble}
    opt::nlopt_opt
    fbc::Cdouble
    sum_apep::Cdouble
    p::Ptr{Cdouble}
    ape::Ptr{Cdouble}
    mat_phi::Ptr{Cdouble}
    mu_Gex::Ptr{Cdouble}
    sf::Ptr{Cdouble}
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
    ElShearMod::Ptr{Cdouble}
    ElBulkMod::Ptr{Cdouble}
    density::Ptr{Cdouble}
    phase_density::Cdouble
    volume::Cdouble
    mass::Cdouble
end

const SS_ref = SS_refs

mutable struct IODATA
    n_phase::Cint
    P::Cdouble
    T::Cdouble
    in_bulk::Ptr{Cdouble}
    phase_names::Ptr{Ptr{Cchar}}
    phase_xeos::Ptr{Ptr{Cdouble}}
    phase_emp::Ptr{Ptr{Cdouble}}
    IODATA() = new()
end

const io_data = IODATA

mutable struct OUTDATA
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
    OUTDATA() = new()
end

const out_data = OUTDATA

mutable struct bulk_infos
    oxName::Ptr{Ptr{Cchar}}
    elName::Ptr{Ptr{Cchar}}
    oxMass::Ptr{Cdouble}
    atPerOx::Ptr{Cint}
    P::Cdouble
    T::Cdouble
    R::Cdouble
    bulk_rock::Ptr{Cdouble}
    bulk_rock_cat::Ptr{Cdouble}
    nzEl_val::Cint
    zEl_val::Cint
    nzEl_array::Ptr{Cint}
    zEl_array::Ptr{Cint}
    id::Ptr{Cint}
    apo::Ptr{Cdouble}
    fbc::Cdouble
    masspo::Ptr{Cdouble}
    opo::Ptr{Cdouble}
    cpo::Ptr{Cdouble}
    ElEntropy::Ptr{Cdouble}
    bulk_infos() = new()
end

const bulk_info = bulk_infos

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
    ss_n_mol::Cdouble
    ss_n_wt::Cdouble
    delta_ss_n::Cdouble
    df::Cdouble
    factor::Cdouble
    factor_norm::Cdouble
    min_time::Cdouble
    sum_xi::Cdouble
    sum_dxi::Cdouble
    p_em::Ptr{Cdouble}
    xi_em::Ptr{Cdouble}
    dguess::Ptr{Cdouble}
    xeos::Ptr{Cdouble}
    xeos_0::Ptr{Cdouble}
    xeos_1::Ptr{Cdouble}
    xeos_r::Ptr{Cdouble}
    dfx::Ptr{Cdouble}
    mu::Ptr{Cdouble}
    delta_mu::Ptr{Cdouble}
    sf::Ptr{Cdouble}
    ss_comp::Ptr{Cdouble}
    ss_comp_mol::Ptr{Cdouble}
    ss_comp_wt::Ptr{Cdouble}
    gbase::Ptr{Cdouble}
    mass::Cdouble
    volume::Cdouble
    phase_density::Cdouble
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_bulkModulus::Cdouble
    phase_isoTbulkModulus::Cdouble
    thetaExp::Cdouble
    phase_shearModulus::Cdouble
    phase_entropy::Cdouble
    phase_enthalpy::Cdouble
end

const csd_phase_set = csd_phase_sets

struct stb_SS_phases
    nOx::Cint
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    n_xeos::Cint
    n_em::Cint
    n_sf::Cint
    Comp::Ptr{Cdouble}
    compVariables::Ptr{Cdouble}
    compVariablesNames::Ptr{Ptr{Cchar}}
    siteFractions::Ptr{Cdouble}
    siteFractionsNames::Ptr{Ptr{Cchar}}
    emNames::Ptr{Ptr{Cchar}}
    emFrac::Ptr{Cdouble}
    emFrac_wt::Ptr{Cdouble}
    emChemPot::Ptr{Cdouble}
    emComp::Ptr{Ptr{Cdouble}}
    Comp_wt::Ptr{Cdouble}
    emComp_wt::Ptr{Ptr{Cdouble}}
    Comp_apfu::Ptr{Cdouble}
    emComp_apfu::Ptr{Ptr{Cdouble}}
end

const stb_SS_phase = stb_SS_phases

struct mstb_SS_phases
    ph_name::Ptr{Cchar}
    ph_type::Ptr{Cchar}
    info::Ptr{Cchar}
    ph_id::Cint
    em_id::Cint
    nOx::Cint
    n_xeos::Cint
    n_em::Cint
    G_Ppc::Cdouble
    DF_Ppc::Cdouble
    comp_Ppc::Ptr{Cdouble}
    p_Ppc::Ptr{Cdouble}
    mu_Ppc::Ptr{Cdouble}
    xeos_Ppc::Ptr{Cdouble}
end

const mstb_SS_phase = mstb_SS_phases

struct stb_PP_phases
    nOx::Cint
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Ptr{Cdouble}
    Comp_wt::Ptr{Cdouble}
    Comp_apfu::Ptr{Cdouble}
end

const stb_PP_phase = stb_PP_phases

struct stb_systems
    MAGEMin_ver::Ptr{Cchar}
    dataset::Ptr{Cchar}
    database::Ptr{Cchar}
    bulk_res_norm::Cdouble
    n_iterations::Cint
    status::Cint
    nOx::Cint
    oxides::Ptr{Ptr{Cchar}}
    elements::Ptr{Ptr{Cchar}}
    P::Cdouble
    T::Cdouble
    X::Cdouble
    bulk::Ptr{Cdouble}
    bulk_wt::Ptr{Cdouble}
    buffer::Ptr{Cchar}
    buffer_n::Cdouble
    gamma::Ptr{Cdouble}
    G::Cdouble
    M_sys::Cdouble
    rho::Cdouble
    fO2::Cdouble
    dQFM::Cdouble
    aH2O::Cdouble
    aSiO2::Cdouble
    aTiO2::Cdouble
    aAl2O3::Cdouble
    aMgO::Cdouble
    aFeO::Cdouble
    alpha::Cdouble
    cp::Cdouble
    s_cp::Cdouble
    cp_wt::Cdouble
    V::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    bulkModulus_M::Cdouble
    bulkModulus_S::Cdouble
    shearModulus_S::Cdouble
    Vp_S::Cdouble
    Vs_S::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    bulk_S::Ptr{Cdouble}
    frac_S::Cdouble
    rho_S::Cdouble
    bulk_M::Ptr{Cdouble}
    frac_M::Cdouble
    rho_M::Cdouble
    bulk_F::Ptr{Cdouble}
    frac_F::Cdouble
    rho_F::Cdouble
    bulk_S_wt::Ptr{Cdouble}
    frac_S_wt::Cdouble
    frac_S_vol::Cdouble
    bulk_M_wt::Ptr{Cdouble}
    frac_M_wt::Cdouble
    frac_M_vol::Cdouble
    bulk_F_wt::Ptr{Cdouble}
    frac_F_wt::Cdouble
    frac_F_vol::Cdouble
    n_ph::Cint
    n_PP::Cint
    n_SS::Cint
    n_mSS::Cint
    ph::Ptr{Ptr{Cchar}}
    sol_name::Ptr{Ptr{Cchar}}
    ph_frac::Ptr{Cdouble}
    ph_frac_wt::Ptr{Cdouble}
    ph_frac_1at::Ptr{Cdouble}
    ph_frac_vol::Ptr{Cdouble}
    ph_type::Ptr{Cint}
    ph_id::Ptr{Cint}
    ph_id_db::Ptr{Cint}
    SS::Ptr{stb_SS_phase}
    mSS::Ptr{mstb_SS_phase}
    PP::Ptr{stb_PP_phase}
end

const stb_system = stb_systems

function global_variable_alloc(z_b)
    ccall((:global_variable_alloc, libMAGEMin), global_variable, (Ptr{bulk_info},), z_b)
end

function global_variable_init(gv, z_b)
    ccall((:global_variable_init, libMAGEMin), global_variable, (global_variable, Ptr{bulk_info}), gv, z_b)
end

function get_bulk_igneous(gv)
    ccall((:get_bulk_igneous, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_metapelite(gv)
    ccall((:get_bulk_metapelite, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_ultramafic(gv)
    ccall((:get_bulk_ultramafic, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_ultramafic_ext(gv)
    ccall((:get_bulk_ultramafic_ext, libMAGEMin), global_variable, (global_variable,), gv)
end

mutable struct Database
    PP_ref_db::Ptr{PP_ref}
    SS_ref_db::Ptr{SS_ref}
    cp::Ptr{csd_phase_set}
    sp::Ptr{stb_system}
    EM_names::Ptr{Ptr{Cchar}}
    FS_names::Ptr{Ptr{Cchar}}
    Database() = new()
end

const Databases = Database

function InitializeDatabases(gv, EM_database)
    ccall((:InitializeDatabases, libMAGEMin), Databases, (global_variable, Cint), gv, EM_database)
end

function FreeDatabases(gv, DB, z_b)
    ccall((:FreeDatabases, libMAGEMin), Cvoid, (global_variable, Databases, bulk_info), gv, DB, z_b)
end

function ComputeG0_point(EM_database, z_b, gv, PP_ref_db, SS_ref_db)
    ccall((:ComputeG0_point, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), EM_database, z_b, gv, PP_ref_db, SS_ref_db)
end

function ComputeEquilibrium_Point(EM_database, input_data, z_b, gv, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:ComputeEquilibrium_Point, libMAGEMin), global_variable, (Cint, io_data, bulk_info, global_variable, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, input_data, z_b, gv, splx_data, PP_ref_db, SS_ref_db, cp)
end

function ComputeLevellingOnly(EM_database, input_data, z_b, gv, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:ComputeLevellingOnly, libMAGEMin), global_variable, (Cint, io_data, bulk_info, global_variable, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, input_data, z_b, gv, splx_data, PP_ref_db, SS_ref_db, cp)
end

function ComputePostProcessing(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:ComputePostProcessing, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function ReadCommandLineOptions(gv, z_b, argc, argv)
    ccall((:ReadCommandLineOptions, libMAGEMin), global_variable, (global_variable, Ptr{bulk_info}, Cint, Ptr{Ptr{Cchar}}), gv, z_b, argc, argv)
end

function SetupDatabase(gv, z_b)
    ccall((:SetupDatabase, libMAGEMin), global_variable, (global_variable, Ptr{bulk_info}), gv, z_b)
end

function PrintOutput(gv, rank, l, DB, time_taken, z_b)
    ccall((:PrintOutput, libMAGEMin), Cvoid, (global_variable, Cint, Cint, Databases, Cdouble, bulk_info), gv, rank, l, DB, time_taken, z_b)
end

function PrintStatus(status)
    ccall((:PrintStatus, libMAGEMin), Cvoid, (Cint,), status)
end

struct UT_hash_bucket
    hh_head::Ptr{Cvoid} # hh_head::Ptr{UT_hash_handle}
    count::Cuint
    expand_mult::Cuint
end

function Base.getproperty(x::UT_hash_bucket, f::Symbol)
    f === :hh_head && return Ptr{UT_hash_handle}(getfield(x, f))
    return getfield(x, f)
end

struct UT_hash_table
    buckets::Ptr{UT_hash_bucket}
    num_buckets::Cuint
    log2_num_buckets::Cuint
    num_items::Cuint
    tail::Ptr{Cvoid} # tail::Ptr{UT_hash_handle}
    hho::Cptrdiff_t
    ideal_chain_maxlen::Cuint
    nonideal_items::Cuint
    ineff_expands::Cuint
    noexpand::Cuint
    signature::UInt32
end

function Base.getproperty(x::UT_hash_table, f::Symbol)
    f === :tail && return Ptr{UT_hash_handle}(getfield(x, f))
    return getfield(x, f)
end

struct UT_hash_handle
    tbl::Ptr{UT_hash_table}
    prev::Ptr{Cvoid}
    next::Ptr{Cvoid}
    hh_prev::Ptr{UT_hash_handle}
    hh_next::Ptr{UT_hash_handle}
    key::Ptr{Cvoid}
    keylen::Cuint
    hashv::Cuint
end

mutable struct metapelite_aq_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{12, NTuple{20, Cchar}}
    PP::NTuple{24, NTuple{20, Cchar}}
    SS::NTuple{18, NTuple{20, Cchar}}
    verifyPC::NTuple{18, Cint}
    n_SS_PC::NTuple{18, Cint}
    SS_PC_stp::NTuple{18, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    metapelite_aq_datasets() = new()
end

const metapelite_aq_dataset = metapelite_aq_datasets

mutable struct metapelite_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{11, NTuple{20, Cchar}}
    PP::NTuple{24, NTuple{20, Cchar}}
    SS::NTuple{18, NTuple{20, Cchar}}
    verifyPC::NTuple{18, Cint}
    n_SS_PC::NTuple{18, Cint}
    SS_PC_stp::NTuple{18, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    metapelite_datasets() = new()
end

const metapelite_dataset = metapelite_datasets

mutable struct metabasite_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{10, NTuple{20, Cchar}}
    PP::NTuple{24, NTuple{20, Cchar}}
    SS::NTuple{17, NTuple{20, Cchar}}
    verifyPC::NTuple{17, Cint}
    n_SS_PC::NTuple{17, Cint}
    SS_PC_stp::NTuple{17, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    metabasite_datasets() = new()
end

const metabasite_dataset = metabasite_datasets

mutable struct igneous_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{11, NTuple{20, Cchar}}
    PP::NTuple{24, NTuple{20, Cchar}}
    SS::NTuple{15, NTuple{20, Cchar}}
    verifyPC::NTuple{15, Cint}
    n_SS_PC::NTuple{15, Cint}
    SS_PC_stp::NTuple{15, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    igneous_datasets() = new()
end

const igneous_dataset = igneous_datasets

mutable struct igneous_igad_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{10, NTuple{20, Cchar}}
    PP::NTuple{23, NTuple{20, Cchar}}
    SS::NTuple{12, NTuple{20, Cchar}}
    verifyPC::NTuple{12, Cint}
    n_SS_PC::NTuple{12, Cint}
    SS_PC_stp::NTuple{12, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    igneous_igad_datasets() = new()
end

const igneous_igad_dataset = igneous_igad_datasets

mutable struct ultramafic_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{7, NTuple{20, Cchar}}
    PP::NTuple{21, NTuple{20, Cchar}}
    SS::NTuple{12, NTuple{20, Cchar}}
    verifyPC::NTuple{12, Cint}
    n_SS_PC::NTuple{12, Cint}
    SS_PC_stp::NTuple{12, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    ultramafic_datasets() = new()
end

const ultramafic_dataset = ultramafic_datasets

mutable struct ultramafic_ext_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{9, NTuple{20, Cchar}}
    PP::NTuple{21, NTuple{20, Cchar}}
    SS::NTuple{15, NTuple{20, Cchar}}
    verifyPC::NTuple{15, Cint}
    n_SS_PC::NTuple{15, Cint}
    SS_PC_stp::NTuple{15, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    ultramafic_ext_datasets() = new()
end

const ultramafic_ext_dataset = ultramafic_ext_datasets

mutable struct mantle_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{6, NTuple{20, Cchar}}
    PP::NTuple{8, NTuple{20, Cchar}}
    SS::NTuple{14, NTuple{20, Cchar}}
    verifyPC::NTuple{14, Cint}
    n_SS_PC::NTuple{14, Cint}
    SS_PC_stp::NTuple{14, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    mantle_datasets() = new()
end

const mantle_dataset = mantle_datasets

mutable struct metapelite_datasets_ext
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{13, NTuple{20, Cchar}}
    PP::NTuple{26, NTuple{20, Cchar}}
    SS::NTuple{23, NTuple{20, Cchar}}
    verifyPC::NTuple{23, Cint}
    n_SS_PC::NTuple{23, Cint}
    SS_PC_stp::NTuple{23, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    metapelite_datasets_ext() = new()
end

const metapelite_dataset_ext = metapelite_datasets_ext

function global_variable_TC_init(gv, z_b)
    ccall((:global_variable_TC_init, libMAGEMin), global_variable, (global_variable, Ptr{bulk_info}), gv, z_b)
end

function get_bulk_metabasite(gv)
    ccall((:get_bulk_metabasite, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_igneous_igad(gv)
    ccall((:get_bulk_igneous_igad, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_mantle(gv)
    ccall((:get_bulk_mantle, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_bulk_metapelite_ext(gv)
    ccall((:get_bulk_metapelite_ext, libMAGEMin), global_variable, (global_variable,), gv)
end

mutable struct stx11_datasets
    ds_version::Cint
    n_ox::Cint
    n_pp::Cint
    n_ss::Cint
    ox::NTuple{6, NTuple{20, Cchar}}
    PP::NTuple{10, NTuple{20, Cchar}}
    SS::NTuple{14, NTuple{20, Cchar}}
    verifyPC::NTuple{14, Cint}
    n_SS_PC::NTuple{14, Cint}
    SS_PC_stp::NTuple{14, Cdouble}
    PC_df_add::Cdouble
    solver_switch_T::Cdouble
    min_melt_T::Cdouble
    inner_PGE_ite::Cdouble
    max_n_phase::Cdouble
    max_g_phase::Cdouble
    max_fac::Cdouble
    merge_value::Cdouble
    re_in_n::Cdouble
    obj_tol::Cdouble
    stx11_datasets() = new()
end

const stx11_dataset = stx11_datasets

function global_variable_SB_init(gv, z_b)
    ccall((:global_variable_SB_init, libMAGEMin), global_variable, (global_variable, Ptr{bulk_info}), gv, z_b)
end

function get_bulk_stx11(gv)
    ccall((:get_bulk_stx11, libMAGEMin), global_variable, (global_variable,), gv)
end

# typedef SS_ref ( * SS_init_type ) ( SS_ref SS_ref_db , global_variable gv )
const SS_init_type = Ptr{Cvoid}

function get_EM_DB_names_tc(gv)
    ccall((:get_EM_DB_names_tc, libMAGEMin), Ptr{Ptr{Cchar}}, (global_variable,), gv)
end

function get_EM_DB_names_sb(gv)
    ccall((:get_EM_DB_names_sb, libMAGEMin), Ptr{Ptr{Cchar}}, (global_variable,), gv)
end

function get_FS_DB_names(gv)
    ccall((:get_FS_DB_names, libMAGEMin), Ptr{Ptr{Cchar}}, (global_variable,), gv)
end

mutable struct em_datas
    C::NTuple{14, Cdouble}
    ElShearMod::Cdouble
    ElBulkMod::Cdouble
    gb::Cdouble
    charge::Cdouble
    em_datas() = new()
end

const em_data = em_datas

function get_em_data(research_group, EM_dataset, len_ox, z_b, P, T, name, state)
    ccall((:get_em_data, libMAGEMin), em_data, (Ptr{Cchar}, Cint, Cint, bulk_info, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), research_group, EM_dataset, len_ox, z_b, P, T, name, state)
end

function get_fs_data(len_ox, z_b, wat, P, T, name, state)
    ccall((:get_fs_data, libMAGEMin), em_data, (Cint, bulk_info, Ptr{solvent_prop}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), len_ox, z_b, wat, P, T, name, state)
end

function CP_INIT_function(cp, gv)
    ccall((:CP_INIT_function, libMAGEMin), csd_phase_set, (csd_phase_set, global_variable), cp, gv)
end

function SP_INIT_function(sp, gv)
    ccall((:SP_INIT_function, libMAGEMin), stb_system, (stb_system, global_variable), sp, gv)
end

function G_SS_init_EM_function(SS_init, ph_id, SS_ref_db, name, gv)
    ccall((:G_SS_init_EM_function, libMAGEMin), SS_ref, (Ptr{SS_init_type}, Cint, SS_ref, Ptr{Cchar}, global_variable), SS_init, ph_id, SS_ref_db, name, gv)
end

function reset_gv(gv, z_b, PP_ref_db, SS_ref_db)
    ccall((:reset_gv, libMAGEMin), global_variable, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}), gv, z_b, PP_ref_db, SS_ref_db)
end

function reset_sp(gv, sp)
    ccall((:reset_sp, libMAGEMin), Cvoid, (global_variable, Ptr{stb_system}), gv, sp)
end

function reset_z_b_bulk(gv, z_b)
    ccall((:reset_z_b_bulk, libMAGEMin), bulk_info, (global_variable, bulk_info), gv, z_b)
end

function reset_cp(gv, z_b, cp)
    ccall((:reset_cp, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{csd_phase_set}), gv, z_b, cp)
end

function reset_SS(gv, z_b, SS_ref_db)
    ccall((:reset_SS, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{SS_ref}), gv, z_b, SS_ref_db)
end

function init_simplex_A(splx_data, gv)
    ccall((:init_simplex_A, libMAGEMin), Cvoid, (Ptr{simplex_data}, global_variable), splx_data, gv)
end

function init_simplex_B_em(splx_data, gv)
    ccall((:init_simplex_B_em, libMAGEMin), Cvoid, (Ptr{simplex_data}, global_variable), splx_data, gv)
end

function reset_simplex_A(splx_data, z_b, gv)
    ccall((:reset_simplex_A, libMAGEMin), Cvoid, (Ptr{simplex_data}, bulk_info, global_variable), splx_data, z_b, gv)
end

function reset_simplex_B_em(splx_data, gv)
    ccall((:reset_simplex_B_em, libMAGEMin), Cvoid, (Ptr{simplex_data}, global_variable), splx_data, gv)
end

function TC_SS_init_mp(SS_init, gv)
    ccall((:TC_SS_init_mp, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_mb(SS_init, gv)
    ccall((:TC_SS_init_mb, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_ig(SS_init, gv)
    ccall((:TC_SS_init_ig, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_igad(SS_init, gv)
    ccall((:TC_SS_init_igad, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_um(SS_init, gv)
    ccall((:TC_SS_init_um, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_um_ext(SS_init, gv)
    ccall((:TC_SS_init_um_ext, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_mp_ext(SS_init, gv)
    ccall((:TC_SS_init_mp_ext, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init_mtl(SS_init, gv)
    ccall((:TC_SS_init_mtl, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function TC_SS_init(SS_init, gv)
    ccall((:TC_SS_init, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function G_SS_mp_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_mp_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_mb_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_mb_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_ig_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_ig_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_igad_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_igad_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_um_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_um_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_um_ext_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_um_ext_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_mtl_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_mtl_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function G_SS_mpe_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_mpe_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function TC_mb_objective_init_function(SS_objective, gv)
    ccall((:TC_mb_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_ig_objective_init_function(SS_objective, gv)
    ccall((:TC_ig_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_igad_objective_init_function(SS_objective, gv)
    ccall((:TC_igad_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_mp_objective_init_function(SS_objective, gv)
    ccall((:TC_mp_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_um_objective_init_function(SS_objective, gv)
    ccall((:TC_um_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_um_ext_objective_init_function(SS_objective, gv)
    ccall((:TC_um_ext_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_mtl_objective_init_function(SS_objective, gv)
    ccall((:TC_mtl_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_mpe_objective_init_function(SS_objective, gv)
    ccall((:TC_mpe_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_SS_objective_init_function(SS_objective, gv)
    ccall((:TC_SS_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function TC_PC_init(PC_read, gv)
    ccall((:TC_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

# typedef void ( * P2X_type ) ( void * SS_ref_db , double eps )
const P2X_type = Ptr{Cvoid}

function TC_P2X_init(P2X_read, gv)
    ccall((:TC_P2X_init, libMAGEMin), Cvoid, (Ptr{P2X_type}, global_variable), P2X_read, gv)
end

function obj_mb_liq(n, x, grad, SS_ref_db)
    ccall((:obj_mb_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_amp(n, x, grad, SS_ref_db)
    ccall((:obj_mb_amp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_aug(n, x, grad, SS_ref_db)
    ccall((:obj_mb_aug, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_dio(n, x, grad, SS_ref_db)
    ccall((:obj_mb_dio, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_opx(n, x, grad, SS_ref_db)
    ccall((:obj_mb_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_g(n, x, grad, SS_ref_db)
    ccall((:obj_mb_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_ol(n, x, grad, SS_ref_db)
    ccall((:obj_mb_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_fsp(n, x, grad, SS_ref_db)
    ccall((:obj_mb_fsp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_abc(n, x, grad, SS_ref_db)
    ccall((:obj_mb_abc, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_k4tr(n, x, grad, SS_ref_db)
    ccall((:obj_mb_k4tr, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_sp(n, x, grad, SS_ref_db)
    ccall((:obj_mb_sp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_spl(n, x, grad, SS_ref_db)
    ccall((:obj_mb_spl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_mb_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_ilmm(n, x, grad, SS_ref_db)
    ccall((:obj_mb_ilmm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_ep(n, x, grad, SS_ref_db)
    ccall((:obj_mb_ep, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_bi(n, x, grad, SS_ref_db)
    ccall((:obj_mb_bi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_mu(n, x, grad, SS_ref_db)
    ccall((:obj_mb_mu, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mb_chl(n, x, grad, SS_ref_db)
    ccall((:obj_mb_chl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_fper(n, x, grad, SS_ref_db)
    ccall((:obj_ig_fper, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_bi(n, x, grad, SS_ref_db)
    ccall((:obj_ig_bi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_cd(n, x, grad, SS_ref_db)
    ccall((:obj_ig_cd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_cpx(n, x, grad, SS_ref_db)
    ccall((:obj_ig_cpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_ep(n, x, grad, SS_ref_db)
    ccall((:obj_ig_ep, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_fl(n, x, grad, SS_ref_db)
    ccall((:obj_ig_fl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_g(n, x, grad, SS_ref_db)
    ccall((:obj_ig_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_amp(n, x, grad, SS_ref_db)
    ccall((:obj_ig_amp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_ig_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_liq(n, x, grad, SS_ref_db)
    ccall((:obj_ig_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_mu(n, x, grad, SS_ref_db)
    ccall((:obj_ig_mu, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_ol(n, x, grad, SS_ref_db)
    ccall((:obj_ig_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_opx(n, x, grad, SS_ref_db)
    ccall((:obj_ig_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_fsp(n, x, grad, SS_ref_db)
    ccall((:obj_ig_fsp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ig_spl(n, x, grad, SS_ref_db)
    ccall((:obj_ig_spl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_liq(n, x, grad, SS_ref_db)
    ccall((:obj_igad_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_fsp(n, x, grad, SS_ref_db)
    ccall((:obj_igad_fsp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_spl(n, x, grad, SS_ref_db)
    ccall((:obj_igad_spl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_g(n, x, grad, SS_ref_db)
    ccall((:obj_igad_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_ol(n, x, grad, SS_ref_db)
    ccall((:obj_igad_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_opx(n, x, grad, SS_ref_db)
    ccall((:obj_igad_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_cpx(n, x, grad, SS_ref_db)
    ccall((:obj_igad_cpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_igad_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_nph(n, x, grad, SS_ref_db)
    ccall((:obj_igad_nph, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_lct(n, x, grad, SS_ref_db)
    ccall((:obj_igad_lct, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_kals(n, x, grad, SS_ref_db)
    ccall((:obj_igad_kals, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_igad_mel(n, x, grad, SS_ref_db)
    ccall((:obj_igad_mel, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_liq(n, x, grad, SS_ref_db)
    ccall((:obj_mp_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_fsp(n, x, grad, SS_ref_db)
    ccall((:obj_mp_fsp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_bi(n, x, grad, SS_ref_db)
    ccall((:obj_mp_bi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_g(n, x, grad, SS_ref_db)
    ccall((:obj_mp_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_ep(n, x, grad, SS_ref_db)
    ccall((:obj_mp_ep, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_ma(n, x, grad, SS_ref_db)
    ccall((:obj_mp_ma, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_mu(n, x, grad, SS_ref_db)
    ccall((:obj_mp_mu, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_opx(n, x, grad, SS_ref_db)
    ccall((:obj_mp_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_sa(n, x, grad, SS_ref_db)
    ccall((:obj_mp_sa, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_cd(n, x, grad, SS_ref_db)
    ccall((:obj_mp_cd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_st(n, x, grad, SS_ref_db)
    ccall((:obj_mp_st, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_chl(n, x, grad, SS_ref_db)
    ccall((:obj_mp_chl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_ctd(n, x, grad, SS_ref_db)
    ccall((:obj_mp_ctd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_sp(n, x, grad, SS_ref_db)
    ccall((:obj_mp_sp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_mp_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_ilmm(n, x, grad, SS_ref_db)
    ccall((:obj_mp_ilmm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mp_mt(n, x, grad, SS_ref_db)
    ccall((:obj_mp_mt, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_fluid(n, x, grad, SS_ref_db)
    ccall((:obj_um_fluid, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_ol(n, x, grad, SS_ref_db)
    ccall((:obj_um_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_br(n, x, grad, SS_ref_db)
    ccall((:obj_um_br, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_ch(n, x, grad, SS_ref_db)
    ccall((:obj_um_ch, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_atg(n, x, grad, SS_ref_db)
    ccall((:obj_um_atg, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_g(n, x, grad, SS_ref_db)
    ccall((:obj_um_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_ta(n, x, grad, SS_ref_db)
    ccall((:obj_um_ta, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_chl(n, x, grad, SS_ref_db)
    ccall((:obj_um_chl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_anth(n, x, grad, SS_ref_db)
    ccall((:obj_um_anth, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_spi(n, x, grad, SS_ref_db)
    ccall((:obj_um_spi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_opx(n, x, grad, SS_ref_db)
    ccall((:obj_um_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_um_po(n, x, grad, SS_ref_db)
    ccall((:obj_um_po, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ume_pl4tr(n, x, grad, SS_ref_db)
    ccall((:obj_ume_pl4tr, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ume_amp(n, x, grad, SS_ref_db)
    ccall((:obj_ume_amp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_ume_aug(n, x, grad, SS_ref_db)
    ccall((:obj_ume_aug, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_g(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_fp(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_fp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_mpv(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_mpv, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_cpv(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_cpv, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_crn(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_crn, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_cf(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_cf, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_nal(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_nal, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_aki(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_aki, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_ol(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_wad(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_wad, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_ring(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_ring, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_cpx(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_cpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_opx(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mtl_hpx(n, x, grad, SS_ref_db)
    ccall((:obj_mtl_hpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_liq(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_liq, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_fsp(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_fsp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_bi(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_bi, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_g(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_g, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_ep(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_ep, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_ma(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_ma, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_mu(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_mu, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_opx(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_sa(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_sa, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_cd(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_cd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_st(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_st, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_chl(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_chl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_ctd(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_ctd, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_sp(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_sp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_ilm(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_ilm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_ilmm(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_ilmm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_mt(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_mt, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_fl(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_fl, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_occm(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_occm, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_dio(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_dio, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_aug(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_aug, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_amp(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_amp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_mpe_po(n, x, grad, SS_ref_db)
    ccall((:obj_mpe_po, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_aq17(n, x, grad, SS_ref_db)
    ccall((:obj_aq17, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function PC_function(gv, PC_read, SS_ref_db, z_b, ph_id)
    ccall((:PC_function, libMAGEMin), SS_ref, (global_variable, Ptr{PC_type}, SS_ref, bulk_info, Cint), gv, PC_read, SS_ref_db, z_b, ph_id)
end

function TC_mp_PC_init(PC_read, gv)
    ccall((:TC_mp_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_mb_PC_init(PC_read, gv)
    ccall((:TC_mb_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_ig_PC_init(PC_read, gv)
    ccall((:TC_ig_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_igad_PC_init(PC_read, gv)
    ccall((:TC_igad_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_um_PC_init(PC_read, gv)
    ccall((:TC_um_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_um_ext_PC_init(PC_read, gv)
    ccall((:TC_um_ext_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_mtl_PC_init(PC_read, gv)
    ccall((:TC_mtl_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function TC_mpe_PC_init(PC_read, gv)
    ccall((:TC_mpe_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

# typedef SS_ref ( * NLopt_type ) ( global_variable gv , SS_ref SS_ref_db )
const NLopt_type = Ptr{Cvoid}

function TC_mp_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_mp_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_mb_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_mb_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_ig_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_ig_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_igad_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_igad_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_um_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_um_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_um_ext_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_um_ext_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_mtl_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_mtl_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function TC_mpe_NLopt_opt_init(NLopt_opt, gv)
    ccall((:TC_mpe_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function NLopt_opt_ig_spl_function(gv, SS_ref_db)
    ccall((:NLopt_opt_ig_spl_function, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function NLopt_opt_ig_cpx_function(gv, SS_ref_db)
    ccall((:NLopt_opt_ig_cpx_function, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function NLopt_opt_ig_amp_function(gv, SS_ref_db)
    ccall((:NLopt_opt_ig_amp_function, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function SS_mp_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_mp_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_mb_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_mb_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_ig_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_ig_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_igad_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_igad_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_um_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_um_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_mtl_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_mtl_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SS_mpe_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SS_mpe_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function SB_SS_init(SS_init, gv)
    ccall((:SB_SS_init, libMAGEMin), Cvoid, (Ptr{SS_init_type}, global_variable), SS_init, gv)
end

function G_SS_sb11_EM_function(gv, SS_ref_db, EM_dataset, z_b, name)
    ccall((:G_SS_sb11_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_dataset, z_b, name)
end

function SB_sb11_objective_init_function(SS_objective, gv)
    ccall((:SB_sb11_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function SB_SS_objective_init_function(SS_objective, gv)
    ccall((:SB_SS_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function SB_PC_init(PC_read, gv)
    ccall((:SB_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function obj_sb11_plg(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_plg, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_sp(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_sp, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_ol(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_ol, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_wa(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_wa, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_ri(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_ri, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_opx(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_opx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_cpx(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_cpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_hpcpx(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_hpcpx, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_ak(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_ak, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_gtmj(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_gtmj, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_pv(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_pv, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_ppv(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_ppv, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_mw(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_mw, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function obj_sb11_cf(n, x, grad, SS_ref_db)
    ccall((:obj_sb11_cf, libMAGEMin), Cdouble, (Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}), n, x, grad, SS_ref_db)
end

function SB_PC__function(gv, PC_read, SS_ref_db, z_b, ph_id)
    ccall((:SB_PC__function, libMAGEMin), SS_ref, (global_variable, Ptr{PC_type}, SS_ref, bulk_info, Cint), gv, PC_read, SS_ref_db, z_b, ph_id)
end

function SB_sb11_PC_init(PC_read, gv)
    ccall((:SB_sb11_PC_init, libMAGEMin), Cvoid, (Ptr{PC_type}, global_variable), PC_read, gv)
end

function SB_NLopt_opt_init(NLopt_opt, gv)
    ccall((:SB_NLopt_opt_init, libMAGEMin), Cvoid, (Ptr{NLopt_type}, global_variable), NLopt_opt, gv)
end

function SB_sb11_pc_init_function(SS_pc_xeos, iss, name)
    ccall((:SB_sb11_pc_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_pc_xeos, iss, name)
end

function PGE(z_b, gv, PC_read, SS_objective, NLopt_opt, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:PGE, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PC_type}, Ptr{obj_type}, Ptr{NLopt_type}, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PC_read, SS_objective, NLopt_opt, splx_data, PP_ref_db, SS_ref_db, cp)
end

function init_LP(z_b, splx_data, gv, PC_read, P2X_read, PP_ref_db, SS_ref_db, cp)
    ccall((:init_LP, libMAGEMin), global_variable, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PC_type}, Ptr{P2X_type}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, splx_data, gv, PC_read, P2X_read, PP_ref_db, SS_ref_db, cp)
end

function run_LP(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:run_LP, libMAGEMin), global_variable, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function run_LP_ig(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:run_LP_ig, libMAGEMin), global_variable, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function LP(z_b, gv, PC_read, P2X_read, SS_objective, NLopt_opt, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:LP, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PC_type}, Ptr{P2X_type}, Ptr{obj_type}, Ptr{NLopt_type}, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PC_read, P2X_read, SS_objective, NLopt_opt, splx_data, PP_ref_db, SS_ref_db, cp)
end

function norm_vector(array, n)
    ccall((:norm_vector, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, n)
end

function dump_init(gv)
    ccall((:dump_init, libMAGEMin), Cvoid, (global_variable,), gv)
end

function fill_output_struct(gv, splx_data, z_b, PP_ref_db, SS_ref_db, cp, sp)
    ccall((:fill_output_struct, libMAGEMin), Cvoid, (global_variable, Ptr{simplex_data}, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}, Ptr{stb_system}), gv, splx_data, z_b, PP_ref_db, SS_ref_db, cp, sp)
end

function save_results_function(gv, z_b, PP_ref_db, SS_ref_db, cp, sp)
    ccall((:save_results_function, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}, Ptr{stb_system}), gv, z_b, PP_ref_db, SS_ref_db, cp, sp)
end

function mergeParallelFiles(gv)
    ccall((:mergeParallelFiles, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_matlab(gv)
    ccall((:mergeParallel_matlab, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_residual_Files(gv)
    ccall((:mergeParallel_residual_Files, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_LocalMinima_Files(gv)
    ccall((:mergeParallel_LocalMinima_Files, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_LevellingGamma_Files(gv)
    ccall((:mergeParallel_LevellingGamma_Files, libMAGEMin), Cvoid, (global_variable,), gv)
end

mutable struct EM2id_
    EM_tag::NTuple{20, Cchar}
    id::Cint
    hh::UT_hash_handle
    EM2id_() = new()
end

const EM2id = EM2id_

mutable struct FS2id_
    FS_tag::NTuple{20, Cchar}
    id::Cint
    hh::UT_hash_handle
    FS2id_() = new()
end

const FS2id = FS2id_

mutable struct PP2id_
    PP_tag::NTuple{20, Cchar}
    id::Cint
    hh::UT_hash_handle
    PP2id_() = new()
end

const PP2id = PP2id_

function find_PP_id(PP_tag)
    ccall((:find_PP_id, libMAGEMin), Cint, (Ptr{Cchar},), PP_tag)
end

function read_in_data(gv, input_data, n_points)
    ccall((:read_in_data, libMAGEMin), Cvoid, (global_variable, Ptr{io_data}, Cint), gv, input_data, n_points)
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

function check_PC(z_b, gv, PC_read, PP_ref_db, SS_ref_db, cp)
    ccall((:check_PC, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PC_type}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PC_read, PP_ref_db, SS_ref_db, cp)
end

function check_PC_driving_force(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:check_PC_driving_force, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function phase_update_function(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:phase_update_function, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function phase_merge_function(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:phase_merge_function, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

mutable struct str
    value::Cdouble
    index::Cint
    str() = new()
end

function cmp_dbl(a, b)
    ccall((:cmp_dbl, libMAGEMin), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), a, b)
end

function cmp_int(a, b)
    ccall((:cmp_int, libMAGEMin), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), a, b)
end

function pp_min_function(gv, z_b, PP_ref_db)
    ccall((:pp_min_function, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{PP_ref}), gv, z_b, PP_ref_db)
end

function init_em_db(EM_database, z_b, gv, PP_ref_db)
    ccall((:init_em_db, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{PP_ref}), EM_database, z_b, gv, PP_ref_db)
end

function init_em_db_sb(EM_database, z_b, gv, PP_ref_db)
    ccall((:init_em_db_sb, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{PP_ref}), EM_database, z_b, gv, PP_ref_db)
end

function update_dG(splx_data)
    ccall((:update_dG, libMAGEMin), Cvoid, (Ptr{simplex_data},), splx_data)
end

function update_local_gamma(A1, g0_A, gam, n)
    ccall((:update_local_gamma, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), A1, g0_A, gam, n)
end

function update_global_gamma(z_b, splx_data)
    ccall((:update_global_gamma, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}), z_b, splx_data)
end

function update_global_gamma_LU(z_b, splx_data)
    ccall((:update_global_gamma_LU, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}), z_b, splx_data)
end

function swap_pure_phases(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:swap_pure_phases, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function swap_pure_endmembers(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:swap_pure_endmembers, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function swap_pseudocompounds(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:swap_pseudocompounds, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function swap_PGE_pseudocompounds(z_b, splx_data, gv, PP_ref_db, SS_ref_db)
    ccall((:swap_PGE_pseudocompounds, libMAGEMin), Cvoid, (bulk_info, Ptr{simplex_data}, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, splx_data, gv, PP_ref_db, SS_ref_db)
end

function Levelling(z_b, gv, PC_read, P2X_read, SS_objective, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:Levelling, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PC_type}, Ptr{P2X_type}, Ptr{obj_type}, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PC_read, P2X_read, SS_objective, splx_data, PP_ref_db, SS_ref_db, cp)
end

function Initial_guess(z_b, gv, PC_read, P2X_read, splx_data, PP_ref_db, SS_ref_db, cp)
    ccall((:Initial_guess, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PC_type}, Ptr{P2X_type}, Ptr{simplex_data}, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PC_read, P2X_read, splx_data, PP_ref_db, SS_ref_db, cp)
end

function destroy_simplex_A(splx_data)
    ccall((:destroy_simplex_A, libMAGEMin), Cvoid, (Ptr{simplex_data},), splx_data)
end

function destroy_simplex_B(splx_data)
    ccall((:destroy_simplex_B, libMAGEMin), Cvoid, (Ptr{simplex_data},), splx_data)
end

function print_levelling(z_b, gv, PP_ref_db, SS_ref_db)
    ccall((:print_levelling, libMAGEMin), Cvoid, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, gv, PP_ref_db, SS_ref_db)
end

function SS_UPDATE_function(gv, SS_ref_db, z_b, name)
    ccall((:SS_UPDATE_function, libMAGEMin), SS_ref, (global_variable, SS_ref, bulk_info, Ptr{Cchar}), gv, SS_ref_db, z_b, name)
end

function CP_UPDATE_function(gv, SS_ref_db, cp, z_b)
    ccall((:CP_UPDATE_function, libMAGEMin), csd_phase_set, (global_variable, SS_ref, csd_phase_set, bulk_info), gv, SS_ref_db, cp, z_b)
end

function split_cp(gv, SS_ref_db, cp)
    ccall((:split_cp, libMAGEMin), global_variable, (global_variable, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, SS_ref_db, cp)
end

function init_PGE_from_LP(gv, PC_read, SS_objective, z_b, SS_ref_db, cp)
    ccall((:init_PGE_from_LP, libMAGEMin), Cvoid, (global_variable, Ptr{PC_type}, Ptr{obj_type}, bulk_info, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, PC_read, SS_objective, z_b, SS_ref_db, cp)
end

function ss_min_PGE(gv, PC_read, SS_objective, NLopt_opt, z_b, SS_ref_db, cp)
    ccall((:ss_min_PGE, libMAGEMin), Cvoid, (global_variable, Ptr{PC_type}, Ptr{obj_type}, Ptr{NLopt_type}, bulk_info, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, PC_read, SS_objective, NLopt_opt, z_b, SS_ref_db, cp)
end

function ss_min_LP(gv, PC_read, SS_objective, NLopt_opt, z_b, SS_ref_db, cp)
    ccall((:ss_min_LP, libMAGEMin), Cvoid, (global_variable, Ptr{PC_type}, Ptr{obj_type}, Ptr{NLopt_type}, bulk_info, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, PC_read, SS_objective, NLopt_opt, z_b, SS_ref_db, cp)
end

function copy_to_Ppc(pc_check, add, ph_id, gv, SS_objective, SS_ref_db)
    ccall((:copy_to_Ppc, libMAGEMin), Cvoid, (Cint, Cint, Cint, global_variable, Ptr{obj_type}, Ptr{SS_ref}), pc_check, add, ph_id, gv, SS_objective, SS_ref_db)
end

function copy_to_cp(i, ph_id, gv, SS_ref_db, cp)
    ccall((:copy_to_cp, libMAGEMin), Cvoid, (Cint, Cint, global_variable, Ptr{SS_ref}, Ptr{csd_phase_set}), i, ph_id, gv, SS_ref_db, cp)
end

function init_ss_db(EM_database, z_b, gv, SS_ref_db)
    ccall((:init_ss_db, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{SS_ref}), EM_database, z_b, gv, SS_ref_db)
end

function init_ss_db_sb(EM_database, z_b, gv, SS_ref_db)
    ccall((:init_ss_db_sb, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{SS_ref}), EM_database, z_b, gv, SS_ref_db)
end

function print_help(gv)
    ccall((:print_help, libMAGEMin), Cvoid, (global_variable,), gv)
end

function retrieve_bulk_PT(gv, input_data, sgleP, z_b)
    ccall((:retrieve_bulk_PT, libMAGEMin), bulk_info, (global_variable, Ptr{io_data}, Cint, bulk_info), gv, input_data, sgleP, z_b)
end

function convert_system_comp(gv, sys_in, z_b)
    ccall((:convert_system_comp, libMAGEMin), Cvoid, (global_variable, Ptr{Cchar}, bulk_info), gv, sys_in, z_b)
end

function get_tests_bulks(gv)
    ccall((:get_tests_bulks, libMAGEMin), global_variable, (global_variable,), gv)
end

function vector_matrix_multiplication(v, M, result, vector_size, matrix_cols)
    ccall((:vector_matrix_multiplication, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Ptr{Cdouble}, Cint, Cint), v, M, result, vector_size, matrix_cols)
end

function matrix_vector_multiplication(M, v, result, matrix_rows, matrix_cols)
    ccall((:matrix_vector_multiplication, libMAGEMin), Cvoid, (Ptr{Ptr{Cdouble}}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint), M, v, result, matrix_rows, matrix_cols)
end

function get_act_sf_id(result, A, n)
    ccall((:get_act_sf_id, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Cint), result, A, n)
end

function inverseMatrix(ipiv, A1, n, work, lwork)
    ccall((:inverseMatrix, libMAGEMin), Cvoid, (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), ipiv, A1, n, work, lwork)
end

function MatMatMul(A, nrowA, B, ncolB, common, C)
    ccall((:MatMatMul, libMAGEMin), Cvoid, (Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Cint, Ptr{Ptr{Cdouble}}), A, nrowA, B, ncolB, common, C)
end

function VecMatMul(B1, A1, B, n)
    ccall((:VecMatMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), B1, A1, B, n)
end

function MatVecMul(A1, br, n_vec, n)
    ccall((:MatVecMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), A1, br, n_vec, n)
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

function print_1D_double_array(nx, array, title)
    ccall((:print_1D_double_array, libMAGEMin), Cvoid, (Cdouble, Ptr{Cdouble}, Ptr{Cchar}), nx, array, title)
end

function print_2D_double_array(nx, ny, array, title)
    ccall((:print_2D_double_array, libMAGEMin), Cvoid, (Cdouble, Cdouble, Ptr{Ptr{Cdouble}}, Ptr{Cchar}), nx, ny, array, title)
end

function print_1D_int_array(nx, array, title)
    ccall((:print_1D_int_array, libMAGEMin), Cvoid, (Cdouble, Ptr{Cint}, Ptr{Cchar}), nx, array, title)
end

function rnd(a)
    ccall((:rnd, libMAGEMin), Cdouble, (Cdouble,), a)
end

function SUPCRT_to_HSC(ElH, comp, size)
    ccall((:SUPCRT_to_HSC, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), ElH, comp, size)
end

function HSC_to_SUPCRT(ElH, comp, size)
    ccall((:HSC_to_SUPCRT, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), ElH, comp, size)
end

function norm_array(array, size)
    ccall((:norm_array, libMAGEMin), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), array, size)
end

function sum_norm_xipi(xi, pi, size)
    ccall((:sum_norm_xipi, libMAGEMin), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint), xi, pi, size)
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
    ccall((:print_cp, libMAGEMin), Cvoid, (global_variable, Ptr{csd_phase_set}), gv, cp)
end

function print_SS_informations(gv, SS_ref_db, iss)
    ccall((:print_SS_informations, libMAGEMin), Cvoid, (global_variable, SS_ref, Cint), gv, SS_ref_db, iss)
end

function rotate_hyperplane(gv, SS_ref_db)
    ccall((:rotate_hyperplane, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function non_rot_hyperplane(gv, SS_ref_db)
    ccall((:non_rot_hyperplane, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function raw_hyperplane(gv, SS_ref_db, gb)
    ccall((:raw_hyperplane, libMAGEMin), SS_ref, (global_variable, SS_ref, Ptr{Cdouble}), gv, SS_ref_db, gb)
end

function restrict_SS_HyperVolume(gv, SS_ref_db, box_size)
    ccall((:restrict_SS_HyperVolume, libMAGEMin), SS_ref, (global_variable, SS_ref, Cdouble), gv, SS_ref_db, box_size)
end

function check_SS_bounds(gv, SS_ref_db)
    ccall((:check_SS_bounds, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
end

function get_pp_id(gv)
    ccall((:get_pp_id, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_ss_id(gv, cp)
    ccall((:get_ss_id, libMAGEMin), global_variable, (global_variable, Ptr{csd_phase_set}), gv, cp)
end

function wave_melt_correction(gv, z_b, aspectRatio)
    ccall((:wave_melt_correction, libMAGEMin), global_variable, (global_variable, bulk_info, Cdouble), gv, z_b, aspectRatio)
end

function anelastic_correction(water, Vs0, P, T)
    ccall((:anelastic_correction, libMAGEMin), Cdouble, (Cint, Cdouble, Cdouble, Cdouble), water, Vs0, P, T)
end

function compute_phase_mol_fraction(gv, z_b, PP_ref_db, SS_ref_db, cp)
    ccall((:compute_phase_mol_fraction, libMAGEMin), global_variable, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, z_b, PP_ref_db, SS_ref_db, cp)
end

function compute_activities(EM_database, gv, PP_ref_db, z_b)
    ccall((:compute_activities, libMAGEMin), global_variable, (Cint, global_variable, Ptr{PP_ref}, bulk_info), EM_database, gv, PP_ref_db, z_b)
end

function compute_density_volume_modulus(EM_database, z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:compute_density_volume_modulus, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, z_b, gv, PP_ref_db, SS_ref_db, cp)
end

const NLOPT_MINF_MAX_REACHED = NLOPT_STOPVAL_REACHED

# Skipping MacroDefinition: NLOPT_DEPRECATED __attribute__ ( ( deprecated ) )

const n_ox_all = 16

# Skipping MacroDefinition: UTHASH_VERSION 2.1.0

const HASH_NONFATAL_OOM = 0

const HASH_INITIAL_NUM_BUCKETS = Cuint(32)

const HASH_INITIAL_NUM_BUCKETS_LOG2 = Cuint(5)

const HASH_BKT_CAPACITY_THRESH = Cuint(10)

const HASH_FCN = HASH_JEN

const HASH_BLOOM_BYTELEN = Cuint(0)

const HASH_SIGNATURE = Cuint(0xa0111fe1)

const HASH_BLOOM_SIGNATURE = Cuint(0xb12220f2)

const n_ox_mpf = 12

const n_ss_mpf = 18

const n_pp_mpf = 24

const n_ss_mp = 18

const n_pp_mp = 24

const n_ox_mp = 11

const ko_no_argument = 0

const ko_required_argument = 1

const ko_optional_argument = 2

#
# START OF EPILOGUE
#

# stable phases
struct SS_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    Comp_wt::Vector{Cdouble}
    Comp_apfu::Vector{Cdouble}
    compVariables::Vector{Cdouble}
    compVariablesNames::Vector{String}
    siteFractions::Vector{Cdouble}
    siteFractionsNames::Vector{String}
    emNames::Vector{String}
    emFrac::Vector{Cdouble}
    emFrac_wt::Vector{Cdouble}
    emChemPot::Vector{Cdouble}
    emComp::Vector{Vector{Float64}}
    emComp_wt::Vector{Vector{Float64}}
    emComp_apfu::Vector{Vector{Float64}}
end

function Base.convert(::Type{SS_data}, a::stb_SS_phases) 
    return SS_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.entropy, a.enthalpy, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                                    unsafe_wrap( Vector{Cdouble},        a.Comp,             a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.Comp_wt,          a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.Comp_apfu,        a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.compVariables,    a.n_xeos),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.compVariablesNames,a.n_xeos)),
                                    unsafe_wrap( Vector{Cdouble},        a.siteFractions,    a.n_sf),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.siteFractionsNames,a.n_sf)),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.emNames,          a.n_em)),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac,           a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac_wt,        a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emChemPot,        a.n_em),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp, a.n_em),  a.nOx),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp_wt, a.n_em),  a.nOx),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp_apfu, a.n_em),  a.nOx)   )
end

# metastable phases
struct mSS_data
    ph_name::String
    ph_type::String
    info::String
    ph_id::Cint
    em_id::Cint
    n_xeos::Cint
    n_em::Cint
    G_Ppc::Cdouble
    DF_Ppc::Cdouble
    comp_Ppc::Vector{Cdouble}
    p_Ppc::Vector{Cdouble}
    mu_Ppc::Vector{Cdouble}
    xeos_Ppc::Vector{Cdouble}
end

function Base.convert(::Type{mSS_data}, a::mstb_SS_phases) 
    return  mSS_data(   unsafe_string(a.ph_name),
                        unsafe_string(a.ph_type),
                        unsafe_string(a.info),
                        a.ph_id, a.em_id, a.n_xeos, a.n_em, a.G_Ppc, a.DF_Ppc,
                        unsafe_wrap( Vector{Cdouble},        a.comp_Ppc,           a.nOx),
                        unsafe_wrap( Vector{Cdouble},        a.p_Ppc,              a.n_em),
                        unsafe_wrap( Vector{Cdouble},        a.mu_Ppc,             a.n_em),
                        unsafe_wrap( Vector{Cdouble},        a.xeos_Ppc,           a.n_xeos)    )
end

# pure phases
struct PP_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    Comp_wt::Vector{Cdouble}
    Comp_apfu::Vector{Cdouble}
end

function Base.convert(::Type{PP_data}, a::stb_PP_phases) 
    return PP_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.entropy, a.enthalpy, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                    unsafe_wrap(Vector{Cdouble},a.Comp, a.nOx),
                    unsafe_wrap(Vector{Cdouble},a.Comp_wt, a.nOx),
                    unsafe_wrap(Vector{Cdouble},a.Comp_apfu, a.nOx))
end





#
# END OF EPILOGUE
#

# exports
const PREFIXES = ["CX", "clang_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
