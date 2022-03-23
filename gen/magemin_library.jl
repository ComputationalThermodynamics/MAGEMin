module LibMAGEMin

using MAGEMin_jll
export MAGEMin_jll

using CEnum

#
# START OF PROLOGUE
#


const HASH_JEN = 0;

if isfile("libMAGEMin.dylib")
    libMAGEMin = joinpath(pwd(),"libMAGEMin.dylib")
    println("Loading local libMAGEMin dynamic library")
else
    using MAGEMin_jll
    export MAGEMin_jll
    println("Loading libMAGEMin dynamic library from MAGEMin_jll")
end




#
# END OF PROLOGUE
#

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

function sum_array(array, size)
    ccall((:sum_array, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, size)
end

function check_sign(v1, v2)
    ccall((:check_sign, libMAGEMin), Cint, (Cdouble, Cdouble), v1, v2)
end

struct PP_refs
    Name::NTuple{20, Cchar}
    Comp::NTuple{11, Cdouble}
    gbase::Cdouble
    gb_lvl::Cdouble
    factor::Cdouble
    phase_density::Cdouble
    phase_shearModulus::Cdouble
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_bulkModulus::Cdouble
    volume::Cdouble
    mass::Cdouble
end

const PP_ref = PP_refs

function G_EM_function(EM_database, bulk_rock, P, T, name, state)
    ccall((:G_EM_function, libMAGEMin), PP_ref, (Cint, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cchar}), EM_database, bulk_rock, P, T, name, state)
end

function runMAGEMin(argc, argv)
    ccall((:runMAGEMin, libMAGEMin), Cint, (Cint, Ptr{Ptr{Cchar}}), argc, argv)
end

function find_EM_id(em_tag)
    ccall((:find_EM_id, libMAGEMin), Cint, (Ptr{Cchar},), em_tag)
end

function norm_array(array, size)
    ccall((:norm_array, libMAGEMin), Ptr{Cdouble}, (Ptr{Cdouble}, Cint), array, size)
end

function get_bulk(bulk_rock, test, n_El)
    ccall((:get_bulk, libMAGEMin), Cvoid, (Ptr{Cdouble}, Cint, Cint), bulk_rock, test, n_El)
end

mutable struct EM_db
    Name::NTuple{20, Cchar}
    Comp::NTuple{12, Cdouble}
    input_1::NTuple{3, Cdouble}
    input_2::NTuple{4, Cdouble}
    input_3::NTuple{11, Cdouble}
    input_4::NTuple{3, Cdouble}
    EM_db() = new()
end

function Access_EM_DB(id, EM_database)
    ccall((:Access_EM_DB, libMAGEMin), EM_db, (Cint, Cint), id, EM_database)
end

function get_EM_DB_names(EM_database)
    ccall((:get_EM_DB_names, libMAGEMin), Ptr{Ptr{Cchar}}, (Cint,), EM_database)
end

mutable struct bulk_info
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
    bulk_info() = new()
end

function zeros_in_bulk(bulk_rock, P, T)
    ccall((:zeros_in_bulk, libMAGEMin), bulk_info, (Ptr{Cdouble}, Cdouble, Cdouble), bulk_rock, P, T)
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
    bounds::Ptr{Ptr{Cdouble}}
    bounds_ref::Ptr{Ptr{Cdouble}}
    z_em::Ptr{Cdouble}
    n_guess::Cint
    iguess::Ptr{Cdouble}
    dguess::Ptr{Cdouble}
    check_df::Cdouble
    forced_stop::Cint
    xeos_sf_ok_saved::Cint
    status::Cint
    nlopt_verb::Cint
    tol_sf::Ptr{Cdouble}
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
    ElShearMod::Ptr{Cdouble}
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
    bulk::Ptr{Cdouble}
    in_gam::Ptr{Cdouble}
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
    phase_bulkModulus::Cdouble
    phase_shearModulus::Cdouble
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
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    n_xeos::Cint
    n_em::Cint
    Comp::Ptr{Cdouble}
    compVariables::Ptr{Cdouble}
    emNames::Ptr{Ptr{Cchar}}
    emFrac::Ptr{Cdouble}
    emChemPot::Ptr{Cdouble}
    emComp::Ptr{Ptr{Cdouble}}
end

const stb_SS_phase = stb_SS_phases

struct stb_PP_phases
    nOx::Cint
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Ptr{Cdouble}
end

const stb_PP_phase = stb_PP_phases

struct stb_systems
    MAGEMin_ver::Ptr{Cchar}
    nOx::Cint
    oxides::Ptr{Ptr{Cchar}}
    P::Cdouble
    T::Cdouble
    bulk::Ptr{Cdouble}
    gamma::Ptr{Cdouble}
    G::Cdouble
    bulk_res_norm::Cdouble
    rho::Cdouble
    bulk_S::Ptr{Cdouble}
    frac_S::Cdouble
    rho_S::Cdouble
    bulk_M::Ptr{Cdouble}
    frac_M::Cdouble
    rho_M::Cdouble
    bulk_F::Ptr{Cdouble}
    frac_F::Cdouble
    rho_F::Cdouble
    n_ph::Cint
    n_PP::Cint
    n_SS::Cint
    ph::Ptr{Ptr{Cchar}}
    ph_frac::Ptr{Cdouble}
    ph_type::Ptr{Cint}
    ph_id::Ptr{Cint}
    SS::Ptr{stb_SS_phase}
    PP::Ptr{stb_PP_phase}
end

const stb_system = stb_systems

mutable struct global_variables
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
    system_density::Cdouble
    system_bulkModulus::Cdouble
    system_shearModulus::Cdouble
    system_Vp::Cdouble
    system_Vs::Cdouble
    global_variables() = new()
end

const global_variable = global_variables

function global_variable_init()
    ccall((:global_variable_init, libMAGEMin), global_variable, ())
end

mutable struct Database
    PP_ref_db::Ptr{PP_ref}
    SS_ref_db::Ptr{SS_ref}
    cp::Ptr{csd_phase_set}
    sp::Ptr{stb_system}
    EM_names::Ptr{Ptr{Cchar}}
    Database() = new()
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
    ccall((:ComputePostProcessing, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), EM_database, z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function ReadCommandLineOptions(gv, argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, maxeval_out, get_version_out)
    ccall((:ReadCommandLineOptions, libMAGEMin), global_variable, (global_variable, Cint, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), gv, argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, maxeval_out, get_version_out)
end

function PrintOutput(gv, rank, l, DB, time_taken, z_b)
    ccall((:PrintOutput, libMAGEMin), Cvoid, (global_variable, Cint, Cint, Databases, Cdouble, bulk_info), gv, rank, l, DB, time_taken, z_b)
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

@cenum var"##Ctag#302"::UInt32 begin
    _tc_ds634_ = 0
end

mutable struct EM2id
    EM_tag::NTuple{20, Cchar}
    id::Cint
    hh::UT_hash_handle
    EM2id() = new()
end

mutable struct PP2id
    PP_tag::NTuple{20, Cchar}
    id::Cint
    hh::UT_hash_handle
    PP2id() = new()
end

function find_PP_id(PP_tag)
    ccall((:find_PP_id, libMAGEMin), Cint, (Ptr{Cchar},), PP_tag)
end

# typedef void ( * sf_type ) ( unsigned m , double * result , unsigned n , const double * x , double * grad , void * data )
const sf_type = Ptr{Cvoid}

function NLopt_global_opt_function(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:NLopt_global_opt_function, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function NLopt_opt_function(gv, SS_ref_db, index)
    ccall((:NLopt_opt_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint), gv, SS_ref_db, index)
end

function PGE(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:PGE, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

function norm_vector(array, n)
    ccall((:norm_vector, libMAGEMin), Cdouble, (Ptr{Cdouble}, Cint), array, n)
end

struct ss_pc
    xeos_pc::NTuple{11, Cdouble}
end

mutable struct PC_refs
    ss_pc_xeos::Ptr{ss_pc}
    PC_refs() = new()
end

const PC_ref = PC_refs

function SS_PC_init_function(SS_PC_xeos, iss, name)
    ccall((:SS_PC_init_function, libMAGEMin), Cvoid, (Ptr{PC_ref}, Cint, Ptr{Cchar}), SS_PC_xeos, iss, name)
end

function dump_init(gv)
    ccall((:dump_init, libMAGEMin), Cvoid, (global_variable,), gv)
end

function fill_output_struct(gv, z_b, PP_ref_db, SS_ref_db, cp, sp)
    ccall((:fill_output_struct, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}, Ptr{stb_system}), gv, z_b, PP_ref_db, SS_ref_db, cp, sp)
end

function dump_results_function(gv, z_b, PP_ref_db, SS_ref_db, cp)
    ccall((:dump_results_function, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), gv, z_b, PP_ref_db, SS_ref_db, cp)
end

function mergeParallelFiles(gv)
    ccall((:mergeParallelFiles, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_LocalMinima_Files(gv)
    ccall((:mergeParallel_LocalMinima_Files, libMAGEMin), Cvoid, (global_variable,), gv)
end

function mergeParallel_LevellingGamma_Files(gv)
    ccall((:mergeParallel_LevellingGamma_Files, libMAGEMin), Cvoid, (global_variable,), gv)
end

function G_SS_EM_function(gv, SS_ref_db, EM_database, z_b, name)
    ccall((:G_SS_EM_function, libMAGEMin), SS_ref, (global_variable, SS_ref, Cint, bulk_info, Ptr{Cchar}), gv, SS_ref_db, EM_database, z_b, name)
end

function G_SS_INIT_EM_function(SS_ref_db, EM_database, name, gv)
    ccall((:G_SS_INIT_EM_function, libMAGEMin), SS_ref, (SS_ref, Cint, Ptr{Cchar}, global_variable), SS_ref_db, EM_database, name, gv)
end

function CP_INIT_function(cp, gv)
    ccall((:CP_INIT_function, libMAGEMin), csd_phase_set, (csd_phase_set, global_variable), cp, gv)
end

function SP_INIT_function(sp, gv)
    ccall((:SP_INIT_function, libMAGEMin), stb_system, (stb_system, global_variable), sp, gv)
end

function read_in_data(gv, input_data, file_name, n_points)
    ccall((:read_in_data, libMAGEMin), Cvoid, (global_variable, Ptr{io_data}, Ptr{Cchar}, Cint), gv, input_data, file_name, n_points)
end

function AddResults_output_struct(gv, z_b, P, T, DB, output)
    ccall((:AddResults_output_struct, libMAGEMin), Cvoid, (global_variable, bulk_info, Cdouble, Cdouble, Databases, out_data), gv, z_b, P, T, DB, output)
end

function InitializeOutput(gv, DB)
    ccall((:InitializeOutput, libMAGEMin), out_data, (global_variable, Databases), gv, DB)
end

function FreeOutput(output)
    ccall((:FreeOutput, libMAGEMin), Cvoid, (out_data,), output)
end

mutable struct ketopt_t
    ind::Cint
    opt::Cint
    arg::Ptr{Cchar}
    longidx::Cint
    i::Cint
    pos::Cint
    n_args::Cint
    ketopt_t() = new()
end

mutable struct ko_longopt_t
    name::Ptr{Cchar}
    has_arg::Cint
    val::Cint
    ko_longopt_t() = new()
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
    ccall((:SS_objective_init_function, libMAGEMin), Cvoid, (Ptr{obj_type}, global_variable), SS_objective, gv)
end

function p2x_bi(SS_ref_db, eps)
    ccall((:p2x_bi, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_cd(SS_ref_db, eps)
    ccall((:p2x_cd, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_cpx(SS_ref_db, eps)
    ccall((:p2x_cpx, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_ep(SS_ref_db, eps)
    ccall((:p2x_ep, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_fl(SS_ref_db, eps)
    ccall((:p2x_fl, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_g(SS_ref_db, eps)
    ccall((:p2x_g, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_hb(SS_ref_db, eps)
    ccall((:p2x_hb, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_ilm(SS_ref_db, eps)
    ccall((:p2x_ilm, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_liq(SS_ref_db, eps)
    ccall((:p2x_liq, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_mu(SS_ref_db, eps)
    ccall((:p2x_mu, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_ol(SS_ref_db, eps)
    ccall((:p2x_ol, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_opx(SS_ref_db, eps)
    ccall((:p2x_opx, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_pl4T(SS_ref_db, eps)
    ccall((:p2x_pl4T, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
end

function p2x_spn(SS_ref_db, eps)
    ccall((:p2x_spn, libMAGEMin), Cvoid, (SS_ref, Cdouble), SS_ref_db, eps)
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

function PC_PX_function(SS_ref_db, x, name)
    ccall((:PC_PX_function, libMAGEMin), SS_ref, (SS_ref, Ptr{Cdouble}, Ptr{Cchar}), SS_ref_db, x, name)
end

function PC_function(gv, SS_ref_db, z_b, name)
    ccall((:PC_function, libMAGEMin), SS_ref, (global_variable, SS_ref, bulk_info, Ptr{Cchar}), gv, SS_ref_db, z_b, name)
end

function P2X(gv, SS_ref_db, z_b, name)
    ccall((:P2X, libMAGEMin), SS_ref, (global_variable, SS_ref, bulk_info, Ptr{Cchar}), gv, SS_ref_db, z_b, name)
end

function get_phase_id(gv, name)
    ccall((:get_phase_id, libMAGEMin), Cint, (global_variable, Ptr{Cchar}), gv, name)
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

function inverseMatrix(A1, n)
    ccall((:inverseMatrix, libMAGEMin), Cvoid, (Ptr{Cdouble}, Cint), A1, n)
end

function VecMatMul(B1, A1, B, n)
    ccall((:VecMatMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), B1, A1, B, n)
end

function MatVecMul(A1, br, n_vec, n)
    ccall((:MatVecMul, libMAGEMin), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), A1, br, n_vec, n)
end

function Levelling(z_b, gv, PP_ref_db, SS_ref_db, cp)
    ccall((:Levelling, libMAGEMin), global_variable, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}, Ptr{csd_phase_set}), z_b, gv, PP_ref_db, SS_ref_db, cp)
end

mutable struct simplex_datas
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
    simplex_datas() = new()
end

const simplex_data = simplex_datas

function print_levelling(z_b, gv, PP_ref_db, SS_ref_db)
    ccall((:print_levelling, libMAGEMin), Cvoid, (bulk_info, global_variable, Ptr{PP_ref}, Ptr{SS_ref}), z_b, gv, PP_ref_db, SS_ref_db)
end

function SS_UPDATE_function(gv, SS_ref_db, z_b, name)
    ccall((:SS_UPDATE_function, libMAGEMin), SS_ref, (global_variable, SS_ref, bulk_info, Ptr{Cchar}), gv, SS_ref_db, z_b, name)
end

function CP_UPDATE_function(gv, SS_ref_db, cp, z_b)
    ccall((:CP_UPDATE_function, libMAGEMin), csd_phase_set, (global_variable, SS_ref, csd_phase_set, bulk_info), gv, SS_ref_db, cp, z_b)
end

function split_cp(i, gv, SS_ref_db, cp)
    ccall((:split_cp, libMAGEMin), global_variable, (Cint, global_variable, Ptr{SS_ref}, Ptr{csd_phase_set}), i, gv, SS_ref_db, cp)
end

function ss_min_PGE(mode, i, gv, z_b, SS_ref_db, cp)
    ccall((:ss_min_PGE, libMAGEMin), Cvoid, (Cint, Cint, global_variable, bulk_info, Ptr{SS_ref}, Ptr{csd_phase_set}), mode, i, gv, z_b, SS_ref_db, cp)
end

function reset_SS(gv, z_b, SS_ref_db)
    ccall((:reset_SS, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{SS_ref}), gv, z_b, SS_ref_db)
end

function reset_gv(gv, z_b, PP_ref_db, SS_ref_db)
    ccall((:reset_gv, libMAGEMin), global_variable, (global_variable, bulk_info, Ptr{PP_ref}, Ptr{SS_ref}), gv, z_b, PP_ref_db, SS_ref_db)
end

function reset_cp(gv, z_b, cp)
    ccall((:reset_cp, libMAGEMin), Cvoid, (global_variable, bulk_info, Ptr{csd_phase_set}), gv, z_b, cp)
end

function reset_sp(gv, sp)
    ccall((:reset_sp, libMAGEMin), Cvoid, (global_variable, Ptr{stb_system}), gv, sp)
end

function init_ss_db(EM_database, z_b, gv, SS_ref_db)
    ccall((:init_ss_db, libMAGEMin), global_variable, (Cint, bulk_info, global_variable, Ptr{SS_ref}), EM_database, z_b, gv, SS_ref_db)
end

function SS_ref_destroy(gv, SS_ref_db)
    ccall((:SS_ref_destroy, libMAGEMin), Cvoid, (global_variable, Ptr{SS_ref}), gv, SS_ref_db)
end

function CP_destroy(gv, cp)
    ccall((:CP_destroy, libMAGEMin), Cvoid, (global_variable, Ptr{csd_phase_set}), gv, cp)
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
    ccall((:print_cp, libMAGEMin), Cvoid, (global_variable, Ptr{csd_phase_set}), gv, cp)
end

function print_SS_informations(gv, SS_ref_db, iss)
    ccall((:print_SS_informations, libMAGEMin), Cvoid, (global_variable, SS_ref, Cint), gv, SS_ref_db, iss)
end

function rotate_hyperplane(gv, SS_ref_db)
    ccall((:rotate_hyperplane, libMAGEMin), SS_ref, (global_variable, SS_ref), gv, SS_ref_db)
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

mutable struct TMatrix
    m::Ptr{Ptr{Cdouble}}
    nRows::Cint
    nCols::Cint
    TMatrix() = new()
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

function get_pp_id(gv)
    ccall((:get_pp_id, libMAGEMin), global_variable, (global_variable,), gv)
end

function get_ss_id(gv, cp)
    ccall((:get_ss_id, libMAGEMin), global_variable, (global_variable, Ptr{csd_phase_set}), gv, cp)
end

const n_em_db = 291

const nEl = 11

const NLOPT_MINF_MAX_REACHED = NLOPT_STOPVAL_REACHED

# Skipping MacroDefinition: NLOPT_DEPRECATED __attribute__ ( ( deprecated ) )

# Skipping MacroDefinition: UTHASH_VERSION 2.1.0

const HASH_NONFATAL_OOM = 0

const HASH_INITIAL_NUM_BUCKETS = Cuint(32)

const HASH_INITIAL_NUM_BUCKETS_LOG2 = Cuint(5)

const HASH_BKT_CAPACITY_THRESH = Cuint(10)

const HASH_FCN = HASH_JEN

const HASH_BLOOM_BYTELEN = Cuint(0)

const HASH_SIGNATURE = Cuint(0xa0111fe1)

const HASH_BLOOM_SIGNATURE = Cuint(0xb12220f2)

const ko_no_argument = 0

const ko_required_argument = 1

const ko_optional_argument = 2

#
# START OF EPILOGUE
#




struct SS_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    compVariables::Vector{Cdouble}
    emNames::Vector{String}
    emFrac::Vector{Cdouble}
    emChemPot::Vector{Cdouble}
    emComp::Vector{Vector{Float64}}
end



function Base.convert(::Type{SS_data}, a::stb_SS_phases) 
    return SS_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                                    unsafe_wrap( Vector{Cdouble},        a.Comp,             a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.compVariables,    a.n_xeos),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.emNames,          a.n_em)),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac,           a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emChemPot,        a.n_em),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp, a.n_em),  a.nOx)   )
end

struct PP_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
end

function Base.convert(::Type{PP_data}, a::stb_PP_phases) 
    return PP_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                    unsafe_wrap(Vector{Cdouble},a.Comp, a.nOx))
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
