module LibMAGEMin

using MAGEMin_jll
export MAGEMin_jll

using CEnum

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
    phase_cp::Cdouble
    phase_expansivity::Cdouble
    phase_shearModulus::Cdouble
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

function global_variable_init()
    ccall((:global_variable_init, libMAGEMin), global_variable, ())
end

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

function ReadCommandLineOptions(argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, n_pc_out, maxeval_out)
    ccall((:ReadCommandLineOptions, libMAGEMin), Cint, (Cint, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), argc, argv, Mode_out, Verb_out, test_out, n_points_out, P, T, Bulk, Gam, InitEM_Prop, File, Phase, n_pc_out, maxeval_out)
end

function PrintOutput(gv, rank, l, DB, time_taken, z_b)
    ccall((:PrintOutput, libMAGEMin), Cvoid, (global_variable, Cint, Cint, Databases, Cdouble, bulk_info), gv, rank, l, DB, time_taken, z_b)
end

const n_em_db = 291

const nEl = 11

const NLOPT_MINF_MAX_REACHED = NLOPT_STOPVAL_REACHED

# Skipping MacroDefinition: NLOPT_DEPRECATED __attribute__ ( ( deprecated ) )

# exports
const PREFIXES = ["CX", "clang_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
