/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <string.h>

#include "nlopt.h"
#include "../MAGEMin.h"
#include "GH_NLopt_opt_function.h"
#include "../all_solution_phases.h"

/**
    Same NLopt setup as SB_database's solution phases (e.g.
    NLopt_opt_sb11_ol_function): n = n_em = n_xeos, box bounds, and a
    single linear equality constraint enforcing Sigma(x) = 1 - matching
    "gh"'s LP-only solving path (see SetupDatabase), which reuses the
    same shared LP/PGE machinery as "sb".
*/
static double gh_equality_constraint(unsigned n, const double *x, double *grad, void *data) {
    if (grad) {
        for (unsigned i = 0; i < n; i++) {
            grad[i] = 1.0;
        }
    }
    double sum = 0.0;
    for (unsigned i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum - 1.0;
}

SS_ref NLopt_opt_gh_liq_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int n_em = SS_ref_db.n_em;
    double *x = SS_ref_db.iguess;

    for (int i = 0; i < (int)n_em; i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }

    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, n_em);
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_gh_liq, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, gh_equality_constraint, NULL, 1e-5);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);

    double minf;
    if (gv.maxeval==1){
        minf = obj_gh_liq(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }

    for (int i = 0; i < (int)SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};

/**
    Shared NLopt setup for the 4 "sb-trivial" mineral solid solutions:
    identical structure to NLopt_opt_gh_liq_function, just parameterized
    over the phase's own objective function, since all 4 share the same
    "n = n_em = n_xeos, box bounds, single Sigma(x)=1 equality constraint"
    formulation.
*/
static SS_ref GH_generic_NLopt_opt_function(global_variable gv, SS_ref SS_ref_db, double (*obj)(unsigned, const double*, double*, void*)){
    unsigned int n_em = SS_ref_db.n_em;
    double *x = SS_ref_db.iguess;

    for (int i = 0; i < (int)n_em; i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }

    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, n_em);
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, gh_equality_constraint, NULL, 1e-5);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);

    double minf;
    if (gv.maxeval==1){
        minf = obj(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }

    for (int i = 0; i < (int)SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
}

SS_ref NLopt_opt_gh_ol_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_ol);
}
SS_ref NLopt_opt_gh_fsp_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_fsp);
}
SS_ref NLopt_opt_gh_bi_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_bi);
}
SS_ref NLopt_opt_gh_g_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_g);
}
SS_ref NLopt_opt_gh_hb_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_hb);
}
SS_ref NLopt_opt_gh_lc_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_lc);
}
SS_ref NLopt_opt_gh_mel_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_mel);
}
SS_ref NLopt_opt_gh_cum_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_cum);
}
SS_ref NLopt_opt_gh_spn_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_spn);
}
SS_ref NLopt_opt_gh_cpx_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_cpx);
}
SS_ref NLopt_opt_gh_opx_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_opx);
}
SS_ref NLopt_opt_gh_fluid_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_fluid);
}
SS_ref NLopt_opt_gh_rhm_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_rhm);
}
SS_ref NLopt_opt_gh_nph_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_nph);
}
SS_ref NLopt_opt_gh_kls_function(global_variable gv, SS_ref SS_ref_db){
    return GH_generic_NLopt_opt_function(gv, SS_ref_db, obj_gh_kls);
}

void GH_NLopt_opt_init(             NLopt_type          *NLopt_opt,
                                    global_variable      gv              ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_liq_function;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_ol_function;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_fsp_function;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_bi_function;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_g_function;
        }
        else if (strcmp( gv.SS_list[iss], "hb") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_hb_function;
        }
        else if (strcmp( gv.SS_list[iss], "lc") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_lc_function;
        }
        else if (strcmp( gv.SS_list[iss], "mel") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_mel_function;
        }
        else if (strcmp( gv.SS_list[iss], "cum") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_cum_function;
        }
        else if (strcmp( gv.SS_list[iss], "spn") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_spn_function;
        }
        else if (strcmp( gv.SS_list[iss], "cpx") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_cpx_function;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_opx_function;
        }
        else if (strcmp( gv.SS_list[iss], "fl") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_fluid_function;
        }
        else if (strcmp( gv.SS_list[iss], "rhm") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_rhm_function;
        }
        else if (strcmp( gv.SS_list[iss], "nph") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_nph_function;
        }
        else if (strcmp( gv.SS_list[iss], "kls") == 0 ){
            NLopt_opt[iss] = NLopt_opt_gh_kls_function;
        }
    }
}
