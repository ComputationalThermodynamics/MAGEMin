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
#include "BR_NLopt_opt_function.h"
#include "br_objective_functions.h"

static double br_equality_constraint(unsigned n, const double *x, double *grad, void *data) {
    (void) data;
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

typedef struct br_sf_lin_ {
    int     n_sf;
    double  sf0[16];
    double  J[16*16];
} br_sf_lin;

static void br_sf_constraint(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data){
    br_sf_lin *L = (br_sf_lin *) data;
    for (unsigned i = 0; i < m; i++){
        double r = L->sf0[i];
        for (unsigned k = 0; k < n; k++){ r += L->J[i*n + k]*x[k]; }
        result[i] = -r;
        if (grad){
            for (unsigned k = 0; k < n; k++){ grad[i*n + k] = -L->J[i*n + k]; }
        }
    }
}

static void br_sf_linearize(SS_ref *d, obj_type obj, br_sf_lin *L){
    int     n  = d->n_xeos;
    double  z[n], s1[d->n_sf];
    L->n_sf = d->n_sf;
    for (int k = 0; k < n; k++){
        for (int j = 0; j < n; j++){ z[j] = (j == k) ? 1.0 : 0.0; }
        obj(n, z, NULL, d);
        for (int i = 0; i < L->n_sf; i++){ s1[i] = d->sf[i]; }
        z[k] = 2.0;
        obj(n, z, NULL, d);
        for (int i = 0; i < L->n_sf; i++){
            L->J[i*n + k] = d->sf[i] - s1[i];
            if (k == 0){ L->sf0[i] = 2.0*s1[i] - d->sf[i]; }
        }
    }
}

static SS_ref NLopt_opt_br_generic_sf(global_variable gv, SS_ref SS_ref_db, obj_type obj, int sf_constrained){
    unsigned int n_em = SS_ref_db.n_em;
    double *x = SS_ref_db.iguess;

    for (int i = 0; i < (int)n_em; i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }

    /* SLSQP cold-starting exactly on a simplex vertex (all-but-one x[i] at
       their lower bound) hits a genuine internal QP-subproblem singularity
       in NLopt - reproducibly for fsp near pure-albite, reported as
       nlopt_optimize returning -1 (NLOPT_FAILURE) even though the gradient
       there is analytically correct (FD-verified). Standard mitigation:
       nudge any component sitting within 1e-3 of a bound inward before the
       solve, then renormalize to keep the equality constraint satisfied. */
    {
        double sum_free = 0.0;
        double sum_pin  = 0.0;
        int nudged = 0;
        for (int i = 0; i < (int)n_em; i++){
            double lb = SS_ref_db.lb[i], ub = SS_ref_db.ub[i];
            if (ub - lb <= 2e-3){
                if (x[i] != lb){ x[i] = lb; nudged = 1; }
                sum_pin += x[i];
                continue;
            }
            if (x[i] - lb < 1e-3){ x[i] = lb + 1e-3; nudged = 1; }
            else if (ub - x[i] < 1e-3){ x[i] = ub - 1e-3; nudged = 1; }
            sum_free += x[i];
        }
        if (nudged && sum_free > 0.0){
            for (int i = 0; i < (int)n_em; i++){
                if (SS_ref_db.ub[i] - SS_ref_db.lb[i] > 2e-3){ x[i] *= (1.0 - sum_pin)/sum_free; }
            }
        }
    }

    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, n_em);
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, br_equality_constraint, NULL, 1e-5);
    br_sf_lin L;
    if (sf_constrained){
        double tol[SS_ref_db.n_sf];
        for (int i = 0; i < SS_ref_db.n_sf; i++){ tol[i] = 0.0; }
        br_sf_linearize(&SS_ref_db, obj, &L);
        nlopt_add_inequality_mconstraint(SS_ref_db.opt, SS_ref_db.n_sf, br_sf_constraint, &L, tol);
    }
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

static SS_ref NLopt_opt_br_generic(global_variable gv, SS_ref SS_ref_db, obj_type obj){
    return NLopt_opt_br_generic_sf(gv, SS_ref_db, obj, 0);
}

SS_ref NLopt_opt_br_ctd_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_ctd);
}

SS_ref NLopt_opt_br_car_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_car);
}

SS_ref NLopt_opt_br_chl_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic_sf(gv, SS_ref_db, obj_br_chl, 1);
}

SS_ref NLopt_opt_br_mica_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic_sf(gv, SS_ref_db, obj_br_mica, 1);
}

SS_ref NLopt_opt_br_talc_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_talc);
}

SS_ref NLopt_opt_br_ilm_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_ilm);
}

SS_ref NLopt_opt_br_bt_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_bt);
}

SS_ref NLopt_opt_br_ol_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_ol);
}

SS_ref NLopt_opt_br_ep_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_ep);
}

SS_ref NLopt_opt_br_opx_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_opx);
}

SS_ref NLopt_opt_br_amph_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_amph);
}

SS_ref NLopt_opt_br_spl_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_spl);
}

SS_ref NLopt_opt_br_stau_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_stau);
}

SS_ref NLopt_opt_br_crd_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_crd);
}

SS_ref NLopt_opt_br_grt_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_grt);
}

SS_ref NLopt_opt_br_omph_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_omph);
}

SS_ref NLopt_opt_br_amphx_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_amphx);
}

SS_ref NLopt_opt_br_fsp_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_fsp);
}

void BR_NLopt_opt_init(NLopt_type *NLopt_opt, global_variable gv){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "ctd") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_ctd_function;
        }
        else if (strcmp( gv.SS_list[iss], "car") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_car_function;
        }
        else if (strcmp( gv.SS_list[iss], "chl") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_chl_function;
        }
        else if (strcmp( gv.SS_list[iss], "mica") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_mica_function;
        }
        else if (strcmp( gv.SS_list[iss], "talc") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_talc_function;
        }
        else if (strcmp( gv.SS_list[iss], "ilm") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_ilm_function;
        }
        else if (strcmp( gv.SS_list[iss], "bt") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_bt_function;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_ol_function;
        }
        else if (strcmp( gv.SS_list[iss], "ep") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_ep_function;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_opx_function;
        }
        else if (strcmp( gv.SS_list[iss], "amph") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_amph_function;
        }
        else if (strcmp( gv.SS_list[iss], "spl") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_spl_function;
        }
        else if (strcmp( gv.SS_list[iss], "stau") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_stau_function;
        }
        else if (strcmp( gv.SS_list[iss], "crd") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_crd_function;
        }
        else if (strcmp( gv.SS_list[iss], "grt") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_grt_function;
        }
        else if (strcmp( gv.SS_list[iss], "omph") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_omph_function;
        }
        else if (strcmp( gv.SS_list[iss], "amphx") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_amphx_function;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            NLopt_opt[iss] = NLopt_opt_br_fsp_function;
        }
    }
}
