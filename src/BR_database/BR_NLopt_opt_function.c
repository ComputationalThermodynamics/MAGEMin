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

static SS_ref NLopt_opt_br_generic(global_variable gv, SS_ref SS_ref_db, obj_type obj){
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
    nlopt_add_equality_constraint(SS_ref_db.opt, br_equality_constraint, NULL, 1e-5);
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

SS_ref NLopt_opt_br_ctd_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_ctd);
}

SS_ref NLopt_opt_br_car_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_car);
}

/* obj_br_chl/obj_br_mica now carry full analytic gradients (derived via the
   site-fraction chain rule, same structural approach as TC's own dp_dx
   matrices) - SLSQP throughout, same as every other BR solution phase. */
SS_ref NLopt_opt_br_chl_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_chl);
}

SS_ref NLopt_opt_br_mica_function(global_variable gv, SS_ref SS_ref_db){
    return NLopt_opt_br_generic(gv, SS_ref_db, obj_br_mica);
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
