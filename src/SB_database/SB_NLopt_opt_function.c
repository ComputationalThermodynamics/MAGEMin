/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
	Function to perform local minimization of the solution phases
	-> using NLopt external library
	-> penalty on SF inequality constraints seems good with a small value of 1e-12 (1e-10 was giving problems)                                         
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "nlopt.h"                  // requires specifying this in the makefile
#include "../MAGEMin.h"
#include "SB_NLopt_opt_function.h"
#include "../all_solution_phases.h"
#include "../toolkit.h"

// Equality constraint function: sum(x) == 1
double equality_constraint(unsigned n, const double *x, double *grad, void *data) {
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

SS_ref NLopt_opt_sb11_plg_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_plg, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_plg(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_sp_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_sp, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_sp(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_ol_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_ol, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_ol(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_wa_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_wa, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_wa(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_ri_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_ri, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_ri(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_opx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_opx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_opx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_cpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_cpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_cpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }

    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_hpcpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_hpcpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_hpcpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_ak_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_ak, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_ak(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_gtmj_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_gtmj, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_gtmj(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_pv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_pv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_pv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_ppv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_ppv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_ppv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_mw_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_mw, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_mw(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb11_cf_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb11_cf, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-7);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb11_cf(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};

void SB_sb11_NLopt_opt_init(        NLopt_type              *NLopt_opt,
                                    global_variable          gv         ){
                                        
    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_plg_function;              }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_sp_function;               }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_ol_function;               }
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_wa_function;               }
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_ri_function;               }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_opx_function;              }
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_cpx_function;              }
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_hpcpx_function;            }
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_ak_function;               }
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_gtmj_function;             }
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_pv_function;               }
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_ppv_function;              }
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_mw_function;               }
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb11_cf_function;               }
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }
    };
}


SS_ref NLopt_opt_sb21_plg_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_plg, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_plg(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_sp_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_sp, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_sp(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_ol_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_ol, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_ol(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_wa_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_wa, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_wa(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_ri_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_ri, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_ri(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_opx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_opx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_opx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_cpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_cpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_cpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_hpcpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_hpcpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_hpcpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_ak_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_ak, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_ak(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_gtmj_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_gtmj, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_gtmj(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_pv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_pv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_pv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_ppv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_ppv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_ppv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_cf_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_cf, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_cf(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_mw_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_mw, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_mw(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb21_nal_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb21_nal, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb21_nal(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};


void SB_sb21_NLopt_opt_init(        NLopt_type                  *NLopt_opt,
                        global_variable          gv                             ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_plg_function;              }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_sp_function;               }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_ol_function;               }
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_wa_function;               }
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_ri_function;               }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_opx_function;              }
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_cpx_function;              }
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_hpcpx_function;            }
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_ak_function;               }
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_gtmj_function;             }
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_pv_function;               }
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_ppv_function;              }
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_cf_function;               }
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_mw_function;               }
        else if (strcmp( gv.SS_list[iss], "nal")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb21_nal_function;              }
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }
    };
}

SS_ref NLopt_opt_sb24_plg_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_plg, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_plg(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_sp_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_sp, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_sp(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_ol_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_ol, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_ol(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_wa_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_wa, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_wa(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_ri_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_ri, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_ri(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_opx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_opx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_opx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_cpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_cpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_cpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_hpcpx_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_hpcpx, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_hpcpx(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_ak_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_ak, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_ak(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_gtmj_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_gtmj, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_gtmj(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_pv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_pv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_pv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_ppv_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_ppv, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_ppv(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_cf_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_cf, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_cf(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_mw_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_mw, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_mw(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};
SS_ref NLopt_opt_sb24_nal_function(global_variable gv, SS_ref SS_ref_db){
    unsigned int    n_em     = SS_ref_db.n_em;
    double *x  = SS_ref_db.iguess;
    for (int i = 0; i < (n_em); i++){
        SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];
        SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];
    }
    SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); 
    nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);
    nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);
    nlopt_set_min_objective(SS_ref_db.opt, obj_sb24_nal, &SS_ref_db);
    nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);
    nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);
    nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);
    double minf;
    if (gv.maxeval==1){  
        minf = obj_sb24_nal(n_em, x, NULL, &SS_ref_db);
    }
    else{
        SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);
    }
    /* Send back needed local solution parameters */
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.xeos[i] = x[i];
    }

    SS_ref_db.df   = minf;
    nlopt_destroy(SS_ref_db.opt);

    return SS_ref_db;
};


void SB_sb24_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,
                        global_variable 	 gv				){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_plg_function; 		}
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_sp_function; 		}
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_ol_function; 		}
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_wa_function; 		}
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_ri_function; 		}
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_opx_function; 		}
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_cpx_function; 		}
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_hpcpx_function; 		}
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_ak_function; 		}
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_gtmj_function; 		}
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_pv_function; 		}
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_ppv_function; 		}
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_cf_function; 		}
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_mw_function; 		}
        else if (strcmp( gv.SS_list[iss], "nal")  == 0 ){
            NLopt_opt[iss]  = NLopt_opt_sb24_nal_function; 		}
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }	
    };
}

void SB_NLopt_opt_init(	        	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				){


	if (gv.EM_database == 0){				// metapelite database //
		SB_sb11_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 1){				// metapelite database //
		SB_sb21_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
	else if (gv.EM_database == 2){				// stixrude 2024 database //
		SB_sb24_NLopt_opt_init(	 				NLopt_opt,
												gv							);
	}
}

