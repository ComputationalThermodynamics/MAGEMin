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
#ifndef __BR_OBJECTIVE_FUNCTIONS_H_
#define __BR_OBJECTIVE_FUNCTIONS_H_

#include "../MAGEMin.h"

double obj_br_ctd(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_car(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_chl(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_mica(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_talc(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_bt(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_ol(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_ep(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_amph(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_spl(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_stau(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_crd(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_grt(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_omph(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_amphx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_br_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db);

void BR_SS_objective_init_function(obj_type *SS_objective, global_variable gv);

void BR_PC_init(PC_type *PC_read, global_variable gv);

SS_ref BR_PC_function(global_variable gv, PC_type *PC_read, SS_ref SS_ref_db, bulk_info z_b, int ph_id);

typedef void (*P2X_type) (				void 				*SS_ref_db,
										double 				 eps				);

void p2x_br_generic(void *SS_ref_db, double eps);

void BR_P2X_init(P2X_type *P2X_read, global_variable gv);

#endif
