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
#ifndef __GH_OBJECTIVE_FUNCTIONS_H_
#define __GH_OBJECTIVE_FUNCTIONS_H_

#include "../MAGEMin.h"

double obj_gh_liq(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_ol(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_bi(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_g(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_hb(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_lc(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_mel(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_cum(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_spn(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_fluid(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_rhm(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_nph(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_kls(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_oox(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_gh_amp(unsigned n, const double *x, double *grad, void *SS_ref_db);
int GH_sf_constrained(const char *name);
int GH_n_sf(const char *name, int n_em);
void GH_sf_of_p(const char *name, const double *p, double *sf);

void GH_SS_objective_init_function(    obj_type            *SS_objective,
                                        global_variable      gv                  );

void GH_PC_init(                       PC_type             *PC_read,
                                        global_variable      gv                  );

SS_ref GH_PC_function(     global_variable      gv,
                            PC_type             *PC_read,
                            SS_ref               SS_ref_db,
                            bulk_info            z_b,
                            int                  ph_id                  );

typedef void (*P2X_type) (				void 				*SS_ref_db,
										double 				 eps				);

void p2x_gh_generic(void *SS_ref_db, double eps);

void GH_P2X_init(                      P2X_type            *P2X_read,
                                        global_variable      gv                  );

#endif
