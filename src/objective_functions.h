#ifndef __OBJECTIVE_FUNCTIONS_H_
#define __OBJECTIVE_FUNCTIONS_H_

#include "MAGEMin.h"


void SS_objective_init_function(	obj_type 		*SS_objective,
									global_variable  gv					);

void p2x_bi(  SS_ref SS_ref_db, double eps);
void p2x_cd(  SS_ref SS_ref_db, double eps);
void p2x_cpx( SS_ref SS_ref_db, double eps);
void p2x_ep(  SS_ref SS_ref_db, double eps);
void p2x_fl(  SS_ref SS_ref_db, double eps);
void p2x_g(   SS_ref SS_ref_db, double eps);
void p2x_hb(  SS_ref SS_ref_db, double eps);
void p2x_ilm( SS_ref SS_ref_db, double eps);
void p2x_liq( SS_ref SS_ref_db, double eps);
void p2x_mu(  SS_ref SS_ref_db, double eps);
void p2x_ol(  SS_ref SS_ref_db, double eps);
void p2x_opx( SS_ref SS_ref_db, double eps);
void p2x_pl4T(SS_ref SS_ref_db, double eps);
void p2x_spn( SS_ref SS_ref_db, double eps);

double obj_bi(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_cd(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_cpx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ep(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_fl(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_g(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_hb(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ilm(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_liq(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_mu(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ol(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_opx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_pl4T(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_spn(unsigned  n, const double *x, double *grad, void *SS_ref_db);

SS_ref PC_PX_function(		SS_ref 				 SS_ref_db, 
							double  			*x,
							char    			*name					);	
													
SS_ref PC_function(			global_variable 	 gv,
							SS_ref 				 SS_ref_db, 
							struct bulk_info 	 z_b,
							char    			*name					);
													
SS_ref P2X(					global_variable 	 gv,
							SS_ref 				 SS_ref_db, 
							struct bulk_info 	 z_b,
							char    			*name					);	
							
int get_phase_id(			global_variable 	 gv,
							char    			*name					);
#endif
