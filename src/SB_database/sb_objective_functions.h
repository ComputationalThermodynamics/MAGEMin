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
#ifndef __SB_OBJECTIVE_FUNCTIONS_H_
#define __SB_OBJECTIVE_FUNCTIONS_H_

#include "../MAGEMin.h"



void SB_sb11_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable  	 gv					);
void SB_SS_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv					);

void SB_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv					);

//headers for objectives functions
double obj_sb11_plg(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_sp(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_ol(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_wa(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_ri(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_hpcpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_ak(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_gtmj(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_pv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_ppv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_mw(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb11_cf(unsigned n, const double *x, double *grad, void *SS_ref_db);


SS_ref SB_PC__function(		global_variable 	 gv,
							PC_type             *PC_read,
							
							SS_ref 				 SS_ref_db, 
							bulk_info 	 		 z_b,
							int    			 	ph_id					);
																		
void SB_sb11_PC_init(	               PC_type 			*PC_read,
									global_variable 	 gv				);

#endif
