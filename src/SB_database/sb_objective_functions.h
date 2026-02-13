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

void SB_sb21_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable  	 gv					);
									
void SB_sb24_objective_init_function(	obj_type 		 	*SS_objective,
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

//headers for objectives functions
double obj_sb21_plg(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_sp(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_ol(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_wa(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_ri(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_hpcpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_ak(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_gtmj(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_pv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_ppv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_cf(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_mw(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb21_nal(unsigned n, const double *x, double *grad, void *SS_ref_db);

//headers for objectives functions
double obj_sb24_plg(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_sp(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_ol(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_wa(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_ri(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_hpcpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_ak(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_gtmj(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_pv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_ppv(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_cf(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_mw(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_sb24_nal(unsigned n, const double *x, double *grad, void *SS_ref_db);

SS_ref SB_PC_function(		global_variable 	 gv,
							PC_type             *PC_read,
							
							SS_ref 				 SS_ref_db, 
							bulk_info 	 		 z_b,
							int    			 	ph_id					);
																		
void SB_sb11_PC_init(	               PC_type 			*PC_read,
									global_variable 	 gv				);

void SB_sb21_PC_init(	               PC_type 			*PC_read,
									global_variable 	 gv				);

void SB_sb24_PC_init(	               PC_type 			*PC_read,
									global_variable 	 gv				);

#endif
