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
#ifndef __OBJECTIVE_FUNCTIONS_H_
#define __OBJECTIVE_FUNCTIONS_H_

#include "../MAGEMin.h"



void TC_mb_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable  	 gv					);

void TC_ig_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable   	 gv					);

void TC_igad_objective_init_function(	obj_type 		 	*SS_objective,
										global_variable   	 gv					);

void TC_mp_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable 	 gv					);

void TC_um_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable  	 gv					);

void TC_um_ext_objective_init_function(	obj_type 		 	*SS_objective,
									global_variable  	 gv					);
void TC_SS_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv					);

void TC_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv					);


typedef void (*P2X_type) (			void 				*SS_ref_db,
									double 				 eps				);

void TC_P2X_init(	                P2X_type 			*P2X_read,
									global_variable 	 gv					);
									
double obj_mb_liq(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_hb(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_aug(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_dio(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_opx(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_g(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_ol(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_fsp(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_abc(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_k4tr(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_sp(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_spn(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_ilm(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_ilmm(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_ep(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_bi(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_mu(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mb_chl(unsigned   n, const double *x, double *grad, void *SS_ref_db);

double obj_ig_fper(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_bi(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_cd(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_cpx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_ep(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_fl(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_g(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_hb(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_ilm(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_liq(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_mu(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_ol(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_opx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_fsp(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_spn(unsigned  n, const double *x, double *grad, void *SS_ref_db);

double obj_igad_liq(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_fsp(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_spn(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_g(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_ol(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_opx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_cpx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_ilm(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_ness(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_lct(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_kals(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igad_mel(unsigned  n, const double *x, double *grad, void *SS_ref_db);

double obj_mp_liq(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_fsp(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_bi(unsigned  	 n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_g(unsigned   	 n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_ep(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_ma(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_mu(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_opx(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_sa(unsigned  	 n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_cd(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_st(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_chl(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_ctd(unsigned 	 n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_sp(unsigned  	 n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_ilm(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_ilmm(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_mt(unsigned  	 n, const double *x, double *grad, void *SS_ref_db);
	
double obj_um_fluid(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_um_ol(unsigned     n, const double *x, double *grad, void *SS_ref_db);
double obj_um_br(unsigned     n, const double *x, double *grad, void *SS_ref_db);
double obj_um_ch(unsigned     n, const double *x, double *grad, void *SS_ref_db);
double obj_um_atg(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_um_g(unsigned      n, const double *x, double *grad, void *SS_ref_db);
double obj_um_ta(unsigned     n, const double *x, double *grad, void *SS_ref_db);
double obj_um_chl(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_um_anth(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_um_spi(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_um_opx(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_um_po(unsigned     n, const double *x, double *grad, void *SS_ref_db);
double obj_ume_pl4tr(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_ume_hb(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_ume_aug(unsigned   n, const double *x, double *grad, void *SS_ref_db);

double obj_aq17(unsigned 	  n, const double *x, double *grad, void *SS_ref_db);

SS_ref PC_function(			global_variable 	 gv,
							PC_type             *PC_read,
							
							SS_ref 				 SS_ref_db, 
							bulk_info 	 		 z_b,
							int    			 	ph_id					);
																		
void TC_mp_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv				);
void TC_mb_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv				);
void TC_ig_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv				);
void TC_igad_PC_init(	            PC_type 			*PC_read,
									global_variable 	 gv				);
void TC_um_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv				);
void TC_um_ext_PC_init(	            PC_type 			*PC_read,
									global_variable 	 gv				);
#endif
