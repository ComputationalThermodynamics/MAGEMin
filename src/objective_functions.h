#ifndef __OBJECTIVE_FUNCTIONS_H_
#define __OBJECTIVE_FUNCTIONS_H_

#include "MAGEMin.h"


void SS_ig_objective_init_function(	obj_type 		*SS_objective,
									global_variable  gv					);

void SS_igd_objective_init_function(obj_type 		*SS_objective,
									global_variable  gv					);

void SS_alk_objective_init_function(obj_type 		*SS_objective,
									global_variable  gv					);
					
void SS_mp_objective_init_function(	obj_type 		*SS_objective,
									global_variable  gv					);

void SS_um_objective_init_function(	obj_type 		*SS_objective,
									global_variable  gv					);


void p2x_alk_liq(SS_ref SS_ref_db, double eps);
void p2x_alk_fl(  SS_ref SS_ref_db, double eps);
void p2x_alk_fsp( SS_ref SS_ref_db, double eps);
void p2x_alk_spn( SS_ref SS_ref_db, double eps);
void p2x_alk_g(   SS_ref SS_ref_db, double eps);
void p2x_alk_ol(  SS_ref SS_ref_db, double eps);
void p2x_alk_opx( SS_ref SS_ref_db, double eps);
void p2x_alk_cpx( SS_ref SS_ref_db, double eps);
void p2x_alk_ilm( SS_ref SS_ref_db, double eps);
void p2x_alk_ness(SS_ref SS_ref_db, double eps);
void p2x_alk_lct( SS_ref SS_ref_db, double eps);
void p2x_alk_kals(SS_ref SS_ref_db, double eps);
void p2x_alk_mel( SS_ref SS_ref_db, double eps);
void p2x_alk_hb(  SS_ref SS_ref_db, double eps);
void p2x_alk_bi(  SS_ref SS_ref_db, double eps);
void p2x_alk_ep(  SS_ref SS_ref_db, double eps);
void p2x_alk_cd(  SS_ref SS_ref_db, double eps);

void p2x_igd_liq( SS_ref SS_ref_db, double eps);
void p2x_igd_fl(  SS_ref SS_ref_db, double eps);
void p2x_igd_fsp( SS_ref SS_ref_db, double eps);
void p2x_igd_spn( SS_ref SS_ref_db, double eps);
void p2x_igd_g(   SS_ref SS_ref_db, double eps);
void p2x_igd_ol(  SS_ref SS_ref_db, double eps);
void p2x_igd_opx( SS_ref SS_ref_db, double eps);
void p2x_igd_cpx( SS_ref SS_ref_db, double eps);
void p2x_igd_ilm( SS_ref SS_ref_db, double eps);
void p2x_igd_hb(  SS_ref SS_ref_db, double eps);
void p2x_igd_bi(  SS_ref SS_ref_db, double eps);
void p2x_igd_ep(  SS_ref SS_ref_db, double eps);
void p2x_igd_cd(  SS_ref SS_ref_db, double eps);

void p2x_ig_bi(  SS_ref SS_ref_db, double eps);
void p2x_ig_cd(  SS_ref SS_ref_db, double eps);
void p2x_ig_cpx( SS_ref SS_ref_db, double eps);
void p2x_ig_ep(  SS_ref SS_ref_db, double eps);
void p2x_ig_fl(  SS_ref SS_ref_db, double eps);
void p2x_ig_g(   SS_ref SS_ref_db, double eps);
void p2x_ig_hb(  SS_ref SS_ref_db, double eps);
void p2x_ig_ilm( SS_ref SS_ref_db, double eps);
void p2x_ig_liq( SS_ref SS_ref_db, double eps);
void p2x_ig_mu(  SS_ref SS_ref_db, double eps);
void p2x_ig_ol(  SS_ref SS_ref_db, double eps);
void p2x_ig_opx( SS_ref SS_ref_db, double eps);
void p2x_ig_pl4T(SS_ref SS_ref_db, double eps);
void p2x_ig_spn( SS_ref SS_ref_db, double eps);

void p2x_mp_liq(  	SS_ref SS_ref_db, double eps);
void p2x_mp_pl4tr(  SS_ref SS_ref_db, double eps);
void p2x_mp_bi( 	SS_ref SS_ref_db, double eps);
void p2x_mp_g(  	SS_ref SS_ref_db, double eps);
void p2x_mp_ep(  	SS_ref SS_ref_db, double eps);
void p2x_mp_ma(   	SS_ref SS_ref_db, double eps);
void p2x_mp_mu(  	SS_ref SS_ref_db, double eps);
void p2x_mp_opx( 	SS_ref SS_ref_db, double eps);
void p2x_mp_sa( 	SS_ref SS_ref_db, double eps);
void p2x_mp_cd(  	SS_ref SS_ref_db, double eps);
void p2x_mp_st(  	SS_ref SS_ref_db, double eps);
void p2x_mp_chl( 	SS_ref SS_ref_db, double eps);
void p2x_mp_ctd(	SS_ref SS_ref_db, double eps);
void p2x_mp_sp( 	SS_ref SS_ref_db, double eps);
void p2x_mp_ilm(	SS_ref SS_ref_db, double eps);
void p2x_mp_mt( 	SS_ref SS_ref_db, double eps);

void p2x_um_fluid(  SS_ref SS_ref_db, double eps);
void p2x_um_ol(  	SS_ref SS_ref_db, double eps);
void p2x_um_br( 	SS_ref SS_ref_db, double eps);
void p2x_um_ch(  	SS_ref SS_ref_db, double eps);
void p2x_um_atg( 	SS_ref SS_ref_db, double eps);
void p2x_um_g(   	SS_ref SS_ref_db, double eps);
void p2x_um_ta(  	SS_ref SS_ref_db, double eps);
void p2x_um_chl( 	SS_ref SS_ref_db, double eps);
void p2x_um_anth( 	SS_ref SS_ref_db, double eps);
void p2x_um_spi(  	SS_ref SS_ref_db, double eps);
void p2x_um_opx(  	SS_ref SS_ref_db, double eps);
void p2x_um_po( 	SS_ref SS_ref_db, double eps);

double obj_alk_liq(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_fl(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_fsp(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_spn(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_g(unsigned    n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_ol(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_opx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_cpx(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_ilm(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_ness(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_lct(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_kals(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_mel(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_hb(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_bi(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_ep(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_alk_cd(unsigned   n, const double *x, double *grad, void *SS_ref_db);


double obj_igd_liq(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_fl(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_spn(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_g(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_ol(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_opx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_hb(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_bi(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_ep(unsigned  n, const double *x, double *grad, void *SS_ref_db);
double obj_igd_cd(unsigned  n, const double *x, double *grad, void *SS_ref_db);

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
double obj_ig_pl4T(unsigned n, const double *x, double *grad, void *SS_ref_db);
double obj_ig_spn(unsigned  n, const double *x, double *grad, void *SS_ref_db);

double obj_mp_liq(unsigned   n, const double *x, double *grad, void *SS_ref_db);
double obj_mp_pl4tr(unsigned n, const double *x, double *grad, void *SS_ref_db);
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


SS_ref PC_function(			global_variable 	 gv,
							SS_ref 				 SS_ref_db, 
							bulk_info 	 		 z_b,
							char    			*name					);
													
SS_ref P2X(					global_variable 	 gv,
							SS_ref 				 SS_ref_db, 
							bulk_info 	 		 z_b,
							char    			*name					);	
							
int get_phase_id(			global_variable 	 gv,
							char    			*name					);
#endif
