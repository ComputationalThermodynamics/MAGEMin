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
/**
               Partitioning Gibbs Energy routine            

The routine is the core of MAGEMin algorithm and is constructed around the Gibbs-Duhem constraint. It for a coupled system of 3 equations:

- mass constraint with phase fraction expressed as function of chemical potential of pure components          
- sum of endmember fractions within a solution phase must be equal to 1.0                                            
- delta_G of pure phase must lie on the G_hyperplane     
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#if __APPLE__
	//extern void dgetrf( int* M, int* N, double* A, int* lda, int* ipiv, int* info);
	//extern void dgetrs(char* T, int* N, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
	extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

#else
	#include <lapacke.h> 
#endif 

#include "MAGEMin.h"
#include "simplex_levelling.h"
#include "toolkit.h"
#include "gem_function.h"
#include "dump_function.h"
#include "phase_update_function.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "all_solution_phases.h"

#define FN_MAX_ITE 			8
#define FN_MAX_LS 			12
#define FN_MAX_ATT 			6
#define FN_MAX_DROP 		2
#define FN_TOL 				1e-9
#define FN_H 				1e-6
#define FN_H2 				1e-5
#define FN_S_FIX 			1e-9
#define FN_N_TOL 			1e-12
#define FN_X_TOL 			1e-9
#define FN_X_ACT 			1e-6
#define FN_PROJ_TOL 		1e-8
#define FN_FRAC 			0.5
#define FN_FRAC_S 			0.995
#define FN_ROOM 			0.01
#define FN_G_TOL 			1e-7
#define FN_G_MB 			1e6
#define FN_DF_TOL 			1e-8
#define FN_MU_TOL 			1e-6
#define FN_DEP_TOL 			1e-8
#define FN_SNAP_TOL 		1e-12

/** 
  Partitioning Gibbs Energy function 
*/
void PGE_print(					bulk_info 				z_b,
								global_variable  		gv,

								PP_ref 					*PP_ref_db,
								SS_ref 					*SS_ref_db,
								csd_phase_set  			*cp
){
	printf("\n _________________________________________________________________\n");
	printf("                          PHASE ASSEMBLAGE                        \n");
	printf(" ═════════════════════════════════════════════════════════════════\n\n");
	printf("ON | phase |  Fraction |  delta_G   |  factor   |   sum_xi   |  N(pi-xi) |    Pi - Xi...\n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1 ){

			printf(" %d | %4s | %+10f | %+10f | %+10f | %+10f | ",cp[i].ss_flags[1],cp[i].name,cp[i].ss_n,cp[i].df,cp[i].factor,cp[i].sum_xi);

			printf(" %+6f |", sum_norm_xipi(cp[i].xi_em,cp[i].p_em,cp[i].n_em) );

			for (int k = 0; k < cp[i].n_em; k++) {
				printf(" %+6f",(cp[i].p_em[k]-cp[i].xi_em[k]*cp[i].p_em[k])*SS_ref_db[cp[i].id].z_em[k]);
			}
			printf("\n");
						
		}
	}
	
	printf("\n");
	printf("ON | phase |  xeos\n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].ss_flags[1] == 1 ){

			printf(" %d | %4s |",cp[i].ss_flags[1],cp[i].name);

			for (int k = 0; k < cp[i].n_xeos; k++) {
				printf(" %+6f",cp[i].xeos[k]);
			}
			printf("\n");
						
		}
	}

	if (gv.n_pp_phase > 0){
		printf("\n");
		printf("ON | P. phase |  Fraction  |  delta_G   |  factor   | \n");
		printf("═══════════════════════════════════════════════════════════════════════════════════════\n");
		for (int i = 0; i < gv.len_pp; i++){ 
			if (gv.pp_flags[i][1] == 1){
				printf(" %d | %4s     | %+10f | %+10f | %+10f | \n",1,gv.PP_list[i],gv.pp_n[i],PP_ref_db[i].gb_lvl*PP_ref_db[i].factor,PP_ref_db[i].factor);
			}
		}	
	}

	printf("\n");
	printf("OFF| phase |  Fraction |  delta_G   |  factor   |   sum_xi   |  N(pi-xi) |  Pi - Xi...\n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].ss_flags[2] == 1){

			printf(" %d | %4s | %+10f | %+10f | %+10f | %+10f | ",cp[i].ss_flags[1],cp[i].name,cp[i].ss_n,cp[i].df*cp[i].factor,cp[i].factor,cp[i].sum_xi);

			printf(" %+6f |", sum_norm_xipi(cp[i].xi_em,cp[i].p_em,cp[i].n_em) );

			for (int k = 0; k < cp[i].n_em; k++) {
				printf(" %+6f",(cp[i].p_em[k]-cp[i].xi_em[k]*cp[i].p_em[k])*SS_ref_db[cp[i].id].z_em[k]);
			}
			printf("\n");
						
		}
	}
	printf("\n");
	printf("OFF| P. phase |  Fraction  |  delta_G  (< 5.0) | \n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");
	for (int i = 0; i < gv.len_pp; i++){ 
		if ((gv.pp_flags[i][2] == 1 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor < 5.0) ){
			printf(" %d | %4s     | %+10f | %+10f | \n",0,gv.PP_list[i],gv.pp_n[i],PP_ref_db[i].gb_lvl*PP_ref_db[i].factor);
		}
	}

	printf("\n\n ════════");
	for (int i = 0; i < z_b.nzEl_val; i++){		
		printf("════════════");
	}
	printf("\n");
	printf(" Oxide  |");
	for (int i = 0; i < z_b.nzEl_val; i++){		
		printf(" %11s",gv.ox[z_b.nzEl_array[i]]);
	}
	printf("\n"); 
	printf(" Gamma  |");
	for (int i = 0; i < z_b.nzEl_val; i++){	
		if	(gv.gam_tot[z_b.nzEl_array[i]] <= -1000.0){ 
			printf(" %.5f",gv.gam_tot[z_b.nzEl_array[i]]);
		}
		else{
			printf(" %.6f",gv.gam_tot[z_b.nzEl_array[i]]);
		}
	}
	printf("\n"); 
	printf(" dGamma |");
	for (int i = 0; i < z_b.nzEl_val; i++){	
		printf(" %+11f",gv.dGamma[z_b.nzEl_array[i]]);
	}
	printf("  -> *%.5f",gv.alpha);
	printf("\n\n");
	printf(" [GIBBS SYSTEM (Gibbs-Duhem) %.8f (with mu %.8f)]\n",gv.G_system,gv.G_system_mu);
	printf(" [MASS RESIDUAL NORM  = %+.8f ]\n",gv.BR_norm);
};

/** 
  Ipdate mass-constraint residual
*/
global_variable PGE_residual_update(			bulk_info 				z_b,
												global_variable  		gv,

												PP_ref 					*PP_ref_db,
												SS_ref 					*SS_ref_db,
												csd_phase_set  			*cp
){
	int ss;

	/**  if we are doing linear programming */
	if (gv.LP == 1 && gv.PGE == 0){
		
		for (int j = 0; j < gv.len_ox; j++){
		gv.mass_residual[j] = -z_b.bulk_rock[j];
			for (int i = 0; i < gv.len_pp; i++){
				if (gv.pp_flags[i][1] == 1){
					gv.mass_residual[j] += PP_ref_db[i].Comp[j]*PP_ref_db[i].factor*gv.pp_n[i];
				}
			}	

			/** calculate residual as function of endmember fractions */
			for (int i = 0; i < gv.len_cp; i++){
				if (cp[i].ss_flags[1] == 1 ){
					ss = cp[i].id;
					for (int k = 0; k < cp[i].n_em; k++){
						gv.mass_residual[j] += SS_ref_db[ss].Comp[k][j]*cp[i].factor*cp[i].p_em[k]*SS_ref_db[ss].z_em[k]*cp[i].ss_n;
					}
				}
			}
		}	
	}
	
	if(gv.LP == 0 && gv.PGE == 1){
		for (int j = 0; j < gv.len_ox; j++){
		gv.mass_residual[j] = -z_b.bulk_rock[j];
			for (int i = 0; i < gv.len_pp; i++){
				if (gv.pp_flags[i][1] == 1){
					gv.mass_residual[j] += PP_ref_db[i].Comp[j]*PP_ref_db[i].factor*gv.pp_n[i];
				}
			}	

			/** calculate residual as function xi fraction and not endmember fractions from x-eos */
			for (int i = 0; i < gv.len_cp; i++){
				if (cp[i].ss_flags[1] == 1 ){
					ss = cp[i].id;
					for (int k = 0; k < cp[i].n_em; k++){
						gv.mass_residual[j] += SS_ref_db[ss].Comp[k][j]*cp[i].factor*cp[i].p_em[k]*cp[i].xi_em[k]*SS_ref_db[ss].z_em[k]*cp[i].ss_n;
					}
				}
			}
		}
	}

	double br_norm = 0.0;
	for (int v = 0; v < z_b.nzEl_val; v++){
		int    j = z_b.nzEl_array[v];
		double r = gv.mass_residual[j];
		if (gv.BR_rel_norm && fabs(z_b.bulk_rock[j]) > 1e-8){ r /= fabs(z_b.bulk_rock[j]); }
		br_norm += r*r;
	}
	gv.BR_norm = sqrt(br_norm);

	/* Calculate G-system */
	gv.G_system = 0.0;
	for (int j = 0; j < gv.len_ox; j++){ gv.G_system += z_b.bulk_rock[j]*gv.gam_tot[j]; }
	
	gv.G_system_mu = gv.G_system;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			for (int j = 0; j < cp[i].n_em; j++){
				gv.G_system_mu +=  cp[i].ss_n*cp[i].p_em[j]*cp[i].mu[j]*cp[i].factor;
			}
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			gv.G_system_mu +=  gv.pp_n[i]*PP_ref_db[i].gb_lvl*PP_ref_db[i].factor;
		}
	}
	gv.gibbs_ev[gv.global_ite] = gv.G_system;

   return gv;
};

/** 
  Function to update chemical potential of endmembers (mui)
*/
global_variable PGE_update_mu(		bulk_info 	z_b,
									global_variable  	gv,

									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp
){
	int 	ss;
	double  cv;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1){
			ss  = cp[i].id;

			/** rotate gbase with respect to the G-hyperplane (change of base) */
			for (int k = 0; k < cp[i].n_em; k++) {
				cp[i].delta_mu[k] = 0.0;
				for (int j = 0; j < gv.len_ox; j++) {
					cp[i].delta_mu[k] 	-= SS_ref_db[ss].Comp[k][j]*gv.delta_gam_tot[j];
				}
				cp[i].mu[k] += cp[i].delta_mu[k];
				cp[i].df 	+= cp[i].delta_mu[k]*cp[i].p_em[k];
			}
		}
	}

   return gv;
};

/** 
  Partitioning Gibbs Energy function 
*/
global_variable PGE_update_xi(		bulk_info 	z_b,
									global_variable  	gv,

									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp
){
	int ss;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1){
			ss = cp[i].id;

			cp[i] = CP_UPDATE_function(			gv, 
												SS_ref_db[ss], 
												cp[i],
												z_b							);							
		}	
	}

   return gv;
};

/** 
	function to fill LHS (J)
*/
void PGE_build_Jacobian( 	double 			    *A,
							bulk_info 	 		 z_b,
							global_variable  	 gv,
							PP_ref 				*PP_ref_db,
							SS_ref 				*SS_ref_db,
							csd_phase_set  		*cp,
							int 				 nEntry
){
	int i,j,k,l,v,x,ph,ss,ix,ix0;
	
	/* 1) fill the Top Left corner of the matrix with: fv = sum(nl*sum(a_ij*a_iv*xi_l)) [nzEl_val * nzEl_val entries] */
	for (v = 0; v < z_b.nzEl_val; v++){
		/* CONSTRUCT LHS */
		for (j = 0; j < z_b.nzEl_val; j++){
			ix = v*nEntry + j;
			A[ix] = 0.;

			for (i = 0; i < gv.n_cp_phase; i++){
				ph = gv.cp_id[i];
				ss = cp[ph].id;
				for (x = 0; x < cp[ph].n_em; x++){
					/* CONSTRUCT TL CORNER */
					A[ix] 					+= 	SS_ref_db[ss].Comp[x][z_b.nzEl_array[j]] * cp[ph].factor * 
												SS_ref_db[ss].Comp[x][z_b.nzEl_array[v]] * cp[ph].factor * 
												cp[ph].xi_em[x]*cp[ph].p_em[x] * cp[ph].ss_n * SS_ref_db[ss].z_em[x];
				}
			}

		}
	}

	/* 2) fill the Middle Left part of the matrix with: hl = sum(a_ij*xi_l)) [n_ss_phase * nzEl_val entries] */
	for (l = 0; l < gv.n_cp_phase; l++){
		ph      = gv.cp_id[l];
		ss 		= cp[ph].id;
		
		/* CONSTRUCT LHS */
		for (j = 0; j < z_b.nzEl_val; j++){								/** shifts by +z_b.nzEl_val */
			ix = (l+z_b.nzEl_val)*nEntry + j;
			A[ix] =  0.0;
			for (i = 0; i < cp[ph].n_em; i++){
				A[ix] +=    SS_ref_db[ss].Comp[i][z_b.nzEl_array[j]] * cp[ph].factor * (cp[ph].p_em[i]*cp[ph].xi_em[i])  * SS_ref_db[ss].z_em[i];		
			}
		}
	}

	/* 3) fill the Bottom Left part of the matrix with: qk = a_ik [n_pp_phase * nzEl_val entries] */
	/* IS THERE A z_b.nzEl_array[v] */
	for (k = 0; k < gv.n_pp_phase; k++){
		ph    = gv.pp_id[k];

		for (v = 0; v < z_b.nzEl_val; v++){ 
			ix = (k+z_b.nzEl_val+gv.n_cp_phase)*nEntry + v;
			A[ix] = PP_ref_db[ph].Comp[z_b.nzEl_array[v]] * PP_ref_db[ph].factor;
		}
	}

	/* 4) fill the TM */
	for (l = 0; l < gv.n_cp_phase; l++){
		ph      = gv.cp_id[l];
		ss 		= cp[ph].id;
		
		/* CONSTRUCT LHS */
		for (j = 0; j < z_b.nzEl_val; j++){
			ix0    = j*nEntry + l + z_b.nzEl_val;
			A[ix0] =  0.0;
			for (i = 0; i < cp[ph].n_em; i++){
				A[ix0] +=  SS_ref_db[ss].Comp[i][z_b.nzEl_array[j]] * cp[ph].factor * (cp[ph].p_em[i]*cp[ph].xi_em[i])  * SS_ref_db[ss].z_em[i];		
			}
		}
	}

	/* 5) fill the TR corner by symmetry */
	for (i = z_b.nzEl_val+gv.n_cp_phase; i < nEntry; i++){
		for (j = 0; j < z_b.nzEl_val; j++){
			ix  = i*nEntry + j;
			ix0 = j*nEntry + i;
			
			A[ix0] = A[ix];
		}
	}
}

/** 
	function to fill RHS gradient
*/
void PGE_build_gradient( 	double				*b,
							bulk_info 	 		 z_b,
							global_variable  	 gv,
							PP_ref 				*PP_ref_db,
							SS_ref 				*SS_ref_db,
							csd_phase_set  		*cp,
							int 				 nEntry
){
	int i,j,k,l,v,x,ph,ss;
	double fac = -1.0;

	/* 1) fill the Top Left corner of the matrix with: fv = sum(nl*sum(a_ij*a_iv*xi_l)) [nzEl_val * nzEl_val entries] */
	for (v = 0; v < z_b.nzEl_val; v++){
		/* CONSTRUCT RHS */
		b[v]   = - z_b.bulk_rock[z_b.nzEl_array[v]];
		
		for (i = 0; i < gv.n_cp_phase; i++){
			ph = gv.cp_id[i];
			ss = cp[ph].id;
			for (x = 0; x < cp[ph].n_em; x++){
				b[v] 	+= SS_ref_db[ss].Comp[x][z_b.nzEl_array[v]] * cp[ph].factor * (cp[ph].p_em[x]*cp[ph].xi_em[x]) * cp[ph].ss_n * SS_ref_db[ss].z_em[x];
			}
		}
		
		for (k = 0; k < gv.n_pp_phase; k++){
			ph = gv.pp_id[k];
			b[v] 	    += PP_ref_db[ph].Comp[z_b.nzEl_array[v]] * PP_ref_db[ph].factor * gv.pp_n[ph] ;
		}
		b[v] *= fac;
	}

	/* 2) fill the Middle Left part of the matrix with: hl = sum(a_ij*xi_l)) [n_ss_phase * nzEl_val entries] */
	for (l = 0; l < gv.n_cp_phase; l++){
		ph      = gv.cp_id[l];
		ss 		= cp[ph].id;
		
		/* CONSTRUCT RHS */
		b[l+z_b.nzEl_val]    = -1.0;
		for (i = 0; i < cp[ph].n_em; i++){
			b[l+z_b.nzEl_val]    +=  (cp[ph].p_em[i]*cp[ph].xi_em[i])* SS_ref_db[ss].z_em[i];
		}
		b[l+z_b.nzEl_val] *= fac;
	}

	/* 3) fill the Bottom Left part of the matrix with: qk = a_ik [n_pp_phase * nzEl_val entries] */
	/* IS THERE A z_b.nzEl_array[v] */
	for (k = 0; k < gv.n_pp_phase; k++){
		ph    = gv.pp_id[k];

		b[k+z_b.nzEl_val+gv.n_cp_phase]        = -PP_ref_db[ph].gbase;	
		for (v = 0; v < z_b.nzEl_val; v++){ 
			b[k+z_b.nzEl_val+gv.n_cp_phase]   += PP_ref_db[ph].Comp[z_b.nzEl_array[v]] * gv.gam_tot[z_b.nzEl_array[v]];
		}
		b[k+z_b.nzEl_val+gv.n_cp_phase]  *= fac;
	}
}

/** 
  Partitioning Gibbs Energy function 
*/
global_variable PGE_update_solution(	global_variable  	 gv,
										bulk_info 	 		 z_b,
										csd_phase_set  		*cp	){
	int 	i, j, k, ph;										
	double 	n_fac, 
			g_fac, 
			alpha, 
			max_dG_ss,
			max_dG,
			max_dn,
			max_dnss,
			max_dnpp;

	/**
		calculate under relaxing factor
	*/
	for (i = 0; i < z_b.nzEl_val; i++){
		gv.dGamma[i] = gv.b_PGE[i];
	}
	for (i = 0; i < gv.n_cp_phase; i++){
		gv.dn_cp[i] = gv.b_PGE[i+z_b.nzEl_val];
	}
	for (i = 0; i < gv.n_pp_phase; i++){
		gv.dn_pp[i] = gv.b_PGE[i+z_b.nzEl_val + gv.n_cp_phase];
	}
	
	max_dG 		= norm_vector(gv.dGamma,z_b.nzEl_val);
	max_dnss 	= norm_vector(gv.dn_cp,gv.n_cp_phase);
	max_dnpp 	= norm_vector(gv.dn_pp,gv.n_pp_phase);
	max_dn 		= ((max_dnss < max_dnpp) ? (max_dnpp) : (max_dnss) );
	max_dG_ss   = gv.relax_PGE_val*exp(-8.0*pow(gv.BR_norm,0.28))+1.0;

	g_fac       = (gv.max_g_phase/max_dG_ss)/max_dG;
	n_fac   	= (gv.max_n_phase/max_dG_ss)/max_dn;
	alpha 		= ((n_fac < g_fac) 	 	? 	(n_fac) 		: (g_fac)	);
	alpha 		= ((alpha > gv.max_fac) ? 	(gv.max_fac) 	: (alpha)	);

	gv.alpha	= alpha;

	/* Update Gamma */
	for (i = 0; i < z_b.nzEl_val; i++){
		gv.delta_gam_tot[z_b.nzEl_array[i]]  = gv.dGamma[i]*gv.alpha;
		gv.gam_tot[z_b.nzEl_array[i]] 		+= gv.dGamma[i]*gv.alpha;
	}

	gv.gamma_norm[gv.global_ite] = norm_vector(gv.dGamma, z_b.nzEl_val);

	/* Update solution phase (SS) fractions */
	for (i = 0; i < gv.n_cp_phase; i++){
		 cp[gv.cp_id[i]].delta_ss_n  = gv.dn_cp[i]*gv.alpha;
		 cp[gv.cp_id[i]].ss_n 		+= gv.dn_cp[i]*gv.alpha;
	}
	
	/* Update pure phase (PP) fractions */
	if (gv.n_pp_phase > 0){
		for (i = 0; i < gv.n_pp_phase; i++){
			 gv.pp_n[gv.pp_id[i]] 		+= gv.dn_pp[i]*gv.alpha;
			 gv.delta_pp_n[gv.pp_id[i]]  = gv.dn_pp[i]*gv.alpha;
		}
	}
	
	return gv;
};

/** 
  Partitioning Gibbs Energy function 
*/
global_variable PGE_solver(		bulk_info 	 		 z_b,
								global_variable  	 gv,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp
){
	/* allocate */
	int 	i,j,k,l,v,x,ph,ss;

	/* extract the number of entries in the matrix */
	int 	nEntry = z_b.nzEl_val + gv.n_phase;
	
	/* LAPACKE memory allocation */
	int 	nrhs   = 1;													/** number of rhs columns, 1 is vector*/
	int 	lda    = nEntry;											/** leading dimesion of A*/
	int 	ldb    = 1;													/** leading dimension of b*/
	int 	info;														/** get info from lapacke function*/

	for (i = 0; i < z_b.nzEl_val;  i++){ gv.dGamma[i] 	= 0.0;	}		/** initialize dGamma to 0.0 */
	for (i = 0; i < gv.n_cp_phase; i++){ gv.dn_cp[i]  	= 0.0;	}		/** initialize dGamma to 0.0 */
	for (i = 0; i < gv.n_pp_phase; i++){ gv.dn_pp[i]  	= 0.0;	}		/** initialize dGamma to 0.0 */
    for (i = 0; i < nEntry*nEntry; i++){ gv.A_PGE[i]  	= 0.0;	}
	for (i = 0; i < nEntry; i++){ 		 gv.b_PGE[i]  	= 0.0;	}

	/**
		get id of active pure phases
	*/
	gv = get_pp_id(				gv					);
	
	/**
		get id of active solution phases
	*/
	gv = get_ss_id(				gv,
								cp					);
	
	/** 
		function to fill Jacobian
	*/
	PGE_build_Jacobian( 		gv.A_PGE,
								z_b,
								gv,

								PP_ref_db,
								SS_ref_db,
								cp,
								nEntry				);

	/** 
		function to fill gradient
	*/
	PGE_build_gradient( 		gv.b_PGE,
								z_b,
								gv,

								PP_ref_db,
								SS_ref_db,
								cp,
								nEntry				);

	// if (gv.verbose == 1){
	// 	int ix;
	// 	int v,j;
	// 	for (v = 0; v < nEntry; v++){
	// 		/* CONSTRUCT LHS */
	// 		double sum = 0.0;
	// 		for (j = 0; j < nEntry; j++){
	// 			ix = v*nEntry + j;
	// 			printf("%.3f ",gv.A_PGE[ix]);;
	// 		}
	// 		printf(" | %.3f\n",gv.b_PGE[v]);	
	// 	}
	// }

	/**
		save RHS vector 
	*/
	gv.fc_norm_t1 = norm_vector(	gv.b_PGE,
									nEntry		);
		

	/**
		call lapacke to solve system of linear equation using LU 
	*/
	#if __APPLE__
		/*
		// Factorisation
		dgetrf(&nEntry, &nEntry, gv.A_PGE, &nEntry, gv.ipiv, &info);

		// Solution (with transpose!)
		char T = 'T';
		dgetrs(						&T,
									&nEntry, 
									&nrhs, 
									gv.A_PGE,
									&nEntry, 
									gv.ipiv, 
									gv.b_PGE, 
									&nEntry,
									&info	);
		*/			

		// remark: apple accelerate uses column-major ordering, whereas lapacke below uses row-major. 
		// As long as the matrix is strictly symmetric this is fine; if not we have to reorder gv.A_PGE				
		dgesv(						&nEntry, 
									&nrhs, 
									gv.A_PGE,
									&nEntry, 
									gv.ipiv, 
									gv.b_PGE, 
									&nEntry,
									&info	);
		

	#else
		info = LAPACKE_dgesv(		LAPACK_ROW_MAJOR, 
									nEntry, 
									nrhs, 
									gv.A_PGE, 
									lda, 
									gv.ipiv,
									gv.b_PGE,
									ldb					);
	#endif
	if (info != 0){
		fprintf(stderr, "MAGEMin: PGE linear solve failed (dgesv info = %d, system size = %d); terminating this point as non-converged\n", info, nEntry);
		gv.div    = 1;
		gv.status = 4;
		return gv;
	}
	/**
		get solution and max values for the set of variables
	*/
	gv = PGE_update_solution(	gv,
								z_b,
								cp					);

	return gv;
};

/** 
  Partitioning Gibbs Energy function 
*/
global_variable PGE_inner_loop(		bulk_info 			 z_b,
									simplex_data	    *splx_data,
									global_variable  	 gv,

									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp
){
	clock_t u; 
	int 	PGEi   			= 0;
	double 	fc_norm_t0 		= 0.0;
	double 	delta_fc_norm 	= 1.0;

	/* transform to while if delta_phase fraction < val */
	while (PGEi < gv.inner_PGE_ite && delta_fc_norm > 1e-10 && gv.div == 0){
		u = clock();

		gv =	PGE_solver(					z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp							); 
				
								
		delta_fc_norm 	= fabs(gv.fc_norm_t1 - fc_norm_t0);
		fc_norm_t0 		= gv.fc_norm_t1;
							
		/**
			calculate delta_G of pure phases 
		*/
		pp_min_function(					gv,
											z_b,
											PP_ref_db				);
										
										
							
		/* Update mu of solution phase  */
		gv =	PGE_update_mu(				z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp						); 

		gv =	PGE_update_xi(				z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp						);  

		gv = 	phase_update_function(		z_b,								/** bulk rock constraint 				*/
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/
											SS_ref_db,							/** solution phase database 			*/ 
											cp						); 

		/** 
			Update mass constraint residual
		*/
		gv = PGE_residual_update(			z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp						);  
		
		u = clock() - u; 
		gv.inner_PGE_ite_time =(((double)u)/CLOCKS_PER_SEC*1000);
		PGEi += 1;
	} 
		
   return gv;
};

global_variable compute_xi_SD(				global_variable  		 gv,
											csd_phase_set  			*cp				){
	gv.mean_sum_xi 	= 0.0;
	gv.sigma_sum_xi = 0.0;
	for (int iss = 0; iss < gv.len_cp; iss++){ 
		if (cp[iss].ss_flags[0] == 1){
			gv.mean_sum_xi += cp[iss].sum_xi/gv.n_cp_phase;
		}
	}
	for (int iss = 0; iss < gv.len_cp; iss++){ 
		if (cp[iss].ss_flags[0] == 1){
			gv.sigma_sum_xi += pow(cp[iss].sum_xi-gv.mean_sum_xi,2.0);
		}
	}
	gv.sigma_sum_xi = sqrt(gv.sigma_sum_xi/gv.mean_sum_xi);
	if (gv.verbose ==1){
		printf("\n mean sum_xi: %+10f [sd: %+10f]\n",gv.mean_sum_xi,gv.sigma_sum_xi);
	}
	
	return gv;
}


/**
  function to run simplex linear programming during PGE with pseudocompounds 
*/	
global_variable run_LP(								bulk_info 			 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db
){

	if (gv.verbose == 1){
		printf("\n");
		printf("Linear-Programming stage [PGE pseudocompounds]\n");	
		printf("══════════════════════════════════════════════\n");	
	}

	simplex_data *d  = (simplex_data *) splx_data;

	int  k 		= 0;
	d->swp 		= 1;
	d->n_swp 	= 0;
	while (d->swp == 1 && k < 32){					/** as long as a phase can be added to the guessed assemblage, go on */
		k 		  += 1;
		d->swp     = 0;
		

		swap_PGE_pseudocompounds(			z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);	
		
		swap_pure_phases(					z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);											
	}
	if (gv.verbose == 1){
		printf("\n");
		printf("  -> number of swap loops: %d\n",k);	
	}

	/* update gamma of SS */
	update_local_gamma(						d->A1,
											d->g0_A,
											d->gamma_ss,
											d->n_Ox			);

	/* update global variable gamma */
	update_global_gamma_LU(					z_b,
											splx_data		);	

	/* copy gamma total to the global variables */
	for (int i = 0; i < gv.len_ox; i++){
		gv.dGamma[i]  = d->gamma_tot[i] - gv.gam_tot[i];
		gv.gam_tot[i] = d->gamma_tot[i];
	}
	gv.gamma_norm[gv.global_ite] = norm_vector(gv.dGamma, z_b.nzEl_val);

	if (gv.verbose == 1){
		printf("\n Total number of LP iterations: %d\n",k);	
		printf(" [----------------------------------------]\n");
		printf(" [  Ph  |   Ph PROP  |   g0_Ph    |  ix   ]\n");
		printf(" [----------------------------------------]\n");

		for (int i = 0; i < d->n_Ox; i++){
			if (d->ph_id_A[i][0] == 0){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]", "F.OX", d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 1){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]", gv.PP_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 2){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]\n", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
			}
			if (d->ph_id_A[i][0] == 3){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
				if (d->stage[i] == 1){
					for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
						printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_Ppc[d->ph_id_A[i][3]][ii] );
					}
				}
				else{
					for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
						printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_pc[d->ph_id_A[i][3]][ii] );
					}
				}
				printf("\n");
			}
		}
		printf(" [----------------------------------------]\n");
		printf(" [  OXIDE      GAMMA                      ]\n");
		printf(" [----------------------------------------]\n");
		for (int i = 0; i < d->n_Ox; i++){
			printf(" [ %5s %+15f                  ]\n", gv.ox[z_b.nzEl_array[i]], d->gamma_tot[z_b.nzEl_array[i]]);
		}
		printf(" [----------------------------------------]\n");
		printf(" [             %4d swaps                 ]\n", d->n_swp);
		printf(" [----------------------------------------]\n");
		
	}

	return gv;
}



/**
  function to run simplex linear programming during PGE with pseudocompounds 
*/	
global_variable run_LP_ig(							bulk_info 			 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db
){

	if (gv.verbose == 1){
		printf("\n");
		printf("Linear-Programming initial guess computation\n");	
		printf("══════════════════════════════════════════════\n");	
	}

	simplex_data *d  = (simplex_data *) splx_data;

	int  k 		= 0;
	d->swp 		= 1;
	d->n_swp 	= 0;
	while (d->swp == 1 && k < 9){					/** as long as a phase can be added to the guessed assemblage, go on */
		k 		  += 1;
		d->swp     = 0;
		
		swap_PGE_pseudocompounds(			z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);	
		
		swap_pure_phases(					z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);											
	}
	if (gv.verbose == 1){
		printf("\n");
		printf("  -> number of swap loops: %d\n",k);	
	}

	/* update gamma of SS */
	update_local_gamma(						d->A1,
											d->g0_A,
											d->gamma_ss,
											d->n_Ox			);

	/* update global variable gamma */
	update_global_gamma_LU(					z_b,
											splx_data		);	

	if (gv.verbose == 1){
		printf("\n Total number of LP_ig iterations: %d\n",k);	
		printf(" [----------------------------------------]\n");
		printf(" [  Ph  |   Ph PROP  |   g0_Ph    |  ix   ]\n");
		printf(" [----------------------------------------]\n");

		for (int i = 0; i < d->n_Ox; i++){
			if (d->ph_id_A[i][0] == 0){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]", "F.OX", d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 1){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]", gv.PP_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 2){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]\n", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
			}
			if (d->ph_id_A[i][0] == 3){
				printf(" ['%5s' %+10f  %+12.4f  %2d %2d ]", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0], d->stage[i]);
				if (d->stage[i] == 1){
					for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
						printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_Ppc[d->ph_id_A[i][3]][ii] );
					}
				}
				else{
					for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
						printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_pc[d->ph_id_A[i][3]][ii] );
					}
				}
				printf("\n");
			}
		}
		printf(" [----------------------------------------]\n");
		printf(" [  OXIDE      GAMMA IG                   ]\n");
		printf(" [----------------------------------------]\n");
		for (int i = 0; i < d->n_Ox; i++){
			printf(" [ %5s %+15f                  ]\n", gv.ox[z_b.nzEl_array[i]], d->gamma_tot[z_b.nzEl_array[i]]);
		}
		printf(" [----------------------------------------]\n");
		printf(" [             %4d swaps ig              ]\n", d->n_swp);
		printf(" [----------------------------------------]\n");
		
	}

	return gv;
}

/**
  function to run simplex linear programming during PGE with pseudocompounds 
*/	
global_variable init_LP(							bulk_info 	 		 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													PC_type				*PC_read,
													P2X_type			*P2X_read,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db,
													csd_phase_set  		*cp	
){
	simplex_data *d  = (simplex_data *) splx_data;

	double distance;
	double min_distance;
	double mid_dG;
	
	int ph_id, npc, id, id_min_distance;
	int id_cp = 0;
	int pc_id;
	int em_id;
	int add_phase;
	int i, j, k, ii;
	int m_pc;
	int n = gv.max_ss_size_cp;

	/**
	   reset variables
	*/
	for (i = 0; i < gv.len_pp; i++){
		gv.pp_flags[i][1]   = 0;
	}

	/* reset pure phases fractions and xi */
	for (i = 0; i < gv.len_pp; i++){		
		gv.pp_n[i] 		  = 0.0;
		gv.pp_n_mol[i]	  = 0.0;
		gv.pp_n_wt[i]	  = 0.0;
		gv.pp_n_vol[i]	  = 0.0;
		gv.delta_pp_n[i]  = 0.0;
		gv.pp_xi[i] 	  = 0.0;
		gv.delta_pp_xi[i] = 0.0;
	}

	gv.len_cp 		  	  = 0;
	gv.ph_change  	      = 0;
	gv.n_cp_phase         = 0;					/** reset the number of ss phases to start with */
	gv.n_pp_phase         = 0;					/** reset the number of pp phases to start with */
	gv.n_phase            = 0;

	/* reset solvi */
    for (i = 0; i < gv.len_ss; i++){	
        gv.n_solvi[i] = 0;
    }
	
	// for (int i = 0; i < gv.max_n_cp; i++){		
	for (int i = 0; i < gv.len_ox; i++){		
		strcpy(cp[i].name,"");						/* get phase name */	
		cp[i].in_iter			=  0;
		cp[i].split				=  0;
		cp[i].id 				= -1;				/* get phaseid */
		cp[i].n_xeos			=  0;				/* get number of compositional variables */
		cp[i].n_em				=  0;				/* get number of endmembers */
		cp[i].n_sf				=  0;			
		cp[i].df 				=  0.0;
		cp[i].factor 			=  0.0;
		
		for (int ii = 0; ii < gv.n_flags; ii++){
			cp[i].ss_flags[ii] 	= 0;
		}

		cp[i].ss_n        		= 0.0;				/* get initial phase fraction */
		cp[i].ss_n_mol      	= 0.0;				/* get initial phase mol fraction */
		cp[i].ss_n_wt			= 0.0;
		cp[i].ss_n_vol   		= 0.0;				/* get initial phase mol fraction */
		cp[i].delta_ss_n    	= 0.0;				/* get initial phase fraction */
		
		for (int ii = 0; ii < n; ii++){
			cp[i].p_em[ii]      = 0.0;
			cp[i].xi_em[ii]     = 0.0;
			cp[i].dguess[ii]    = 0.0;
			cp[i].xeos[ii]      = 0.0;
			cp[i].delta_mu[ii]  = 0.0;
			cp[i].dfx[ii]       = 0.0;
			cp[i].mu[ii]        = 0.0;
			cp[i].gbase[ii]     = 0.0;
			cp[i].ss_comp[ii]   = 0.0;
		}
		 
		for (int ii = 0; ii < n*2; ii++){
			cp[i].sf[ii]    	= 0.0;
		}
		cp[i].mass 				= 0.0;
		cp[i].volume 			= 0.0;
		cp[i].phase_density 	= 0.0;
		cp[i].phase_cp 			= 0.0;
	}

	/** 
		get initial conditions for active phases
	*/
	for (i = 0; i < d->n_Ox; i++){
		add_phase 	= 0;
		ph_id 		= d->ph_id_A[i][1];
			
		
		if (d->ph_id_A[i][0] == 0){						/* if phase if a fake oxide, do nothing! */

		}
		else if (d->ph_id_A[i][0] == 1 ){ 				/* if phase is a pure species */
			gv.pp_flags[ph_id][1] 	= 1;
			gv.pp_flags[ph_id][2] 	= 0;
			gv.pp_n[ph_id]          = d->n_vec[i];
			gv.n_pp_phase		   += 1;
			gv.n_phase 			   += 1;
		}
		else {											/* pure endmembers as solution phase */
			
			if (d->ph_id_A[i][0] == 2){
			em_id 					= d->ph_id_A[i][3];

			for (j = 0; j < SS_ref_db[ph_id].n_em; j++) {	
				SS_ref_db[ph_id].p[j] = gv.em2ss_shift;
			}
			SS_ref_db[ph_id].p[em_id] = 1.0 - gv.em2ss_shift*SS_ref_db[ph_id].n_em;
			
			(*P2X_read[ph_id])(		&SS_ref_db[ph_id],
									gv.bnd_val					);		
			}
		
			/* solution phase */
			if (d->ph_id_A[i][0] == 3 && d->stage[i] == 1){
				pc_id 					= d->ph_id_A[i][3];

				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					SS_ref_db[ph_id].iguess[ii]  = SS_ref_db[ph_id].xeos_Ppc[pc_id][ii];
				}
			}
			if (d->ph_id_A[i][0] == 3 && d->stage[i] == 0){
				pc_id 					= d->ph_id_A[i][3];

				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					SS_ref_db[ph_id].iguess[ii]  = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
				}
			}

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			SS_ref_db[ph_id] = PC_function(			gv,
													PC_read,
													SS_ref_db[ph_id], 
													z_b,
													ph_id						);

			SS_ref_db[ph_id] = SS_UPDATE_function(	gv, 
													SS_ref_db[ph_id], 
													z_b, 
													gv.SS_list[ph_id]			);

			strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */

			cp[id_cp].split 		= 0;							
			cp[id_cp].id 			= ph_id;						/* get phase id */
			cp[id_cp].n_xeos		= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
			cp[id_cp].n_em			= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
			cp[id_cp].n_sf			= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
			
			cp[id_cp].df			= 0.0;
			cp[id_cp].factor		= SS_ref_db[ph_id].factor;	
			
			cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
			cp[id_cp].ss_flags[1] 	= 1;
			cp[id_cp].ss_flags[2] 	= 0;
			cp[id_cp].sum_xi		= SS_ref_db[ph_id].sum_xi;		
			cp[id_cp].ss_n          = d->n_vec[i];			/* get initial phase fraction */

			for (ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
				cp[id_cp].p_em[ii]  = SS_ref_db[ph_id].p[ii];
				cp[id_cp].xi_em[ii]	= SS_ref_db[ph_id].xi_em[ii];
				cp[id_cp].mu[ii]	= SS_ref_db[ph_id].mu[ii];
			}

			for (ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
				cp[id_cp].dguess[ii]  = SS_ref_db[ph_id].iguess[ii];
				cp[id_cp].xeos[ii]    = SS_ref_db[ph_id].iguess[ii];
			}

			for (int ii = 0; ii < gv.len_ox; ii++){
				cp[id_cp].ss_comp[ii]	= SS_ref_db[ph_id].ss_comp[ii];
			}
			
			for (int ii = 0; ii < SS_ref_db[ph_id].n_sf; ii++){
				cp[id_cp].sf[ii]		= SS_ref_db[ph_id].sf[ii];
			}
			gv.n_solvi[ph_id] 	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;
		}
	}

	/* reinitialize the number of SS instances */
	for (i = 0; i < gv.len_ss; i++){
		gv.n_solvi[i] = 0;
	}

	/* get number of duplicated phases and their cp id */
	for (i = 0; i < gv.len_cp; i++){
		
		if (cp[i].ss_flags[0] == 1 ){
			ph_id = cp[i].id;
			SS_ref_db[ph_id].solvus_id[gv.n_solvi[ph_id]] = i;
			gv.n_solvi[ph_id] += 1;
		}
	}

	return gv;
}


/**
  function to run simplex linear programming during PGE with pseudocompounds 
*/	
global_variable update_cp_after_LP(					bulk_info 	 		 z_b,
													global_variable 	 gv,
													PC_type 			*PC_read,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db,
													csd_phase_set  		*cp	
){
	int 	ph_id;
	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[1] == 1){

			ph_id = cp[i].id;

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(		gv, 
														SS_ref_db[ph_id]		);

			/**
				establish a set of conditions to update initial guess for next round of local minimization 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos[k];
			}
			
			SS_ref_db[ph_id] = PC_function(				gv,
														PC_read,
														SS_ref_db[ph_id], 
														z_b,
														ph_id				);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
														SS_ref_db[ph_id],
														z_b,
														gv.SS_list[ph_id]		);

			/* DEW-only: unlike every other phase, xeos alone does not determine a
			   valid DEW state - it must also satisfy the internal charge-balance
			   condition Sum(m_i*z_i)=0, which only DEW_aq_min_multistart/_warmstart
			   (via NLopt_opt_DEW_function) actually enforces. This function copies
			   whatever xeos cp[i] already holds (ss_flags[1]==1 phases can reach here
			   without ever having been through ss_min_PGE/ss_min_LP's NLopt_opt call in
			   the same run, if they never graduate to ss_flags[0]==1), re-scored only
			   by the cheap PC_function G-evaluation used above - which, for DEW,
			   just reads whatever composition it is given and cannot detect or correct
			   a charge imbalance. Confirmed via a user-reported reproduction (mbe,
			   11.26kbar/581.25C, rm_list=[1,4,6,15,19,-13]) where a DEW instance
			   with chargeResidual=9.25 (raw, never-refined seed-derived: Ca+2/K+/Na+
			   elevated with no compensating anion) still reached the final output even
			   after NLopt_opt_DEW_function itself was hardened to never return an
			   unconverged candidate - because this function's copy_to_cp call never
			   goes through NLopt_opt_DEW_function at all. mole-fraction-weighted
			   Sum(x_i*z_i) is proportional to the true molality-weighted charge residual
			   (same positive Omega/x_water scale factor multiplies every species), so it
			   is equally valid as a charge-balance check and needs no extra molality
			   recomputation. 1e-3 safely separates a genuinely bad, unrefined candidate
			   (this case: ~0.2) from the worst-case innocent trace-level asymmetry of a
			   never-elevated baseline seed (~1e-4 at most, ~100 suppressed-or-not
			   species at 1e-6 each). */
			if ((strcmp(gv.SS_list[ph_id], "DEW") == 0 || strcmp(gv.SS_list[ph_id], "DEW_S24") == 0) && SS_ref_db[ph_id].sf_ok == 1){
				/* SS_ref_db[ph_id] is a shared per-PHASE-TYPE scratch struct, not
				   per-cp[]-instance - .xeos in particular is only ever written by
				   NLopt_opt_DEW_function (S.x), never by PC_function/obj_DEW
				   (which only write .iguess/.sf/.mu), so it still holds whatever the
				   LAST NLopt_opt call for this ph_id left there regardless of which
				   cp[i] this loop iteration is on - reading it here would silently
				   check the wrong instance whenever DEW has multiple simultaneous
				   cp[] entries (a real solvus). .iguess is the correct field: it was
				   just set from THIS cp[i].xeos a few lines above, and is what
				   PC_function actually read and what copy_to_cp will copy back. */
				double z_res = 0.0;
				int n_sp = SS_ref_db[ph_id].n_em - 1;   /* water excluded: z=0, see G_SS_DEW_function */
				for (int k = 0; k < n_sp; k++){
					z_res += SS_ref_db[ph_id].iguess[k]*SS_ref_db[ph_id].mat_phi[k];
				}
				if (fabs(z_res) > 1e-3){
					SS_ref_db[ph_id].sf_ok = 0;
				}
			}

			/**
				print solution phase informations (print has to occur before saving PC)
			*/
			if (gv.verbose == 1){
				print_SS_informations(  				gv,
														SS_ref_db[ph_id],
														ph_id					);
			}

			/* if site fractions are respected then save the minimized point */
			if (SS_ref_db[ph_id].sf_ok == 1){
				/**
					copy the minimized phase informations to cp structure
				*/
				copy_to_cp(								i,
														ph_id,
														gv,
														SS_ref_db,
														cp						);
			}
			else{
				if (gv.verbose == 1){
					printf(" !> SF [:%d] not respected for %4s (SS not updated)\n",SS_ref_db[ph_id].sf_id,gv.SS_list[ph_id]);
				}
			}
		}
	}

	return gv;
}

static void fn_bounds(				SS_ref 				*d,
									double 				 widen			){
	for (int j = 0; j < d->n_xeos; j++){
		d->bounds[j][0] = d->bounds_ref[j][0] - widen;
		d->bounds[j][1] = d->bounds_ref[j][1] + widen;
	}
}

static double fn_eval(				global_variable 	 gv,
									SS_ref 				*d,
									const double 		*gam,
									const double 		*s,
									double 				*g,
									double 				*c				){
	for (int k = 0; k < d->n_em; k++){
		d->gb_lvl[k] = d->gbase[k];
		for (int j = 0; j < gv.len_ox; j++){ d->gb_lvl[k] -= d->Comp[k][j]*gam[j]; }
	}
	double D = ns_eval(gv, d, s, g);
	if (c != NULL){
		for (int j = 0; j < gv.len_ox; j++){
			double a = 0.0;
			for (int k = 0; k < d->n_em; k++){ a += d->Comp[k][j]*d->p[k]*d->z_em[k]; }
			c[j] = d->factor*a;
		}
	}
	return D;
}

static double fn_step(				SS_ref 				*d,
									const double 		*s,
									const double 		*Nq,
									int 				 k,
									int 				 kk,
									double 				 h0				){
	double h = h0;
	for (int i = 0; i < d->ns_nc; i++){
		double v = fabs(Nq[i*k + kk]);
		if (v == 0.0){ continue; }
		double room = (s[i] < d->ns_ub[i] - s[i]) ? s[i] : d->ns_ub[i] - s[i];
		if (FN_ROOM*room/v < h){ h = FN_ROOM*room/v; }
	}
	return h;
}

static double fn_norm(				const double 		*R,
									int 				 nz				){
	double r = 0.0;
	for (int i = 0; i < nz; i++){
		double a = fabs(R[i]);
		if (!(a <= r)){ r = a; }
	}
	return r;
}

static int fn_orth_add(				double 				*Q,
									int 				 n_q,
									int 				 n,
									double 				*v				){
	double nv0 = 0.0;
	for (int i = 0; i < n; i++){ nv0 += v[i]*v[i]; }
	nv0 = sqrt(nv0);
	if (!(nv0 > 0.0)){ return 0; }
	for (int pass = 0; pass < 2; pass++){
		for (int r = 0; r < n_q; r++){
			double a = 0.0;
			for (int i = 0; i < n; i++){ a += Q[r*n + i]*v[i]; }
			for (int i = 0; i < n; i++){ v[i] -= a*Q[r*n + i]; }
		}
	}
	double nv = 0.0;
	for (int i = 0; i < n; i++){ nv += v[i]*v[i]; }
	nv = sqrt(nv);
	if (!(nv > FN_DEP_TOL*nv0)){ return 0; }
	for (int i = 0; i < n; i++){ Q[n_q*n + i] = v[i]/nv; }
	return 1;
}

static void fn_shift(				SS_ref 				*d,
									const double 		*s,
									const double 		*Nq,
									int 				 k,
									const double 		*dy,
									double 				*out			){
	for (int i = 0; i < d->ns_nc; i++){
		double a = 0.0;
		for (int kk = 0; kk < k; kk++){ a += Nq[i*k + kk]*dy[kk]; }
		out[i] = s[i] + a;
	}
}

static double fn_xval(				global_variable 	 gv,
									SS_ref 				*d,
									const double 		*s,
									int 				 j				){
	ns_x_of_sf(gv, d, s);
	return d->ns_x[j];
}

static void fn_xgrad(				global_variable 	 gv,
									SS_ref 				*d,
									const double 		*s,
									const double 		*Nq,
									int 				 k,
									int 				 j,
									double 				*g,
									double 				*H				){
	int 	nc = d->ns_nc;
	double 	sp[nc], dy[k > 0 ? k : 1];

	for (int kk = 0; kk < k; kk++){
		double h = FN_H;
		for (int l = 0; l < k; l++){ dy[l] = (l == kk) ? h : 0.0; }
		fn_shift(d, s, Nq, k, dy, sp);
		double xp = fn_xval(gv, d, sp, j);
		dy[kk] = -h;
		fn_shift(d, s, Nq, k, dy, sp);
		double xm = fn_xval(gv, d, sp, j);
		g[kk] = (xp - xm)/(2.0*h);
	}
	if (H == NULL){ return; }
	for (int kk = 0; kk < k*k; kk++){ H[kk] = 0.0; }
	if (d->ns_mode != 1 && d->ns_mode != 3){ return; }
	for (int a = 0; a < k; a++){
		for (int b = a; b < k; b++){
			double ha = FN_H2;
			double hb = FN_H2;
			double v[4];
			for (int q = 0; q < 4; q++){
				for (int l = 0; l < k; l++){ dy[l] = 0.0; }
				dy[a] += (q < 2) ? ha : -ha;
				dy[b] += (q % 2 == 0) ? hb : -hb;
				fn_shift(d, s, Nq, k, dy, sp);
				v[q] = fn_xval(gv, d, sp, j);
			}
			H[a*k + b] = (v[0] - v[1] - v[2] + v[3])/(4.0*ha*hb);
			H[b*k + a] = H[a*k + b];
		}
	}
}

static double fn_sf_min(			SS_ref 				*d,
									const double 		*x				){
	double xc[d->n_xeos];
	for (int j = 0; j < d->n_xeos; j++){ xc[j] = x[j]; }
	d->ns_obj(d->n_xeos, xc, NULL, d);
	double r = 1.0;
	for (int i = 0; i < d->n_sf; i++){
		if (!(d->sf[i] >= r)){ r = d->sf[i]; }
	}
	return r;
}

static int fn_sf_fit(				SS_ref 				*d,
									double 				*x,
									const double 		*x0				){
	int 	n_x = d->n_xeos;
	double 	xt[n_x];
	double 	dx  = 0.0;

	for (int j = 0; j < n_x; j++){
		if (fabs(x[j] - x0[j]) <= FN_SNAP_TOL){ x[j] = x0[j]; }
		if (fabs(x[j] - x0[j]) > dx){ dx = fabs(x[j] - x0[j]); }
	}
	if (fn_sf_min(d, x) >= 0.0){ return 1; }
	if (!(dx > 0.0)){ return 0; }

	double t = (FN_X_TOL/dx < 1.0) ? FN_X_TOL/dx : 1.0;
	for (int k = 0; k < 24; k++){
		double tk = t*pow(0.5, (double)(23 - k));
		for (int j = 0; j < n_x; j++){ xt[j] = (1.0 - tk)*x[j] + tk*x0[j]; }
		if (fn_sf_min(d, xt) >= 0.0){
			for (int j = 0; j < n_x; j++){ x[j] = xt[j]; }
			return 1;
		}
	}
	return 0;
}

static int fn_system(				global_variable 	 gv,
									bulk_info 	 		 z_b,
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp,
									int 				 n_ph,
									int 				 n_pp,
									const int 			*ph_cp,
									const int 			*pp_id,
									const int 			*off_y,
									const int 			*off_s,
									const int 			*off_N,
									const double 		*Nb,
									const int 			*off_c,
									const int 			*c_id,
									const int 			*cx_j,
									const double 		*cx_b,
									const double 		*gam,
									const double 		*s,
									const double 		*mu,
									const double 		*n,
									const int 			*drop,
									double 				*R,
									double 				*A				){
	int 	m    = z_b.nzEl_val;
	int 	n_y  = off_y[n_ph];
	int 	n_cc = off_c[n_ph];
	int 	nz   = m + n_y + n_cc + n_ph + n_pp;
	int 	rc   = n_y;
	int 	r2   = n_y + n_cc;
	int 	r3   = r2 + n_ph;
	int 	r4   = r3 + n_pp;
	int 	cm   = m + n_y;
	int 	cn   = m + n_y + n_cc;

	if (A != NULL){
		for (int i = 0; i < nz*nz; i++){ A[i] = 0.0; }
	}
	for (int jj = 0; jj < m; jj++){ R[r4 + jj] = -z_b.bulk_rock[z_b.nzEl_array[jj]]; }

	for (int q = 0; q < n_ph; q++){
		SS_ref 		 *d  = &SS_ref_db[cp[ph_cp[q]].id];
		int 		  nc = d->ns_nc;
		int 		  k  = off_y[q + 1] - off_y[q];
		const double *Nq = &Nb[off_N[q]];
		const double *sq = &s[off_s[q]];
		double 		  g[nc], gp[nc], gm[nc], sp[nc], dy[k > 0 ? k : 1];
		double 		  c[gv.len_ox], cpl[gv.len_ox], cml[gv.len_ox];
		double 		  cg[k > 0 ? k : 1], cH[k > 0 ? k*k : 1];

		double D = fn_eval(gv, d, gam, sq, g, c);
		if (isnan(D) || isinf(D)){ return 1; }
		for (int kk = 0; kk < k; kk++){
			double a = 0.0;
			for (int i = 0; i < nc; i++){ a += Nq[i*k + kk]*g[i]; }
			R[off_y[q] + kk] = a;
		}
		R[r2 + q] = (drop[q]) ? n[q] : D;
		for (int jj = 0; jj < m; jj++){ R[r4 + jj] += n[q]*c[z_b.nzEl_array[jj]]; }

		for (int ci = off_c[q]; ci < off_c[q + 1]; ci++){
			int id = c_id[ci];
			fn_xgrad(gv, d, sq, Nq, k, cx_j[id], cg, (A != NULL) ? cH : NULL);
			R[rc + ci] = fn_xval(gv, d, sq, cx_j[id]) - cx_b[id];
			for (int kk = 0; kk < k; kk++){ R[off_y[q] + kk] -= mu[ci]*cg[kk]; }
			if (A != NULL){
				for (int kk = 0; kk < k; kk++){
					A[(cm + ci)*nz + off_y[q] + kk] 	  = -cg[kk];
					A[(m + off_y[q] + kk)*nz + rc + ci] =  cg[kk];
					for (int ll = 0; ll < k; ll++){
						A[(m + off_y[q] + ll)*nz + off_y[q] + kk] -= mu[ci]*cH[kk*k + ll];
					}
				}
			}
		}

		if (A == NULL){ continue; }
		for (int jj = 0; jj < m; jj++){
			int ox = z_b.nzEl_array[jj];
			A[jj*nz + r2 + q] 		 = (drop[q]) ? 0.0 : -c[ox];
			A[(cn + q)*nz + r4 + jj] =  c[ox];
		}
		if (drop[q]){ A[(cn + q)*nz + r2 + q] = 1.0; }
		for (int kk = 0; kk < k; kk++){
			double h  = fn_step(d, sq, Nq, k, kk, FN_H);
			int    cy = m + off_y[q] + kk;
			if (!(h > 0.0)){ return 1; }

			for (int l = 0; l < k; l++){ dy[l] = (l == kk) ? h : 0.0; }
			fn_shift(d, sq, Nq, k, dy, sp);
			double Dp = fn_eval(gv, d, gam, sp, gp, cpl);
			dy[kk] = -h;
			fn_shift(d, sq, Nq, k, dy, sp);
			double Dm = fn_eval(gv, d, gam, sp, gm, cml);
			if (isnan(Dp) || isnan(Dm) || isinf(Dp) || isinf(Dm)){ return 1; }

			for (int ll = 0; ll < k; ll++){
				double a = 0.0;
				for (int i = 0; i < nc; i++){ a += Nq[i*k + ll]*(gp[i] - gm[i]); }
				A[cy*nz + off_y[q] + ll] += a/(2.0*h);
			}
			A[cy*nz + r2 + q] = (Dp - Dm)/(2.0*h);
			for (int jj = 0; jj < m; jj++){
				int    ox = z_b.nzEl_array[jj];
				double dc = (cpl[ox] - cml[ox])/(2.0*h);
				A[cy*nz + r4 + jj] 		 =  n[q]*dc;
				A[jj*nz + off_y[q] + kk] = -dc;
			}
		}
	}

	for (int q = 0; q < n_pp; q++){
		PP_ref *pp = &PP_ref_db[pp_id[q]];
		double  G  = pp->gbase*pp->factor;
		for (int jj = 0; jj < m; jj++){
			int    ox = z_b.nzEl_array[jj];
			double cq = pp->Comp[ox]*pp->factor;
			G 		   -= gam[ox]*cq;
			R[r4 + jj] += n[n_ph + q]*cq;
			if (A != NULL){
				A[jj*nz + r3 + q] 				= (drop[n_ph + q]) ? 0.0 : -cq;
				A[(cn + n_ph + q)*nz + r4 + jj] =  cq;
			}
		}
		R[r3 + q] = (drop[n_ph + q]) ? n[n_ph + q] : G;
		if (A != NULL && drop[n_ph + q]){ A[(cn + n_ph + q)*nz + r3 + q] = 1.0; }
	}

	return 0;
}

static double fn_G_system(			global_variable 	 gv,
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp,
									int 				 n_ph,
									int 				 n_pp,
									const int 			*ph_cp,
									const int 			*pp_id,
									const int 			*off_s,
									const double 		*s,
									const double 		*n				){
	double zero[gv.len_ox];
	double G = 0.0;

	for (int j = 0; j < gv.len_ox; j++){ zero[j] = 0.0; }
	for (int q = 0; q < n_ph; q++){
		G += n[q]*fn_eval(gv, &SS_ref_db[cp[ph_cp[q]].id], zero, &s[off_s[q]], NULL, NULL);
	}
	for (int q = 0; q < n_pp; q++){
		G += n[n_ph + q]*PP_ref_db[pp_id[q]].gbase*PP_ref_db[pp_id[q]].factor;
	}
	return G;
}

static void fn_restore_ss(			global_variable 	 gv,
									bulk_info 	 		 z_b,
									PC_type 			*PC_read,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp,
									int 				 n_ph,
									const int 			*ph_cp			){
	for (int q = 0; q < n_ph; q++){
		int i     = ph_cp[q];
		int ph_id = cp[i].id;

		fn_bounds(&SS_ref_db[ph_id], 0.0);
		SS_ref_db[ph_id] = rotate_hyperplane(		gv,
													SS_ref_db[ph_id]		);
		for (int k = 0; k < cp[i].n_xeos; k++){
			SS_ref_db[ph_id].iguess[k] = cp[i].xeos[k];
		}
		SS_ref_db[ph_id] = PC_function(				gv,
													PC_read,
													SS_ref_db[ph_id],
													z_b,
													ph_id					);
		SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
													SS_ref_db[ph_id],
													z_b,
													gv.SS_list[ph_id]		);
	}
}

static int fn_basis(				SS_ref 				*d,
									const double 		*sq,
									double 				*Nq				){
	int 	nc  = d->ns_nc;
	int 	k0  = d->ns_n_dir0;
	double 	v[k0 > 0 ? k0 : 1], Q[k0 > 0 ? k0*k0 : 1];
	int 	act[nc];
	int 	n_q = 0;

	for (int i = 0; i < nc; i++){
		double nr = 0.0;
		for (int kk = 0; kk < k0; kk++){ nr += fabs(d->ns_N0[i][kk]); }
		act[i] = (d->ns_sf_state[i] == 0 && nr > FN_N_TOL && (sq[i] <= FN_S_FIX || sq[i] >= d->ns_ub[i] - FN_S_FIX));
	}
	for (int i = 0; i < nc; i++){
		if (act[i] == 0){ continue; }
		for (int kk = 0; kk < k0; kk++){ v[kk] = d->ns_N0[i][kk]; }
		n_q += fn_orth_add(Q, n_q, k0, v);
	}
	int r_a = n_q;
	for (int kk = 0; kk < k0; kk++){
		for (int l = 0; l < k0; l++){ v[l] = (l == kk) ? 1.0 : 0.0; }
		n_q += fn_orth_add(Q, n_q, k0, v);
	}
	int k = n_q - r_a;
	for (int i = 0; i < nc; i++){
		for (int c = 0; c < k; c++){
			double a = 0.0;
			if (act[i] == 0 && d->ns_sf_state[i] == 0){
				for (int kk = 0; kk < k0; kk++){ a += d->ns_N0[i][kk]*Q[(r_a + c)*k0 + kk]; }
			}
			Nq[i*k + c] = a;
		}
	}
	return k;
}

static int fn_add_xbounds(			global_variable 	 gv,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp,
									int 				 n_ph,
									const int 			*ph_cp,
									const int 			*off_s,
									const double 		*s,
									double 				 tol,
									int 				*cx_q,
									int 				*cx_j,
									double 				*cx_b,
									double 				*cx_sg,
									int 				*cx_on,
									int 				*n_cx,
									int 				 n_max			){
	int n_add = 0;
	for (int q = 0; q < n_ph; q++){
		SS_ref *d = &SS_ref_db[cp[ph_cp[q]].id];
		ns_x_of_sf(gv, d, &s[off_s[q]]);
		for (int j = 0; j < d->n_xeos; j++){
			double lo   = d->bounds_ref[j][0];
			double hi   = d->bounds_ref[j][1];
			int    side = (d->ns_x[j] < lo + tol) ? 1 : (d->ns_x[j] > hi - tol) ? -1 : 0;
			if (side == 0 || !(hi - lo > 1e-6)){ continue; }
			int known = 0;
			for (int c = 0; c < *n_cx; c++){
				if (cx_q[c] == q && cx_j[c] == j){ known = (cx_on[c] == 1) ? 1 : 2; }
			}
			if (known != 0 || *n_cx >= n_max){ continue; }
			cx_q[*n_cx]  = q;
			cx_j[*n_cx]  = j;
			cx_b[*n_cx]  = (side == 1) ? lo : hi;
			cx_sg[*n_cx] = (double)side;
			cx_on[*n_cx] = 1;
			*n_cx 		+= 1;
			n_add 		+= 1;
		}
	}
	return n_add;
}

static int fn_kkt_reserve(			fn_kkt 				*k,
									int 				 nz,
									int 				 n_s,
									int 				 n_Nb			){
	if (nz > k->cap_nz){
		free(k->c_x);
		free(k->c_b);
		free(k->mu);
		free(k->R);
		free(k->J);
		k->c_x    = malloc ((nz) 					* sizeof(int)	);
		k->c_b    = malloc ((nz) 					* sizeof(double));
		k->mu     = malloc ((nz) 					* sizeof(double));
		k->R      = malloc ((nz) 					* sizeof(double));
		k->J      = malloc ((nz*nz) 				* sizeof(double));
		k->cap_nz = nz;
	}
	if (n_s > k->cap_s){
		free(k->s);
		k->s      = malloc ((n_s) 					* sizeof(double));
		k->cap_s  = n_s;
	}
	if (n_Nb > k->cap_Nb){
		free(k->Nb);
		k->Nb     = malloc ((n_Nb) 					* sizeof(double));
		k->cap_Nb = n_Nb;
	}
	if (k->c_x == NULL || k->c_b == NULL || k->mu == NULL || k->R == NULL || k->J == NULL || k->s == NULL || k->Nb == NULL){
		k->cap_nz = 0;
		k->cap_s  = 0;
		k->cap_Nb = 0;
		return 1;
	}
	return 0;
}

static int fn_kkt_store(			global_variable 	 gv,
									bulk_info 	 		 z_b,
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp,
									int 				 n_ph,
									int 				 n_pp,
									const int 			*ph_cp,
									const int 			*pp_id,
									const int 			*off_y,
									const int 			*off_s,
									const int 			*off_N,
									const double 		*Nb,
									const int 			*off_c,
									const int 			*c_id,
									const int 			*cx_j,
									const double 		*cx_b,
									const double 		*gam,
									const double 		*s,
									const double 		*mu,
									const double 		*n,
									const int 			*drop,
									double 				 G,
									fn_kkt 				*k				){
	int 	m    = z_b.nzEl_val;
	int 	n_y  = off_y[n_ph];
	int 	n_cc = off_c[n_ph];
	int 	nz   = m + n_y + n_cc + n_ph + n_pp;

	k->stored = 0;
	if (fn_kkt_reserve(k, nz, (off_s[n_ph] > 0) ? off_s[n_ph] : 1, (off_N[n_ph] > 0) ? off_N[n_ph] : 1) != 0){ return 1; }
	if (fn_system(gv, z_b, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_y, off_s, off_N, Nb, off_c, c_id, cx_j, cx_b, gam, s, mu, n, drop, k->R, k->J) != 0){ return 1; }

	k->nz   = nz;
	k->m    = m;
	k->n_y  = n_y;
	k->n_cc = n_cc;
	k->n_ph = n_ph;
	k->n_pp = n_pp;
	k->n_s  = off_s[n_ph];
	k->n_Nb = off_N[n_ph];
	k->G    = G;
	for (int jj = 0; jj < m; jj++){ k->ox_id[jj] = z_b.nzEl_array[jj]; }
	for (int j = 0; j < gv.len_ox; j++){ k->gam[j] = gam[j]; }
	for (int q = 0; q <= n_ph; q++){
		k->off_y[q] = off_y[q];
		k->off_s[q] = off_s[q];
		k->off_N[q] = off_N[q];
		k->off_c[q] = off_c[q];
	}
	for (int q = 0; q < n_ph; q++){
		k->cp_id[q] = ph_cp[q];
		k->ss_id[q] = cp[ph_cp[q]].id;
	}
	for (int q = 0; q < n_pp; q++){ k->pp_id[q] = pp_id[q]; }
	for (int q = 0; q < n_ph + n_pp; q++){
		k->n[q]    = n[q];
		k->drop[q] = drop[q];
	}
	for (int c = 0; c < n_cc; c++){
		k->c_x[c] = cx_j[c_id[c]];
		k->c_b[c] = cx_b[c_id[c]];
		k->mu[c]  = mu[c];
	}
	for (int i = 0; i < off_s[n_ph]; i++){ k->s[i] = s[i]; }
	for (int i = 0; i < off_N[n_ph]; i++){ k->Nb[i] = Nb[i]; }

	return 0;
}

global_variable final_Newton(		bulk_info 	 		 z_b,
									global_variable 	 gv,
									PC_type 			*PC_read,
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp				){
	gv.fn_status = 0;
	gv.fn_ite    = 0;
	if (gv.fn_sys != NULL){ gv.fn_sys->stored = 0; }
	if (gv.final_Newton_step != 1 || gv.fn_A == NULL){ return gv; }

	gv.fn_status = -1;
	if (gv.status != 0 || gv.solver == 3 || gv.n_mu_fix > 0){ return gv; }

	int 	m    = z_b.nzEl_val;
	int 	n_ph = 0;
	int 	n_pp = 0;
	int 	ph_cp[gv.len_ox];
	int 	pp_id[gv.len_ox];
	int 	off_d[gv.len_ox + 1];
	int 	off_s[gv.len_ox + 1];
	int 	off_x[gv.len_ox + 1];
	int 	off_N[gv.len_ox + 1];

	off_d[0] = 0;
	off_s[0] = 0;
	off_x[0] = 0;
	off_N[0] = 0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] != 1){ continue; }
		SS_ref *d = &SS_ref_db[cp[i].id];
		if (n_ph >= gv.len_ox || d->ns_ok != 1 || d->ns_fd != 0){ return gv; }
		ph_cp[n_ph]     = i;
		off_d[n_ph + 1] = off_d[n_ph] + d->ns_n_dir0;
		off_s[n_ph + 1] = off_s[n_ph] + d->ns_nc;
		off_x[n_ph + 1] = off_x[n_ph] + d->n_xeos;
		off_N[n_ph + 1] = off_N[n_ph] + d->ns_nc*d->ns_n_dir0;
		n_ph           += 1;
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] != 1){ continue; }
		if (n_pp >= gv.len_ox || gv.pp_flags[i][4] == 1){ return gv; }
		pp_id[n_pp] = i;
		n_pp       += 1;
	}
	if (n_ph == 0 || m + 2*off_d[n_ph] + n_ph + n_pp > gv.fn_nz_max){ return gv; }

	int 	n_dmax = off_d[n_ph];
	int 	off_y[n_ph + 1];
	int 	off_c[n_ph + 1];
	int 	c_id[n_dmax + 1];
	int 	cx_q[n_dmax + 1];
	int 	cx_j[n_dmax + 1];
	int 	cx_on[n_dmax + 1];
	double 	cx_b[n_dmax + 1];
	double 	cx_sg[n_dmax + 1];
	double 	Nb[off_N[n_ph] > 0 ? off_N[n_ph] : 1];
	double 	s0[off_s[n_ph]], s[off_s[n_ph]], st[off_s[n_ph]];
	double 	xs[off_x[n_ph]];
	double 	n0[n_ph + n_pp], n[n_ph + n_pp], nt[n_ph + n_pp];
	double 	gam0[gv.len_ox], gam[gv.len_ox], gamt[gv.len_ox];
	double 	mu[n_dmax + 1], mut[n_dmax + 1];
	int 	n_cx = 0;
	int 	ok   = 1;
	int 	set  = 1;

	for (int j = 0; j < gv.len_ox; j++){ gam0[j] = gv.gam_tot[j]; }
	for (int q = 0; q < n_ph; q++){ n0[q] = cp[ph_cp[q]].ss_n; }
	for (int q = 0; q < n_pp; q++){ n0[n_ph + q] = gv.pp_n[pp_id[q]]; }

	for (int q = 0; q < n_ph; q++){ fn_bounds(&SS_ref_db[cp[ph_cp[q]].id], 1.0); }

	for (int q = 0; q < n_ph && set; q++){
		SS_ref *d   = &SS_ref_db[cp[ph_cp[q]].id];
		int     nc  = d->ns_nc;
		int     k0  = d->ns_n_dir0;
		double *sq  = &s0[off_s[q]];
		double  x[d->n_xeos], r[nc], z[k0 > 0 ? k0 : 1];

		for (int j = 0; j < d->n_xeos; j++){ x[j] = cp[ph_cp[q]].xeos[j]; }
		d->ns_obj(d->n_xeos, x, NULL, d);
		for (int i = 0; i < nc; i++){
			r[i] = ((d->ns_mode == 3 || d->ns_mode == 4) ? ((i < d->n_em) ? d->p[i] : 0.0) : d->sf[i]) - d->ns_sf0[i];
		}
		for (int kk = 0; kk < k0; kk++){
			z[kk] = 0.0;
			for (int i = 0; i < nc; i++){ z[kk] += d->ns_N0[i][kk]*r[i]; }
		}
		for (int i = 0; i < nc; i++){
			double a = 0.0;
			for (int kk = 0; kk < k0; kk++){ a += d->ns_N0[i][kk]*z[kk]; }
			if (fabs(a - r[i]) > FN_PROJ_TOL){ set = 0; }
			sq[i] = (d->ns_sf_state[i] == 0) ? d->ns_sf0[i] + a : d->ns_sf0[i];
			if (sq[i] < -FN_S_FIX || sq[i] > d->ns_ub[i] + FN_S_FIX){ set = 0; }
		}
	}
	if (set){
		fn_add_xbounds(gv, SS_ref_db, cp, n_ph, ph_cp, off_s, s0, FN_X_TOL, cx_q, cx_j, cx_b, cx_sg, cx_on, &n_cx, n_dmax);
	}

	ok = set;
	double 	r0   = 0.0;
	double 	G0   = 0.0;
	double 	Gc   = 0.0;
	double 	G1   = 0.0;
	int 	ite  = 0;
	int 	done = 0;
	int 	n_y  = 0;
	double 	rmb  = 0.0;

	int 	drop[n_ph + n_pp], push[n_ph + n_pp];
	int 	n_drop = 0;
	int 	n_cx0  = n_cx;
	for (int q = 0; q < n_ph + n_pp; q++){ drop[q] = 0; }

	for (int dr = 0; dr <= FN_MAX_DROP && set; dr++){
		for (int q = 0; q < n_ph + n_pp; q++){ push[q] = 0; }
		n_cx = n_cx0;
		for (int c = 0; c < n_cx; c++){ cx_on[c] = (drop[cx_q[c]]) ? 0 : 1; }
		ok   = 1;
		done = 0;
		for (int j = 0; j < gv.len_ox; j++){ gam[j] = gam0[j]; }
		for (int i = 0; i < off_s[n_ph]; i++){ s[i] = s0[i]; }
		for (int q = 0; q < n_ph + n_pp; q++){ n[q] = (drop[q]) ? 0.0 : n0[q]; }

		for (int att = 0; att < FN_MAX_ATT && ok && done == 0; att++){
			int n_cc = 0;
			off_y[0] = 0;
			off_c[0] = 0;
			for (int q = 0; q < n_ph; q++){
				SS_ref *d  = &SS_ref_db[cp[ph_cp[q]].id];
				double *Nq = &Nb[off_N[q]];
				int     k  = (drop[q]) ? 0 : fn_basis(d, &s[off_s[q]], Nq);
				double  cg[k > 0 ? k : 1], Qx[k > 0 ? k*k : 1];
				int     n_qx = 0;

				off_y[q + 1] = off_y[q] + k;
				for (int c = 0; c < n_cx && k > 0; c++){
					if (cx_q[c] != q || cx_on[c] != 1){ continue; }
					fn_xgrad(gv, d, &s[off_s[q]], Nq, k, cx_j[c], cg, NULL);
					if (fn_orth_add(Qx, n_qx, k, cg) == 0){ continue; }
					n_qx 	   += 1;
					c_id[n_cc]  = c;
					n_cc 	   += 1;
				}
				off_c[q + 1] = n_cc;
			}
			n_y    = off_y[n_ph];
			int nz = m + n_y + n_cc + n_ph + n_pp;
			int cm = m + n_y;
			int cn = m + n_y + n_cc;
			if (nz > gv.fn_nz_max){ ok = 0; break; }

			double R[nz], Rt[nz];
			for (int c = 0; c < n_cc; c++){ mu[c] = 0.0; }

			if (fn_system(gv, z_b, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_y, off_s, off_N, Nb, off_c, c_id, cx_j, cx_b, gam, s, mu, n, drop, R, NULL) != 0){ ok = 0; break; }
			r0 = fn_norm(R, nz);
			if (att == 0 && dr == 0){
				G0 = fn_G_system(gv, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_s, s, n);
				Gc = G0;
				for (int jj = 0; jj < m; jj++){
					double rj = R[n_y + n_cc + n_ph + n_pp + jj];
					Gc  -= gam[z_b.nzEl_array[jj]]*rj;
					rmb += rj*rj;
				}
			}
			ite = 0;

			while (ok && r0 > FN_TOL && ite < FN_MAX_ITE){
				if (fn_system(gv, z_b, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_y, off_s, off_N, Nb, off_c, c_id, cx_j, cx_b, gam, s, mu, n, drop, R, gv.fn_A) != 0){ ok = 0; break; }
				for (int i = 0; i < nz; i++){ gv.fn_b[i] = -R[i]; }

				int 	nrhs = 1;
				int 	lda  = nz;
				int 	info = 0;
				#if __APPLE__
					dgesv(&nz, &nrhs, gv.fn_A, &lda, gv.fn_ipiv, gv.fn_b, &lda, &info);
				#else
					info = LAPACKE_dgesv(LAPACK_COL_MAJOR, nz, nrhs, gv.fn_A, lda, gv.fn_ipiv, gv.fn_b, lda);
				#endif
				if (info != 0){ ok = 0; break; }

				double alpha = 1.0;
				for (int q = 0; q < n_ph; q++){
					SS_ref 		 *d  = &SS_ref_db[cp[ph_cp[q]].id];
					int 		  k  = off_y[q + 1] - off_y[q];
					const double *Nq = &Nb[off_N[q]];
					const double *dy = &gv.fn_b[m + off_y[q]];
					for (int i = 0; i < d->ns_nc; i++){
						double ds = 0.0;
						for (int kk = 0; kk < k; kk++){ ds += Nq[i*k + kk]*dy[kk]; }
						double si = s[off_s[q] + i];
						if (ds < 0.0 && FN_FRAC_S*si/(-ds) < alpha){ alpha = FN_FRAC_S*si/(-ds); }
						if (ds > 0.0 && FN_FRAC_S*(d->ns_ub[i] - si)/ds < alpha){ alpha = FN_FRAC_S*(d->ns_ub[i] - si)/ds; }
					}
				}
				for (int q = 0; q < n_ph + n_pp; q++){
					double dn = gv.fn_b[cn + q];
					if (drop[q]){ continue; }
					if (dn < 0.0 && n[q] + dn <= 0.0){ push[q] += 1; }
					if (dn < 0.0 && FN_FRAC*n[q]/(-dn) < alpha){ alpha = FN_FRAC*n[q]/(-dn); }
				}

				int 	acc = 0;
				double 	r1  = r0;
				for (int ls = 0; ls < FN_MAX_LS && acc == 0; ls++){
					for (int j = 0; j < gv.len_ox; j++){ gamt[j] = gam[j]; }
					for (int jj = 0; jj < m; jj++){ gamt[z_b.nzEl_array[jj]] += alpha*gv.fn_b[jj]; }
					for (int q = 0; q < n_ph; q++){
						SS_ref *d = &SS_ref_db[cp[ph_cp[q]].id];
						int     k = off_y[q + 1] - off_y[q];
						double  dy[k > 0 ? k : 1];
						for (int kk = 0; kk < k; kk++){ dy[kk] = alpha*gv.fn_b[m + off_y[q] + kk]; }
						fn_shift(d, &s[off_s[q]], &Nb[off_N[q]], k, dy, &st[off_s[q]]);
					}
					for (int c = 0; c < n_cc; c++){ mut[c] = mu[c] + alpha*gv.fn_b[cm + c]; }
					for (int q = 0; q < n_ph + n_pp; q++){ nt[q] = n[q] + alpha*gv.fn_b[cn + q]; }

					if (fn_system(gv, z_b, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_y, off_s, off_N, Nb, off_c, c_id, cx_j, cx_b, gamt, st, mut, nt, drop, Rt, NULL) == 0){
						r1 = fn_norm(Rt, nz);
						if (r1 < r0 || r1 <= FN_TOL){ acc = 1; }
					}
					if (acc == 0){ alpha *= 0.5; }
				}
				if (acc == 0){ ok = 0; break; }

				for (int j = 0; j < gv.len_ox; j++){ gam[j] = gamt[j]; }
				for (int i = 0; i < off_s[n_ph]; i++){ s[i] = st[i]; }
				for (int c = 0; c < n_cc; c++){ mu[c] = mut[c]; }
				for (int q = 0; q < n_ph + n_pp; q++){ n[q] = nt[q]; }
				r0   = r1;
				ite += 1;
			}
			if (r0 > FN_TOL){ ok = 0; }
			if (!ok){
				int n_new = fn_add_xbounds(gv, SS_ref_db, cp, n_ph, ph_cp, off_s, s0, FN_X_ACT, cx_q, cx_j, cx_b, cx_sg, cx_on, &n_cx, n_dmax)
						  + fn_add_xbounds(gv, SS_ref_db, cp, n_ph, ph_cp, off_s, s, FN_X_ACT, cx_q, cx_j, cx_b, cx_sg, cx_on, &n_cx, n_dmax);
				int n_fix = 0;
				for (int q = 0; q < n_ph; q++){
					SS_ref *d = &SS_ref_db[cp[ph_cp[q]].id];
					int     k = off_y[q + 1] - off_y[q];
					for (int i = 0; i < d->ns_nc; i++){
						double v = 0.0;
						for (int kk = 0; kk < k; kk++){ v += fabs(Nb[off_N[q] + i*k + kk]); }
						double si = s[off_s[q] + i];
						if (v > 0.0 && (si <= FN_S_FIX || si >= d->ns_ub[i] - FN_S_FIX)){ n_fix += 1; }
					}
				}
				if (n_new > 0 && att < FN_MAX_ATT - 1){
					for (int j = 0; j < gv.len_ox; j++){ gam[j] = gam0[j]; }
					for (int i = 0; i < off_s[n_ph]; i++){ s[i] = s0[i]; }
					for (int q = 0; q < n_ph + n_pp; q++){ n[q] = (drop[q]) ? 0.0 : n0[q]; }
					ok = 1;
					continue;
				}
				if (n_fix > 0 && att < FN_MAX_ATT - 1){ ok = 1; continue; }
				break;
			}

			int n_rel = 0;
			int n_add = 0;
			for (int c = 0; c < n_cc; c++){
				if (cx_sg[c_id[c]]*mu[c] < -FN_MU_TOL){ cx_on[c_id[c]] = 0; n_rel += 1; }
			}
			if (n_rel == 0){
				n_add = fn_add_xbounds(gv, SS_ref_db, cp, n_ph, ph_cp, off_s, s, -FN_X_TOL, cx_q, cx_j, cx_b, cx_sg, cx_on, &n_cx, n_dmax);
			}
			if (n_rel == 0 && n_add == 0){ done = 1; }
			else if (att == FN_MAX_ATT - 1){ ok = 0; }
			else{
				for (int j = 0; j < gv.len_ox; j++){ gam[j] = gam0[j]; }
				for (int i = 0; i < off_s[n_ph]; i++){ s[i] = s0[i]; }
				for (int q = 0; q < n_ph + n_pp; q++){ n[q] = (drop[q]) ? 0.0 : n0[q]; }
			}
		}
		if (done == 0){ ok = 0; }
		if (ok){
			for (int q = 0; q < n_ph + n_pp; q++){
				if (drop[q] == 0 && !(n[q] > 0.0)){ push[q] += 1; ok = 0; }
			}
		}
		if (ok){ break; }

		int qd = -1;
		for (int q = 0; q < n_ph + n_pp; q++){
			if (drop[q] == 1 || push[q] == 0){ continue; }
			if (qd < 0 || push[q] > push[qd] || (push[q] == push[qd] && n[q]/n0[q] < n[qd]/n0[qd])){ qd = q; }
		}
		if (qd < 0 || dr == FN_MAX_DROP){ break; }
		drop[qd] = 1;
		n_drop 	+= 1;
	}

	if (ok){
		G1 = fn_G_system(gv, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_s, s, n);
		if (!(G1 <= Gc + FN_G_TOL + FN_G_MB*rmb)){ ok = 0; }
	}
	if (ok){
		for (int q = 0; q < n_ph && ok; q++){
			if (drop[q]){ continue; }
			SS_ref *d = &SS_ref_db[cp[ph_cp[q]].id];
			ns_x_of_sf(gv, d, &s[off_s[q]]);
			for (int j = 0; j < d->n_xeos; j++){
				if (d->ns_x[j] < d->bounds_ref[j][0] - FN_X_TOL || d->ns_x[j] > d->bounds_ref[j][1] + FN_X_TOL){ ok = 0; }
				xs[off_x[q] + j] = fmin(fmax(d->ns_x[j], d->bounds_ref[j][0]), d->bounds_ref[j][1]);
			}
			if (ok && fn_sf_fit(d, &xs[off_x[q]], cp[ph_cp[q]].xeos) == 0){ ok = 0; }
		}
	}

	int kkt = 0;
	if (ok && gv.fn_store_kkt == 1 && gv.fn_sys != NULL){
		kkt = (fn_kkt_store(gv, z_b, PP_ref_db, SS_ref_db, cp, n_ph, n_pp, ph_cp, pp_id, off_y, off_s, off_N, Nb, off_c, c_id, cx_j, cx_b, gam, s, mu, n, drop, G1, gv.fn_sys) == 0);
	}

	for (int q = 0; q < n_ph + n_pp && ok && n_drop > 0; q++){
		if (drop[q] == 0){ continue; }
		if (q < n_ph){
			int 	i     = ph_cp[q];
			int 	ph_id = cp[i].id;
			SS_ref *d     = &SS_ref_db[ph_id];
			for (int k = 0; k < d->n_em; k++){
				d->gb_lvl[k] = d->gbase[k];
				for (int j = 0; j < gv.len_ox; j++){ d->gb_lvl[k] -= d->Comp[k][j]*gam[j]; }
			}
			for (int k = 0; k < cp[i].n_xeos; k++){ d->iguess[k] = cp[i].xeos[k]; }
			SS_ref_db[ph_id] = NS_opt_function(		gv,
														SS_ref_db[ph_id]		);
			if (SS_ref_db[ph_id].ns_status != 3 || !(SS_ref_db[ph_id].df >= -FN_DF_TOL)){ ok = 0; }
		}
		else{
			PP_ref *pp = &PP_ref_db[pp_id[q - n_ph]];
			double  D  = pp->gbase*pp->factor;
			for (int j = 0; j < gv.len_ox; j++){ D -= gam[j]*pp->Comp[j]*pp->factor; }
			if (!(D >= -FN_DF_TOL)){ ok = 0; }
		}
	}

	if (ok){
		for (int j = 0; j < gv.len_ox; j++){ gv.gam_tot[j] = gam[j]; }
		for (int pass = 0; pass < 2 && ok; pass++){
			for (int q = 0; q < n_ph; q++){
				if (drop[q]){ continue; }
				int i     = ph_cp[q];
				int ph_id = cp[i].id;

				fn_bounds(&SS_ref_db[ph_id], 0.0);
				SS_ref_db[ph_id] = rotate_hyperplane(		gv,
															SS_ref_db[ph_id]		);
				for (int k = 0; k < cp[i].n_xeos; k++){
					SS_ref_db[ph_id].iguess[k] = xs[off_x[q] + k];
				}
				SS_ref_db[ph_id] = PC_function(				gv,
															PC_read,
															SS_ref_db[ph_id],
															z_b,
															ph_id					);
				SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
															SS_ref_db[ph_id],
															z_b,
															gv.SS_list[ph_id]		);
				if (SS_ref_db[ph_id].sf_ok != 1){ ok = 0; break; }
				if (pass == 1){
					copy_to_cp(								i,
															ph_id,
															gv,
															SS_ref_db,
															cp						);
					cp[i].ss_n = n[q];
				}
			}
		}
		if (ok){
			for (int q = 0; q < n_pp; q++){ gv.pp_n[pp_id[q]] = n[n_ph + q]; }
			for (int q = 0; q < n_ph + n_pp; q++){
				if (drop[q] == 0){ continue; }
				if (q < n_ph){
					cp[ph_cp[q]].ss_flags[1] = 0;
					cp[ph_cp[q]].ss_flags[2] = 1;
					cp[ph_cp[q]].ss_n 		 = 0.0;
					gv.n_cp_phase 			-= 1;
				}
				else{
					gv.pp_flags[pp_id[q - n_ph]][1] = 0;
					gv.pp_flags[pp_id[q - n_ph]][2] = 1;
					gv.pp_n[pp_id[q - n_ph]] 		= 0.0;
					gv.n_pp_phase 				   -= 1;
				}
				gv.n_phase -= 1;
			}
		}
		else{
			for (int j = 0; j < gv.len_ox; j++){ gv.gam_tot[j] = gam0[j]; }
		}
	}

	if (ok){
		int LP0  = gv.LP;
		int PGE0 = gv.PGE;
		gv.LP    = 1;
		gv.PGE   = 0;
		gv = PGE_residual_update(		z_b,
										gv,
										PP_ref_db,
										SS_ref_db,
										cp					);
		gv.LP    = LP0;
		gv.PGE   = PGE0;
		gv.fn_status = (n_drop > 0) ? 2 : 1;
		if (kkt){ gv.fn_sys->stored = 1; }
	}
	else{
		fn_restore_ss(gv, z_b, PC_read, SS_ref_db, cp, n_ph, ph_cp);
		gv.fn_status = (set) ? -2 : -1;
	}
	gv.fn_ite = ite;

	if (gv.verbose == 1){
		printf("\n Final Newton step: status %d, %d iterations, %d bound constraints, %d phase(s) removed, |R| %.3e, dG %+.3e\n", gv.fn_status, ite, n_cx, n_drop, r0, G1 - Gc);
	}

	return gv;
}

/**
  function to run simplex linear programming during PGE with pseudocompounds
*/
global_variable LP_pc_composite(					bulk_info 			 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													PC_type 			*PC_read,
													P2X_type 			*P2X_read,

													obj_type 			*SS_objective,
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int i, j, k, l, m, em_id, pc_id, ph_id, n_xeos, n_em;
	int nOcc, m_Ppc;

	double factor_mean, factor_composite, G;
	double p0, p1;

	double sum_n_vec 	= 0.0;
	double sum_n_vec_cor= 0.0;

	if (gv.verbose == 1){
		printf("\nPseudocompounds collapse (intermediate stage) \n");
		printf("══════════════════════════════════════════════\n");
	}

	/* loops through active solution phases and store their information */
	for (ph_id = 0; ph_id < gv.len_ss; ph_id++){
		if (SS_ref_db[ph_id].ss_flags[0] == 1/* && strcmp(gv.SS_list[ph_id],"liq") == 0*/){
			gv.n_ss_ph[ph_id] 	= 0;
			sum_n_vec 			= 0.0;
			sum_n_vec_cor 		= 0.0;
			nOcc 				= 0;

			/* first we retrieve the indexes of the solution phase */
			for (i = 0; i < d->n_Ox; i++){
				if (d->ph_id_A[i][0] == 2 || d->ph_id_A[i][0] == 3){
					k 		= d->ph_id_A[i][1];
					if (ph_id == k){
						gv.pc_id[nOcc]   = i;
						nOcc 			+= 1;
					}
				}
			}
			gv.n_ss_ph[ph_id] = nOcc;

			if (nOcc > 1){

				/* get unrotated gbase */
				SS_ref_db[ph_id] = non_rot_hyperplane(		gv, 
															SS_ref_db[ph_id]		);
				n_xeos 		=  SS_ref_db[ph_id].n_xeos;
				n_em 		=  SS_ref_db[ph_id].n_em;
				
				/* get information of the pseudocompounds for ph_id */
				for (i = 0; i < nOcc; i++){
					k = gv.pc_id[i];
					if (d->ph_id_A[k][0] == 2){
						em_id 			= d->ph_id_A[k][3];

						for (j = 0; j < n_em; j++) {	
							SS_ref_db[ph_id].p[j] = gv.em2ss_shift;
						}
						SS_ref_db[ph_id].p[em_id] = 1.0 - gv.em2ss_shift*n_em;
						
						(*P2X_read[ph_id])(		&SS_ref_db[ph_id],
												gv.bnd_val					);

						G 	= (*SS_objective[ph_id])(SS_ref_db[ph_id].n_xeos, SS_ref_db[ph_id].iguess, 	NULL, &SS_ref_db[ph_id]);

						for (j = 0; j < n_xeos; j++){
							gv.A[i][j] = SS_ref_db[ph_id].iguess[j];
						}

						gv.b[i] 		 = d->n_vec[k];
						gv.tmp1[i] 	     = SS_ref_db[ph_id].factor;
						sum_n_vec_cor 	+= gv.b1[i];
						sum_n_vec 		+= d->n_vec[k];
					}			
					if (d->ph_id_A[k][0] == 3 && d->stage[k] == 1){
						pc_id 			 = d->ph_id_A[k][3];

						for (j = 0; j < n_xeos; j++){
							SS_ref_db[ph_id].iguess[j] = SS_ref_db[ph_id].xeos_Ppc[pc_id][j];
						}
						/* then compute the normalization factor for un-corrected xeos */
						G 	= (*SS_objective[ph_id])(SS_ref_db[ph_id].n_xeos, SS_ref_db[ph_id].iguess, 	NULL, &SS_ref_db[ph_id]);

						for (j = 0; j < n_xeos; j++){
							gv.A[i][j]  = SS_ref_db[ph_id].iguess[j];
						}

						gv.b[i] 		 = d->n_vec[k];
						gv.tmp1[i] 	     = SS_ref_db[ph_id].factor;
						sum_n_vec_cor 	+= gv.b1[i];
						sum_n_vec 		+= d->n_vec[k];
					}
					if (d->ph_id_A[k][0] == 3 && d->stage[k] == 0){
						pc_id 			 = d->ph_id_A[k][3];

						for (j = 0; j < n_xeos; j++){
							SS_ref_db[ph_id].iguess[j] = SS_ref_db[ph_id].xeos_pc[pc_id][j];
						}
						/* then compute the normalization factor for un-corrected xeos */
						G 	= (*SS_objective[ph_id])(SS_ref_db[ph_id].n_xeos, SS_ref_db[ph_id].iguess, 	NULL, &SS_ref_db[ph_id]);

						for (j = 0; j < n_xeos; j++){
							gv.A[i][j]  = SS_ref_db[ph_id].iguess[j];
						}
						
						gv.b[i] 		 = d->n_vec[k];
						gv.tmp1[i] 	     = SS_ref_db[ph_id].factor;
						sum_n_vec_cor 	+= gv.b1[i];
						sum_n_vec 		+= d->n_vec[k];

					}
				}

				// first reset arrays
				for (j = 0; j < n_xeos; j++){
					SS_ref_db[ph_id].iguess[j] = 0.0;
				}

				/* retrieve initial guess */
				for (i = 0; i < nOcc; i++){
					gv.b[i] 	/= sum_n_vec;
					for (j = 0; j < n_xeos; j++){
						SS_ref_db[ph_id].iguess[j] += gv.A[i][j]*gv.b[i];
					}
				}

				/* retrieve normalization factor for correcting fraction vector */
				G 	= (*SS_objective[ph_id])(SS_ref_db[ph_id].n_xeos, SS_ref_db[ph_id].iguess, 	NULL, &SS_ref_db[ph_id]);

				factor_mean = SS_ref_db[ph_id].factor;

				/* correct pseudocompounds factors */
				sum_n_vec = 0.0;
				for (i = 0; i < nOcc; i++){
					gv.b1[i]    = gv.b[i]*gv.tmp1[i]*factor_mean;
					sum_n_vec  += gv.b1[i];
				}

				/* normalized corrected pseudocompounds fractions */
				for (i = 0; i < nOcc; i++){
					gv.b1[i] /= sum_n_vec;
				}

				/* retrieve corrected initial guess (stored in tmp2) */
				for (j = 0; j < n_xeos; j++){
					gv.tmp2[j] = 0.0;
				}
				for (i = 0; i < nOcc; i++){
					for (j = 0; j < n_xeos; j++){
						gv.tmp2[j] += gv.A[i][j]*gv.b1[i];
					}
				}

				/* At this stage: 
					n_vec_uncor are stored in gv.b 			[nOcc]

					xeos_ini 	are stored in gv.A			[nOcc]
					factor_ini  are stored in gv.tmp1 		[nOcc]

					xeos_mean 	is  stored in gv.tmp2		[n_xeos]
					factor_mean is  stored in factor_mean	[scalar]

				*/
				for (i = 0; i < nOcc; i++){

					/* First get un-corrected composite xeos */
					for (j = 0; j < n_xeos; j++){
						SS_ref_db[ph_id].iguess[j] = gv.A[i][j]*(gv.pc_composite_dist) + gv.tmp2[j] * (1.0 - gv.pc_composite_dist);
					}
					/* then compute the normalization factor for un-corrected xeos */
					G 	= (*SS_objective[ph_id])(SS_ref_db[ph_id].n_xeos, SS_ref_db[ph_id].iguess, 	NULL, &SS_ref_db[ph_id]);

					factor_composite = SS_ref_db[ph_id].factor;

					/* Compute corrected composite xeos */
					p0 			= gv.pc_composite_dist*(gv.tmp1[i] * factor_composite);
					p1 			= (1.0 - gv.pc_composite_dist)*(factor_mean * factor_composite);

					sum_n_vec 	= p0+p1;
					p0 		   /= sum_n_vec;
					p1 		   /= sum_n_vec;

					for (j = 0; j < n_xeos; j++){
						SS_ref_db[ph_id].iguess[j] = gv.A[i][j]*p0 + gv.tmp2[j]*p1;
						gv.A2[i][j] = SS_ref_db[ph_id].iguess[j];
					}

					SS_ref_db[ph_id] = PC_function(				gv,
																PC_read,
																SS_ref_db[ph_id], 
																z_b,
																ph_id 					);

					SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																SS_ref_db[ph_id], 
																z_b, 
																gv.SS_list[ph_id]		);

					copy_to_Ppc(								0,
																1,
																ph_id,
																gv,

																SS_objective,
																SS_ref_db				);	
				}
			}
		
		}
	}


	return gv;
}

/**
  Main LP routine
*/ 
global_variable LP(		bulk_info 			z_b,
						global_variable 	gv,
						PC_type				*PC_read,
						P2X_type			*P2X_read,

						obj_type 			*SS_objective,
						NLopt_type			*NLopt_opt,
						simplex_data	    *splx_data,
						PP_ref 				*PP_ref_db,
						SS_ref 				*SS_ref_db,
						csd_phase_set  		*cp					){
		
	clock_t t; 	

	gv.LP 	 = 1;	
	gv.PGE 	 = 0;

	int    mode = 0;
	int    gi   = 0;
	int iterate = 1;

	gv.PC_checked = 0;
	gv = init_LP(			z_b,
							splx_data,
							gv,
							PC_read,
							P2X_read,
									
							PP_ref_db,
							SS_ref_db,
							cp	);

	while (iterate == 1){
		

		t = clock();

		/* gv.global_ite > 1 must be checked FIRST: && short-circuits left-to-right in C, so
		   the original ordering (bound check last) still evaluated gv.gamma_norm[global_ite-1]
		   on every entry to this loop, including the very first one for a point
		   (global_ite==0), reading gv.gamma_norm[-1] - a heap-buffer-overflow confirmed via
		   AddressSanitizer (PGE_function.c:1533, one double before the gv.gamma_norm
		   allocation in TC_init_database.c). Usually silent (reads adjacent heap bytes), but
		   segfaults whenever the allocator happens to place gamma_norm at the start of a
		   fresh page - a heap-layout coincidence, not tied to any particular database/oxide
		   chemistry, though richer speciation (e.g. mpe with both S and CO2 active) changes
		   allocation patterns enough to shift how often that coincidence is hit. */
		if ((gv.global_ite > 1 && gv.gamma_norm[gv.global_ite-1] < 1.0 && gv.PC_checked < 2)){
			gv.PC_checked += 1;
			if (gv.verbose == 1){
				printf(" Checking PC for re-introduction:\n");
				printf(" ════════════════════════════════\n");
			}
			gv = check_PC( 				z_b,						/** bulk rock constraint 				*/ 
										gv,							/** global variables (e.g. Gamma) 		*/
										PC_read,

										PP_ref_db,					/** pure phase database 				*/ 
										SS_ref_db,
										cp					); 
			if (gv.verbose == 1){
				printf("\n");
			}

		}

		if (gv.verbose == 1){
			printf("\n__________________________________________ ‿︵MAGEMin‿︵ "); printf("_ %5s _",gv.version);
			printf("\n                     GLOBAL ITERATION %i\n",gv.global_ite);
			printf("═════════════════════════════════════════════════════════════════\n");
			printf("\nMinimize solution phases\n");
			printf("═════════════════════════\n");
			printf(" phase |  delta_G   | SF |   sum_xi   | time(ms)   |   x-eos ...\n");
			printf("══════════════════════════════════════════════════════════════════\n");
		}

		/** 
			update delta_G of pure phases as function of updated Gamma
		*/
		pp_min_function(				gv,
										z_b,
										PP_ref_db			);		

		/**
			Local minimization of the solution phases
		*/
		ss_min_LP(						gv, 							/** global variables (e.g. Gamma) 		*/
										PC_read,

										SS_objective,	
										NLopt_opt,						
										z_b,							/** bulk-rock, pressure and temperature conditions */
										SS_ref_db,						/** solution phase database 			*/	
										cp 					);

		/**
		   Here the linear programming method is used after the PGE step to get a new Gibbs hyper-plane
		*/
		gv = run_LP(					z_b,
										splx_data,
										gv,
												
										PP_ref_db,
										SS_ref_db			);

		gv = init_LP(					z_b,
										splx_data,
										gv,
										PC_read,
										P2X_read,
										
										PP_ref_db,
										SS_ref_db,
										cp					);	

		gv = compute_xi_SD(				gv,
										cp					);

		if (gv.verbose == 1){
			/* Partitioning Gibbs Energy */
			PGE_print(					z_b,							/** bulk rock constraint 				*/ 
										gv,								/** global variables (e.g. Gamma) 		*/

										PP_ref_db,						/** pure phase database 				*/ 
										SS_ref_db,						/** solution phase database 			*/
										cp					); 
		}

		/** 
			Update mass constraint residual
		*/
		gv = PGE_residual_update(		z_b,							/** bulk rock constraint 				*/ 
										gv,								/** global variables (e.g. Gamma) 		*/

										PP_ref_db,						/** pure phase database 				*/ 
										SS_ref_db,						/** solution phase database 			*/
										cp					);  
		
		/* Increment global iteration value */
		gv.global_ite += 1;

		/* PGE_mass_norm/Alg/gamma_norm/gibbs_ev/ite_time are allocated at gv.it_f*2 entries
		   (TC_init_database.c) - the convergence checks further down only bound global_ite
		   against gv.it_f (half that), on the assumption this loop and the "PGE proper"
		   loop below each stay within that budget. A genuinely hard/slow point (confirmed
		   via AddressSanitizer: ume 1kbar/850C) can still drive global_ite past the array's
		   actual gv.it_f*2 capacity before either soft check fires, overflowing these
		   arrays. Hard-clamp here so that case reports clean non-convergence instead of
		   corrupting the heap. */
		if (gv.global_ite >= gv.it_f*2){
			gv.global_ite = gv.it_f*2 - 1;
			gv.div        = 1;
			gv.status     = 4;
			iterate       = 0;
			if (gv.verbose == 1){
				printf(" >%d iterations (hard cap), not diverging but not converging\n\n", gv.it_f*2);
			}
			break;
		}

		/* check evolution of mass constraint residual */
		gv.PGE_mass_norm[gv.global_ite]  = gv.BR_norm;	/** save norm for the current global iteration */
		gv.Alg[gv.global_ite] 			 = 0;
		t 								 = clock() - t;

		if (gv.verbose == 1){
			printf("\n __ iteration duration: %+4f ms __\n\n\n",((double)t)/CLOCKS_PER_SEC*1000);
		}
		gv.ite_time[gv.global_ite] 		 = ((double)t)/CLOCKS_PER_SEC*1000;
		gi += 1;

		if ((gv.gamma_norm[gv.global_ite-1] < 1e-6 || gi >= gv.max_LP_ite) && gv.PC_checked == 2){
			iterate = 0;

			if (gv.gamma_norm[gv.global_ite-1] < 1e-6){
				gv.status = 0;
			}
			if (gi >= gv.max_LP_ite){
				if (gv.gamma_norm[gv.global_ite-1] < 1e-2){
					gv.status = 1;
				}
				else if ( gv.gamma_norm[gv.global_ite-1] >= 1e-2 && gv.gamma_norm[gv.global_ite-1] < 0.1){
					gv.status = 2;
				}
				else if ( gv.gamma_norm[gv.global_ite-1] >= 0.1 && gv.gamma_norm[gv.global_ite-1] < 1.0){
					gv.status = 3;
				}
				else if ( gv.gamma_norm[gv.global_ite-1] >= 1.0 && gv.gamma_norm[gv.global_ite-1] < 10.0){
					gv.status = 4;
				}
			}

			if (gv.BR_norm > 1e-3){
				gv.status = -1;
			}

			for (int i = 0; i < z_b.nzEl_val; i++){
				if (gv.gam_tot[z_b.nzEl_array[i]] >= 0.0 && strcmp(gv.ox[z_b.nzEl_array[i]], "Fe") != 0 && strcmp(gv.research_group, "tc") == 0){ gv.status = -1; }
			}

		}
	}

	/**
		Merge instances of the same solution phase that are compositionnally close 
	*/
	gv = phase_merge_function(		z_b,							/** bulk rock constraint 				*/
									gv,								/** global variables (e.g. Gamma) 		*/

									PP_ref_db,						/** pure phase database 				*/
									SS_ref_db,						/** solution phase database 			*/ 
									cp					); 

	gv = update_cp_after_LP(		z_b,
									gv,
									PC_read,

									PP_ref_db,
									SS_ref_db,
									cp					);

	return gv;
};



/**
  Main LP routine
*/ 
global_variable LP_metastable(	bulk_info 			z_b,
								global_variable 	gv,
								PC_type				*PC_read,
								P2X_type			*P2X_read,

								obj_type 			*SS_objective,
								NLopt_type			*NLopt_opt,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					){
		
	clock_t t; 	

	gv.LP 	 = 1;	
	gv.PGE 	 = 0;

	gv = init_LP(					z_b,
									splx_data,
									gv,
									PC_read,
									P2X_read,
											
									PP_ref_db,
									SS_ref_db,
									cp					);
	/**
		Merge instances of the same solution phase that are compositionnally close 
	*/
	gv = phase_merge_function(		z_b,							/** bulk rock constraint 				*/
									gv,								/** global variables (e.g. Gamma) 		*/

									PP_ref_db,						/** pure phase database 				*/
									SS_ref_db,						/** solution phase database 			*/ 
									cp					); 

	gv = update_cp_after_LP(		z_b,
									gv,
									PC_read,

									PP_ref_db,
									SS_ref_db,
									cp					);

	return gv;
};




/**
  Main PGE routine
*/ 
global_variable PGE(	bulk_info 			z_b,
						global_variable 	gv,
						PC_type             *PC_read,

						obj_type 			*SS_objective,
						NLopt_type 			*NLopt_opt,
						simplex_data	    *splx_data,
						PP_ref 				*PP_ref_db,
						SS_ref 				*SS_ref_db,
						csd_phase_set  		*cp					){
		
	clock_t t, v; 	

	gv.LP 			 = 0;	
	gv.PGE 			 = 1;

	int mode 	   = 0;
	int iterate    = 1;
	int pc_checked = 0;

	while (iterate == 1){
		gv.PC_checked = 0;
		pc_checked = 0;
		t = clock();
		if (gv.verbose == 1){
			printf("\n__________________________________________ ‿︵MAGEMin‿︵ "); printf("_ %5s _",gv.version);
			printf("\n                     GLOBAL ITERATION %i\n",gv.global_ite);
			printf("═════════════════════════════════════════════════════════════════\n");
		}
		
		/* calculate delta_G of solution phases (including local minimization) */
		if (gv.verbose == 1){
			printf("\nMinimize solution phases\n");
			printf("═════════════════════════\n");
			printf(" phase |  delta_G   | SF |   sum_xi   | time(ms)   |   x-eos ...\n");
			printf("══════════════════════════════════════════════════════════════════\n");
		}
		
		/** 
			update delta_G of pure phases as function of updated Gamma
		*/
		pp_min_function(				gv,
										z_b,
										PP_ref_db			);
			
		/**
			check driving force of PC when getting close to convergence
		*/
		v = clock();
		if (gv.BR_norm < gv.PC_check_val1 && gv.check_PC1 == 0 && pc_checked == 0){
			if (gv.verbose == 1){
				printf("\n Checking PC driving force 1\n");	
				printf("═════════════════════════════\n");	
					
			}
			gv = check_PC( 					z_b,						/** bulk rock constraint 				*/ 
											gv,							/** global variables (e.g. Gamma) 		*/
											PC_read,

											PP_ref_db,					/** pure phase database 				*/ 
											SS_ref_db,
											cp				); 					
			
			gv.check_PC1 		= 1;
			pc_checked 			= 1;				
		}
		/**
			check driving force of PC when getting close to convergence
		*/
		if (gv.BR_norm < gv.PC_check_val2 && gv.check_PC2 == 0 && pc_checked == 0){
			gv.PC_checked = 1;
			if (gv.verbose == 1){
				printf("\n Checking PC driving force 2\n");	
				printf("═════════════════════════════\n");	
					
			}
			gv = check_PC( 					z_b,						/** bulk rock constraint 				*/ 
											gv,							/** global variables (e.g. Gamma) 		*/
											PC_read,

											PP_ref_db,					/** pure phase database 				*/ 
											SS_ref_db,
											cp				); 					
			
			gv.check_PC2 		= 1;					
		}	

		/**
			Split phase if the current xeos is far away from the initial one 
		*/
		gv = split_cp(					gv, 						/** global variables (e.g. Gamma) 		*/
										SS_ref_db,					/** solution phase database 			*/	
										cp 					);		
		/**
			Local minimization of the solution phases
		*/

		ss_min_PGE(						gv, 						/** global variables (e.g. Gamma) 		*/
										PC_read,

										SS_objective,
										NLopt_opt,						
										z_b,						/** bulk-rock, pressure and temperature conditions */
										SS_ref_db,					/** solution phase database 			*/	
										cp 					);	
		v = clock() - v; 
		gv.tot_min_time += ((double)v)/CLOCKS_PER_SEC*1000;

		/**
			Merge instances of the same solution phase that are compositionnally close 
		*/
		gv = phase_merge_function(		z_b,							/** bulk rock constraint 				*/
										gv,								/** global variables (e.g. Gamma) 		*/

										PP_ref_db,						/** pure phase database 				*/
										SS_ref_db,						/** solution phase database 			*/ 
										cp					); 

		/**
			Actual Partitioning Gibbs Energy stage 
		/*/
		gv = PGE_inner_loop(			z_b,							/** bulk rock constraint 				*/ 
										splx_data,
										gv,								/** global variables (e.g. Gamma) 		*/

										PP_ref_db,						/** pure phase database 				*/ 
										SS_ref_db,						/** solution phase database 			*/
										cp					);

		if (gv.div == 1){
			gv.status = 4;
			iterate   = 0;
			break;
		}

		/* dump & print */
		if (gv.verbose == 1){
			/* Partitioning Gibbs Energy */
			PGE_print(					z_b,							/** bulk rock constraint 				*/ 
										gv,								/** global variables (e.g. Gamma) 		*/

										PP_ref_db,						/** pure phase database 				*/ 
										SS_ref_db,						/** solution phase database 			*/
										cp					); 
		}

		t = clock() - t; 
		if (gv.verbose == 1){
			printf("\n __ iteration duration: %+4f ms __\n\n\n",((double)t)/CLOCKS_PER_SEC*1000);
		}

		/* hard-clamp against PGE_mass_norm/Alg/gamma_norm/gibbs_ev/ite_time's actual
		   gv.it_f*2 capacity - see the matching guard in the LP loop above for why the
		   softer gv.it_f-based checks further down aren't sufficient on their own. */
		if (gv.global_ite >= gv.it_f*2){
			gv.global_ite = gv.it_f*2 - 1;
			gv.div        = 1;
			gv.status     = 4;
			iterate       = 0;
			if (gv.verbose == 1){
				printf(" >%d iterations (hard cap), not diverging but not converging\n\n", gv.it_f*2);
			}
			break;
		}

		/* check evolution of mass constraint residual */
		gv.PGE_mass_norm[gv.global_ite]  = gv.BR_norm;				/** save norm for the current global iteration */
		gv.Alg[gv.global_ite] 			 = 1;
		gv.ite_time[gv.global_ite] 		 = ((double)t)/CLOCKS_PER_SEC*1000;
		gv.global_ite 			  		+= 1;

		/*********************************************************/
		/**               CHECK MINIMIZATION STATUS              */
		/*********************************************************/

		/* checks for full convergence  						 */
		/* the second term checks if solution phase have been tested in case convergence is too fast */
		if (gv.BR_norm < gv.br_max_tol && gv.check_PC2 == 1){ gv.status = 0;	iterate = 0;}
		
		/* checks for dampened convergence  */
		if (gv.global_ite > gv.it_1 && gv.BR_norm < gv.br_max_tol*gv.ur_1){		if (gv.verbose == 1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_1, gv.ur_1);}; 	gv.status = 1; iterate = 0;}
		if (gv.global_ite > gv.it_2 && gv.BR_norm < gv.br_max_tol*gv.ur_2){		if (gv.verbose == 1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_2, gv.ur_2);}; 	gv.status = 2; iterate = 0;}
		if (gv.global_ite > gv.it_3 && gv.BR_norm < gv.br_max_tol*gv.ur_3){		if (gv.verbose == 1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_3, gv.ur_3);}; 	gv.status = 3; iterate = 0;}
		
		/* checks for not diverging but non converging cases  */
		if (gv.global_ite >= gv.it_f){  										if (gv.verbose == 1){printf(" >%d iterations, not diverging but not converging\n\n",gv.it_f);}	gv.div = 1; gv.status = 4; iterate = 0;}

		if ((log10(gv.BR_norm) > -1.5 && gv.global_ite > 64)	){	gv.div = 1;	iterate = 0;}
		if ((log10(gv.BR_norm) > -1.5 && gv.global_ite > 64)	){	gv.div = 1;	iterate = 0;}
		if ((log10(gv.BR_norm) > -2.5 && gv.global_ite > 128)	){	gv.div = 1;	iterate = 0;}
		if ((log10(gv.BR_norm) > -3.5 && gv.global_ite > 192)	){	gv.div = 1;	iterate = 0;}
		if (gv.gamma_norm[gv.global_ite-1] > 1e8				){	gv.div = 1;	iterate = 0;}
		
	}

	return gv;
};		