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
	extern void dgetrf( int* M, int* N, double* A, int* lda, int* ipiv, int* info);
	extern void dgetrs(char* T, int* N, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
#else
	#include <lapacke.h> 
#endif 


#include "MAGEMin.h"
#include "simplex_levelling.h"
#include "toolkit.h"
#include "gem_function.h"
#include "gss_function.h"
#include "dump_function.h"
#include "phase_update_function.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "objective_functions.h"
#include "NLopt_opt_function.h"

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
		if ((gv.pp_flags[i][2] == 1 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor < 5.0) || (gv.pp_flags[i][2] == 0 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor > 0.0)){
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

	gv.BR_norm    = norm_vector(	gv.mass_residual,
									z_b.nzEl_val				);

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

	// gsys = 0.0;
	// gref = 0.0;
	// for (int i = 0; i < z_b.nzEl_val; i++){ gref += z_b.bulk_rock[i]*gv.gam_tot[z_b.nzEl_array[i]]; }

	// for (int i = 0; i < z_b.nzEl_val; i++){ gsys += z_b.bulk_rock[i]*(gv.gam_tot[z_b.nzEl_array[i]] + gv.dGamma[i]); }
	// beta = 1.0;
	// if (gv.global_ite > 0){
	// 	if (gsys > gref){
	// 		beta = 0.9;
	// 		grel = gsys;
	// 		while (grel > (0.1 + gref)){
	// 			beta *= 0.5;
	// 			grel  = 0.0;
	// 			for (int i = 0; i < z_b.nzEl_val; i++){ grel += z_b.bulk_rock[i]*(gv.gam_tot[z_b.nzEl_array[i]] + gv.dGamma[i]*beta); }
	// 		}
	// 		printf("gsys: %+10f relax: %+10f grel: %+10f\n",gsys,beta,grel);
	// 	}
	// }
	// if (gv.alpha > beta){gv.alpha = beta;}

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
	while (PGEi < gv.inner_PGE_ite && delta_fc_norm > 1e-10){
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

/** 
  Partitioning Gibbs Energy function 
*/
global_variable PGE_inner_loop2(	bulk_info 			 z_b,
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
	while (PGEi < gv.inner_PGE_ite && delta_fc_norm > 1e-10){
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
global_variable init_LP(							bulk_info 	 		 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													
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
			
		/* if phase is a pure species */
		if (d->ph_id_A[i][0] == 1 ){
			gv.pp_flags[ph_id][1] 	= 1;
			gv.pp_flags[ph_id][2] 	= 0;
			gv.pp_n[ph_id]          = d->n_vec[i];
			gv.n_pp_phase		   += 1;
			gv.n_phase 			   += 1;
		}
		else {
			/* pure endmembers as solution phase */
			if (d->ph_id_A[i][0] == 2){
			em_id 					= d->ph_id_A[i][3];

			for (j = 0; j < SS_ref_db[ph_id].n_em; j++) {	
				SS_ref_db[ph_id].p[j] = gv.em2ss_shift;
			}
			SS_ref_db[ph_id].p[em_id] = 1.0 - gv.em2ss_shift*SS_ref_db[ph_id].n_em;
			
			SS_ref_db[ph_id] = P2X(			gv,
											SS_ref_db[ph_id],
											z_b,
											gv.SS_list[ph_id]		);
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
													SS_ref_db[ph_id], 
													z_b,
													gv.SS_list[ph_id] 			);

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
														SS_ref_db[ph_id], 
														z_b,
														gv.SS_list[ph_id] 		);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
														SS_ref_db[ph_id], 
														z_b, 
														gv.SS_list[ph_id]		);

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

/**
  Main LP routine
*/ 
global_variable LP(		bulk_info 			z_b,
						global_variable 	gv,

						obj_type 			*SS_objective,
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
	int nCheck  = 0;

	gv = init_LP(			z_b,
							splx_data,
							gv,
									
							PP_ref_db,
							SS_ref_db,
							cp	);

	while (iterate == 1){

		t = clock();

		if ((gv.gamma_norm[gv.global_ite-1] < 1.0 && nCheck < 3 && gv.global_ite > 1)){
			if (gv.verbose == 1){
				printf(" Checking PC for re-introduction:\n");
				printf(" ════════════════════════════════\n");
			}
			gv = check_PC( 				z_b,						/** bulk rock constraint 				*/ 
										gv,							/** global variables (e.g. Gamma) 		*/

										PP_ref_db,					/** pure phase database 				*/ 
										SS_ref_db,
										cp					); 
			if (gv.verbose == 1){
				printf("\n");
			}
			nCheck += 1;
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

										SS_objective,							
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

		/* check evolution of mass constraint residual */
		gv.PGE_mass_norm[gv.global_ite]  = gv.BR_norm;	/** save norm for the current global iteration */
		gv.Alg[gv.global_ite] 			 = 0;
		t 								 = clock() - t; 

		if (gv.verbose == 1){
			printf("\n __ iteration duration: %+4f ms __\n\n\n",((double)t)/CLOCKS_PER_SEC*1000);
		}
		gv.ite_time[gv.global_ite] 		 = ((double)t)/CLOCKS_PER_SEC*1000;
		gi += 1;

		if ((gv.gamma_norm[gv.global_ite-1] < 1e-4 || gi >= gv.max_LP_ite) && nCheck > 1){
			iterate = 0;
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

						obj_type 			*SS_objective,
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
			if (gv.verbose == 1){
				printf("\n Checking PC driving force 2\n");	
				printf("═════════════════════════════\n");	
					
			}
			gv = check_PC( 					z_b,						/** bulk rock constraint 				*/ 
											gv,							/** global variables (e.g. Gamma) 		*/

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

										SS_objective,							
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
	}

	return gv;
};		
