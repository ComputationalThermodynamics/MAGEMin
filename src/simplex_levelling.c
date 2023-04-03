/**
The purpose is to use G0 of pure-phases and endmembers 
of solid solutions to find out the potential stable       
phases. The final set of SS at equilibrium are contained  
in the estimated leveling assemblage.  

Levelling occurs in two stages:

- Stage 1: levelling using pure species only and a simplex approach
- Stage 2: levelling using progressive solution phase minimization and simplex approach
                     
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

#include "mpi.h"

#include "MAGEMin.h"
#include "simplex_levelling.h"
#include "gem_function.h"
#include "gss_function.h"
#include "NLopt_opt_function.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "dump_function.h"
#include "objective_functions.h"
#include "nlopt.h"
#include "toolkit.h"
#include "PGE_function.h"
#include "SS_xeos_PC_mp.h" 				//mp is first, it contains the structure definition
#include "SS_xeos_PC_ig.h"

/**
	associate the array of pointer with the right solution phase
*/
void SS_ig_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "bi")  == 0 ){
			SS_objective[iss]  = obj_ig_bi; 		}
		else if (strcmp( gv.SS_list[iss], "cd")  == 0){
			SS_objective[iss]  = obj_ig_cd; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_objective[iss]  = obj_ig_cpx; 		}
		else if (strcmp( gv.SS_list[iss], "ep")  == 0){
			SS_objective[iss]  = obj_ig_ep; 		}
		else if (strcmp( gv.SS_list[iss], "fl")  == 0){
			SS_objective[iss]  = obj_ig_fl; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_objective[iss]  = obj_ig_g; 		}
		else if (strcmp( gv.SS_list[iss], "hb")  == 0){
			SS_objective[iss]  = obj_ig_hb; 		}
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			SS_objective[iss]  = obj_ig_ilm; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			SS_objective[iss]  = obj_ig_liq; 		}
		else if (strcmp( gv.SS_list[iss], "mu")  == 0){
			SS_objective[iss]  = obj_ig_mu; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_objective[iss]  = obj_ig_ol; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_objective[iss]  = obj_ig_opx; 		}
		else if (strcmp( gv.SS_list[iss], "pl4T") == 0){
			SS_objective[iss]  = obj_ig_pl4T; 		}
		else if (strcmp( gv.SS_list[iss], "spn") == 0){
			SS_objective[iss]  = obj_ig_spn; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}

void SS_mp_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "liq")   == 0 ){
			SS_objective[iss]  = obj_mp_liq; 		}
		else if (strcmp( gv.SS_list[iss], "pl4tr") == 0){
			SS_objective[iss]  = obj_mp_pl4tr; 		}
		else if (strcmp( gv.SS_list[iss], "bi")    == 0){
			SS_objective[iss]  = obj_mp_bi; 		}
		else if (strcmp( gv.SS_list[iss], "g")     == 0){
			SS_objective[iss]  = obj_mp_g; 			}
		else if (strcmp( gv.SS_list[iss], "ep")    == 0){
			SS_objective[iss]  = obj_mp_ep; 		}
		else if (strcmp( gv.SS_list[iss], "ma")    == 0){
			SS_objective[iss]  = obj_mp_ma; 		}
		else if (strcmp( gv.SS_list[iss], "mu")    == 0){
			SS_objective[iss]  = obj_mp_mu; 		}
		else if (strcmp( gv.SS_list[iss], "opx")   == 0){
			SS_objective[iss]  = obj_mp_opx; 		}
		else if (strcmp( gv.SS_list[iss], "sa")    == 0){
			SS_objective[iss]  = obj_mp_sa; 		}
		else if (strcmp( gv.SS_list[iss], "cd")    == 0){
			SS_objective[iss]  = obj_mp_cd; 		}
		else if (strcmp( gv.SS_list[iss], "st")    == 0){
			SS_objective[iss]  = obj_mp_st; 		}
		else if (strcmp( gv.SS_list[iss], "chl")   == 0){
			SS_objective[iss]  = obj_mp_chl; 		}
		else if (strcmp( gv.SS_list[iss], "ctd")   == 0){
			SS_objective[iss]  = obj_mp_ctd; 		}
		else if (strcmp( gv.SS_list[iss], "sp")    == 0){
			SS_objective[iss]  = obj_mp_sp; 		}
		else if (strcmp( gv.SS_list[iss], "ilm")   == 0){
			SS_objective[iss]  = obj_mp_ilm; 		}
		else if (strcmp( gv.SS_list[iss], "mt")    == 0){
			SS_objective[iss]  = obj_mp_mt; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}
/**
  function to calculate delta G (position of potential phases with current G hyperplane*)
*/
void update_dG(	simplex_data *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;
	double F;
	double minF = 0.0;

	VecMatMul(	d->B1,
				d->A1,
				d->B,
				d->n_Ox	);
	
	d->dG_B = d->g0_B;
	for (int j = 0; j < d->n_Ox; j++){
		d->dG_B -= d->B1[j]*d->g0_A[j];
	}

	d->ph2swp = -1;	
	if (d->dG_B < d->dG_B_tol){	
		d->min_F  =  d->min_F_tol;												/** max value for F, tentative here because F can tend to +inf */
		for (int i = 0; i < d->n_Ox; i++){
			F = (d->n_vec[i])/d->B1[i];
			if (F < d->min_F && F > minF){
				d->min_F  = F;
				d->ph2swp = i;
			}
		}
	}
}




/**
  update Gamma using LAPACKE dgesv
*/	
void update_local_gamma(	double *A1, 
							double *g0_A, 
							double *gam, 
							int 	n			){
	int k;
	for (int i = 0; i < n; i++){
		gam[i] = 0.0;
		for (int j = 0; j < n; j++){
			k = i + j*n;
			gam[i] += g0_A[j]*A1[k];
		}
	}
};

/**
  update global Gamma 
*/	
void update_global_gamma( 				bulk_info 			z_b,
										simplex_data 	   *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

	/** update gam_tot using solution phase levelling gamma */
	for (int i = 0; i < d->n_Ox; i++){
		d->gamma_delta[z_b.nzEl_array[i]] 	= d->gamma_ss[i] - d->gamma_tot[z_b.nzEl_array[i]];
		d->gamma_tot[z_b.nzEl_array[i]]     = d->gamma_ss[i];
	}	
	
};

/**
  update global Gamma 
*/	
void update_global_gamma_LU( 				bulk_info 			z_b,
											simplex_data 	   *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

	int i,j,k,l;

	/* LAPACKE memory allocation */
	int 	nrhs   = 1;													/** number of rhs columns, 1 is vector*/
	int 	lda    = d->n_Ox;											/** leading dimesion of A*/
	int 	ldb    = 1;													/** leading dimension of b*/
	int 	ipiv[d->n_Ox];												/** pivot indices*/
	int 	info;														/** get info from lapacke function*/	

	for (i = 0; i < d->n_Ox;i++){
		d->gamma_ss[i] = d->g0_A[i];
	}
	for (i = 0; i < d->n_Ox;i++){
		for (j = 0; j < d->n_Ox;j++){
			k = i + j*d->n_Ox;
			l = j + i*d->n_Ox;
			d->Alu[k] = d->A[l];
		}
	}

	/**
		call lapacke to solve system of linear equation using LU 
	*/
#if __APPLE__
	lda    = d->n_Ox;											/** leading dimension of A*/
	ldb    = d->n_Ox;

	// Factorisation
	dgetrf(&d->n_Ox, &d->n_Ox, d->Alu, &lda, &ipiv, &info);

	char T = 'T';
	dgetrs(						&T,
								&d->n_Ox, 
								&nrhs, 
								d->Alu, 
								&lda, 
								ipiv, 
								d->gamma_ss, 
								&ldb,
								&info				);

#else	
	info = LAPACKE_dgesv(		LAPACK_ROW_MAJOR, 
								d->n_Ox, 
								nrhs, 
								d->Alu, 
								lda, 
								ipiv, 
								d->gamma_ss, 
								ldb					);
#endif
	/** update gam_tot using solution phase levelling gamma */
	for (int i = 0; i < d->n_Ox; i++){
		d->gamma_delta[z_b.nzEl_array[i]] 	= d->gamma_ss[i] - d->gamma_tot[z_b.nzEl_array[i]];
		d->gamma_tot[z_b.nzEl_array[i]]     = d->gamma_ss[i];
	}	

};

/**
  function to swp pure phases
*/	
void swap_pure_phases(				bulk_info 	 		 z_b,
									simplex_data 		*splx_data,
									global_variable 	 gv,
									
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;
	int     k;
	
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][0] == 1){
			
			d->g0_B 		= PP_ref_db[i].gbase*PP_ref_db[i].factor;
			d->ph_id_B[0]  	= 1;															/** added phase is a pure species */
			d->ph_id_B[1]	= i;															/** save pure species index */
		
			/* retrieve the composition in the right (reduced) chemical space */
			for (int j = 0; j < z_b.nzEl_val; j++){
				d->B[j] = PP_ref_db[i].Comp[z_b.nzEl_array[j]]*PP_ref_db[i].factor; 
			}
			
			/** update deltaG with respect to G hyperplane */
			update_dG(splx_data);
			
			/** swap phase */
			if (d->ph2swp != -1){															/** if the phase can be added */
				d->swp 						   = 1;
				d->n_swp 					  += 1;
				d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
				d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
				d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
				d->g0_A[d->ph2swp] 	     = d->g0_B;
				
				for (int j = 0; j < d->n_Ox; j++){				
					k = d->ph2swp + j*d->n_Ox;
					d->A[k] = d->B[j];
				}
				
				for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

				/** inverse guessed assemblage stoechiometry matrix */
				inverseMatrix(	gv.ipiv,
								d->A1,
								d->n_Ox,
								gv.work,
								gv.lwork		);

				/** update phase fractions */
				MatVecMul(		d->A1,
								z_b.bulk_rock_cat,
								d->n_vec,
								d->n_Ox		);	
			}
		}
	}
	
}

/**
  function to swp pure endmembers
*/	
void swap_pure_endmembers(				bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){

	simplex_data *d  = (simplex_data *) splx_data;

	int     k;
	double factor;

	for (int i = 0; i < gv.len_ss; i++){												/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){														/** if SS is not filtered out then continue */

			for (int l = 0; l < SS_ref_db[i].n_em; l++){	
				/** if bulk-rock satisfy the compositions of endmembers, retrieve their informations */
				if (SS_ref_db[i].z_em[l] == 1.0){ 
					
					/* update normalizing factor for solution models than need it */
					factor 	= z_b.fbc/SS_ref_db[i].ape[l];	

					d->g0_B		  = SS_ref_db[i].gbase[l]*factor;
					d->ph_id_B[0] = 2;														/** added phase is a pure species */
					d->ph_id_B[1] = i;														/** save pure species index */
					d->ph_id_B[2] = 0;														/** save used initial guess */	
					
					/** retrieve the composition in the right (reduced) chemical space */
					for (int j = 0; j < z_b.nzEl_val; j++){
						d->B[j] 	 = SS_ref_db[i].Comp[l][z_b.nzEl_array[j]]*factor; 
					}

					/** update deltaG with respect to G hyperplane */
					update_dG(splx_data);
									
					/** swap phase */
					if (d->ph2swp != -1){													/** if the phase can be added */
						d->swp 						   = 1;
						d->n_swp 					  += 1;
						d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
						d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
						d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
						d->ph_id_A[d->ph2swp][3] = l;	
						d->g0_A[d->ph2swp] 	   	 = d->g0_B;
						
						for (int j = 0; j < d->n_Ox; j++){				
							int k = d->ph2swp + j*d->n_Ox;
							d->A[k] = d->B[j];
						}
						for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

						/** inverse guessed assemblage stoechiometry matrix */
						inverseMatrix(	gv.ipiv,
										d->A1,
										d->n_Ox,
										gv.work,
										gv.lwork	);

						/** update phase fractions */
						MatVecMul(		d->A1,
										z_b.bulk_rock_cat,
										d->n_vec,
										d->n_Ox		);	
					}
				}
			}
		}
	}
	
}

/**
  function to swp solution phase pseudocompounds
*/	
void swap_pseudocompounds(				bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k, max_n_pc;
	
	for (int i = 0; i < gv.len_ss; i++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){

			max_n_pc = SS_ref_db[i].tot_pc;

			for (int l = 0; l < max_n_pc; l++){

				d->g0_B		  = SS_ref_db[i].G_pc[l];
				d->ph_id_B[0] = 3;														/** added phase is a pure species */
				d->ph_id_B[1] = i;														/** save solution phase index */
				d->ph_id_B[2] = 0;														/** save pseudocompound index */
			
				/* retrieve the composition in the right (reduced) chemical space */
				for (int j = 0; j < z_b.nzEl_val; j++){
					d->B[j] = SS_ref_db[i].comp_pc[l][z_b.nzEl_array[j]]; 
				}

				/** update deltaG with respect to G hyperplane */
				update_dG(splx_data);
				
				/** save updated driving force */
				SS_ref_db[i].DF_pc[l] = d->dG_B;
				
				/** swap phase */
				if (d->ph2swp != -1){													/** if the phase can be added */
					d->swp 				 	 = 1;
					d->n_swp 				+= 1;
					d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
					d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
					d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
					d->ph_id_A[d->ph2swp][3] = l;									/** save pseudocompound number */
					d->g0_A[d->ph2swp] 	     = d->g0_B;
					
					for (int j = 0; j < d->n_Ox; j++){				
						int k = d->ph2swp + j*d->n_Ox;
						d->A[k] = d->B[j];
					}
					
					for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

					/** inverse guessed assemblage stoechiometry matrix */
					inverseMatrix(	gv.ipiv,
									d->A1,
									d->n_Ox,
									gv.work,
									gv.lwork	);

					/** update phase fractions */
					MatVecMul(		d->A1,
									z_b.bulk_rock_cat,
									d->n_vec,
									d->n_Ox		);
				}
			}
		}
	}
}

/**
  function to swp solution phase pseudocompounds during PGE iterations
*/	
void swap_PGE_pseudocompounds(			bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k, max_n_Ppc;
	
	for (int i = 0; i < gv.len_ss; i++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){

			max_n_Ppc = SS_ref_db[i].tot_Ppc;

			for (int l = 0; l < max_n_Ppc; l++){

				d->g0_B		  = SS_ref_db[i].G_Ppc[l];
				d->ph_id_B[0] = 3;														/** added phase is a pure species */
				d->ph_id_B[1] = i;														/** save solution phase index */
				d->ph_id_B[2] = 0;														/** save pseudocompound index */
			
				/* retrieve the composition in the right (reduced) chemical space */
				for (int j = 0; j < z_b.nzEl_val; j++){
					d->B[j] = SS_ref_db[i].comp_Ppc[l][z_b.nzEl_array[j]]; 
				}

				/** update deltaG with respect to G hyperplane */
				update_dG(splx_data);
				
				/** save updated driving force */
				SS_ref_db[i].DF_Ppc[l] = d->dG_B;
				
				/** swap phase */
				if (d->ph2swp != -1){													/** if the phase can be added */
					d->swp 				 	 = 1;
					d->n_swp 				+= 1;
					d->ph_id_A[d->ph2swp][0] = d->ph_id_B[0];
					d->ph_id_A[d->ph2swp][1] = d->ph_id_B[1];
					d->ph_id_A[d->ph2swp][2] = d->ph_id_B[2];
					d->ph_id_A[d->ph2swp][3] = l;										/** save pseudocompound number */
					d->g0_A[d->ph2swp] 	     = d->g0_B;
					d->stage[d->ph2swp] 	 = 1;										/** just to indicate that the phase belongs to stage 2 of LP */
					
					for (int j = 0; j < d->n_Ox; j++){				
						int k = d->ph2swp + j*d->n_Ox;
						d->A[k] = d->B[j];
					}
					
					for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

					/** inverse guessed assemblage stoechiometry matrix */
					inverseMatrix(	gv.ipiv,
									d->A1,
									d->n_Ox,
									gv.work,
									gv.lwork	);

					/** update phase fractions */
					MatVecMul(		d->A1,
									z_b.bulk_rock_cat,
									d->n_vec,
									d->n_Ox		);
				}
			}
		}
	}

}



/**
  function to run simplex linear programming 
*/	
void fill_simplex_arrays_A(				bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){

	simplex_data *d  = (simplex_data *) splx_data;
	/* fill reference assemblage */

	for (int k = 0; k < z_b.nzEl_val; k++) {
		d->g0_A[k]    		    = 1e10;								/** penalty G */
		d->ph_id_A[k][0]  	    = 0;									/** phase_id for penalty phase */
		d->A[k+k*z_b.nzEl_val]  = 1.0;								/** eye matrix for stoichiometry */
		d->A1[k+k*z_b.nzEl_val] = 1.0;
		d->n_vec[k] 		    = z_b.bulk_rock[z_b.nzEl_array[k]];	/** initial phase fraction simply corresponds to the bulk rock composition Ax = br */
	}

}

/**
	print levelling informations 
*/
void print_levelling(		bulk_info 	 		 z_b,
							global_variable 	 gv,
							
							PP_ref 				*PP_ref_db,
							SS_ref 				*SS_ref_db
){
	int max_n_pc;
	
	printf("\n");
	printf("DISPLAY SWAP NUMBER FOR EACH PC\n");
	printf("-------------------------------\n");
	for (int i = 0; i < gv.len_ss; i++){

		int n_em = SS_ref_db[i].n_em;
		max_n_pc = SS_ref_db[i].tot_pc;	

		for (int l = 0; l < max_n_pc; l++){
			SS_ref_db[i].DF_pc[l] = SS_ref_db[i].G_pc[l];
			for (int j = 0; j < gv.len_ox; j++) {
				SS_ref_db[i].DF_pc[l] -= SS_ref_db[i].comp_pc[l][j]*gv.gam_tot[j];
			}
			
			// if (SS_ref_db[i].DF_pc[l] < 1.0){
				printf(" %4s %04d #stage %04d | ",gv.SS_list[i],l,SS_ref_db[i].info[l]);
				printf("DF: %+4f | ", SS_ref_db[i].DF_pc[l]);
				//for (int j = 0; j < SS_ref_db[i].n_xeos; j++){
					//printf(" %+4f",SS_ref_db[i].xeos_pc[l][j]);
				//}
				//for (int k = SS_ref_db[i].n_xeos; k < 10; k++){
					//printf(" %4s","-");
				//}
				//printf(" | ");

				// for (int j = 0; j < SS_ref_db[i].n_em; j++){
				// 	for (int k = 0; k < gv.len_ox; k++) {
				// 		SS_ref_db[i].mu_pc[l][j] -= SS_ref_db[i].comp_pc[l][k]*gv.gam_tot[k];
				// 	}
				// 	printf(" %+4f",SS_ref_db[i].mu_pc[l][j]);
				// }
				for (int k = SS_ref_db[i].n_em; k < 11; k++){
					printf(" %4s","-");
				}
				printf(" | ");
				for (int j = 0; j < SS_ref_db[i].n_em; j++){
					printf(" %+4f",SS_ref_db[i].p_pc[l][j]);
				}
				for (int k = SS_ref_db[i].n_em; k < 11; k++){
					printf(" %4s","-");
				}
				printf("\n");

			// }
		}
	}
}

/**
	Generate pseudocompounds
*/
void generate_pseudocompounds(	int 		 		 ss,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db,
								PC_ref 				*SS_pc_xeos,
								obj_type 			*SS_objective				){

	struct ss_pc get_ss_pv;							
	double sum, G, df, ape;
	int i,j,k,l,p,swp;
	int m_pc    = 0;								/** pseudocompound index to store negative DF vallues */
	
	/* change of base using Gamma */
	for (int k = 0; k < SS_ref_db[ss].n_em; k++) {
		SS_ref_db[ss].gb_lvl[k] = SS_ref_db[ss].gbase[k];
	}
						
	for (int k = 0; k < gv.n_SS_PC[ss]; k++){
		get_ss_pv = SS_pc_xeos[ss].ss_pc_xeos[k]; 
		
		/* TMP, not so elegant way to deal with cases were an oxide of the bulk rock composition = 0.0 */	
		for (int i = 0; i < SS_ref_db[ss].n_xeos; i++){
			if (get_ss_pv.xeos_pc[i] > SS_ref_db[ss].bounds_ref[i][1]){
				get_ss_pv.xeos_pc[i] = SS_ref_db[ss].bounds_ref[i][1];
			}
		}
		G 	= (*SS_objective[ss])(SS_ref_db[ss].n_xeos, get_ss_pv.xeos_pc, 	NULL, &SS_ref_db[ss]);

		/** store pseudocompound */
		if ( G < gv.max_G_pc ){

			/* get composition of solution phase */
			for (j = 0; j < gv.len_ox; j++){
				SS_ref_db[ss].ss_comp[j] = 0.0;
				for (i = 0; i < SS_ref_db[ss].n_em; i++){
					SS_ref_db[ss].ss_comp[j] += SS_ref_db[ss].Comp[i][j]*SS_ref_db[ss].p[i]*SS_ref_db[ss].z_em[i];
				} 
			}

			// if (SS_ref_db[ss].id_pc >= SS_ref_db[ss].n_pc){ SS_ref_db[ss].id_pc = 0; printf("MAXIMUM STORAGE SPACE FOR PC IS REACHED, INCREASED #PC_MAX\n");}
			m_pc = SS_ref_db[ss].id_pc;
			SS_ref_db[ss].info[m_pc]      = 0;
			SS_ref_db[ss].factor_pc[m_pc] = SS_ref_db[ss].factor;
			SS_ref_db[ss].DF_pc[m_pc]     = G;
			
			/* get pseudocompound composition */
			for ( j = 0; j < gv.len_ox; j++){				
				SS_ref_db[ss].comp_pc[m_pc][j] = SS_ref_db[ss].ss_comp[j]*SS_ref_db[ss].factor;	/** composition */
			}
			for ( j = 0; j < SS_ref_db[ss].n_em; j++){												/** save coordinates */
				SS_ref_db[ss].p_pc[m_pc][j]  = SS_ref_db[ss].p[j];												
				// SS_ref_db[ss].mu_pc[m_pc][j] = SS_ref_db[ss].mu[j]*SS_ref_db[ss].z_em[j];										
			}
			/* save xeos */
			for ( j = 0; j < SS_ref_db[ss].n_xeos; j++){		
				SS_ref_db[ss].xeos_pc[m_pc][j] = get_ss_pv.xeos_pc[j];							/** compositional variables */
			}	
			
			SS_ref_db[ss].G_pc[m_pc] = G;
			
			/* add increment to the number of considered phases */
			SS_ref_db[ss].tot_pc += 1;
			SS_ref_db[ss].id_pc  += 1;

		}
	}
}


/** 
  Get bounds for pseudocompounds compositional bounds after levelling
*/
void reduce_ss_list( SS_ref 			*SS_ref_db,
					 global_variable 	 gv				){
						 
	int max_n_pc;			 
	int phase_on;
	for (int iss = 0; iss < gv.len_ss; iss++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[iss].ss_flags[0] == 1){
			phase_on = 0;
			
			max_n_pc = SS_ref_db[iss].tot_pc;
			
			for (int l = 0; l < max_n_pc; l++){
				/* if the driving force of the pseudocompound is lower to the filter then consider it */
				if (SS_ref_db[iss].DF_pc[l]*SS_ref_db[iss].factor_pc[l] < gv.bnd_filter_pc){
					phase_on = 1;
				}	
			}
			/* if no generated pseudocompound are close to the hyperplane then turn off the solution phase */
			if (phase_on == 0){
				if (gv.verbose > 0){
					printf("  -> deleted = %s\n",gv.SS_list[iss]);
				}
				SS_ref_db[iss].ss_flags[0] = 0;
				SS_ref_db[iss].ss_flags[1] = 0;
				SS_ref_db[iss].ss_flags[2] = 0;
				SS_ref_db[iss].ss_flags[3] = 1;
			}
		}
	}				 
}

global_variable update_global_info(		bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp,
										obj_type 			*SS_objective
){	
	simplex_data *d  = (simplex_data *) splx_data;

	/* copy gamma total to the global variables */
	for (int i = 0; i < gv.len_ox; i++){
		gv.gam_tot[i] = d->gamma_tot[i];
	}

	double distance;
	double min_distance;
	double mid_dG;
	
	int ph_id, npc, id, id_min_distance;
	int id_cp = 0;
	int pc_id;
	int em_id;
	int phase_on[gv.len_ss];
	int i, j, k, ii;
	int m_pc;
	
	/* initialiaze phase active in the considered assemblage */
	for (i = 0; i < gv.len_ss; i++){
		phase_on[i] = 0;
	}

	/** 
		get initial conditions for active phases
	*/
	for (i = 0; i < d->n_Ox; i++){
		ph_id 		= d->ph_id_A[i][1];
			
		/* if phase is a pure species */
		if (d->ph_id_A[i][0] == 1 ){
			gv.pp_flags[ph_id][1] 	= 1;
			gv.pp_flags[ph_id][2] 	= 0;
			gv.pp_n[ph_id]          = d->n_vec[i];
			gv.n_pp_phase		   += 1;
			gv.n_phase 			   += 1;
		}
		/* pure endmembers as solution phase */
		if (d->ph_id_A[i][0] == 2){
			phase_on[ph_id] 		= 1;
			em_id 					= d->ph_id_A[i][3];

			for (j = 0; j < SS_ref_db[ph_id].n_em; j++) {	
				SS_ref_db[ph_id].p[j] = gv.em2ss_shift;
			}
			SS_ref_db[ph_id].p[em_id] = 1.0 - gv.em2ss_shift*SS_ref_db[ph_id].n_em;
			
			SS_ref_db[ph_id] = P2X(			gv,
											SS_ref_db[ph_id],
											z_b,
											gv.SS_list[ph_id]		);

			SS_ref_db[ph_id] = PC_function(	gv,
											SS_ref_db[ph_id], 
											z_b,
											gv.SS_list[ph_id] 		);
											
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
			
			cp[id_cp].ss_n          = d->n_vec[i];			/* get initial phase fraction */
			
			for (ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
				cp[id_cp].p_em[ii]  = SS_ref_db[ph_id].p[ii];
			}
			for (ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
				cp[id_cp].dguess[ii]  = SS_ref_db[ph_id].iguess[ii];
				cp[id_cp].xeos[ii]    = SS_ref_db[ph_id].iguess[ii];
			}
	
			gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
			gv.n_solvi[ph_id] 	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;

			/** 
			   add PC to Ppc list 
			*/
			if (SS_ref_db[ph_id].id_Ppc >= SS_ref_db[ph_id].n_Ppc){ SS_ref_db[ph_id].id_Ppc = 0; printf("MAXIMUM STORAGE SPACE FOR PC IS REACHED, INCREASED #PC_MAX\n");}
			m_pc = SS_ref_db[ph_id].id_Ppc;

			SS_ref_db[ph_id].info_Ppc[m_pc]   = 2;
			SS_ref_db[ph_id].factor_Ppc[m_pc] = SS_ref_db[ph_id].factor;
			SS_ref_db[ph_id].DF_Ppc[m_pc]     = SS_ref_db[ph_id].df;
			
			/* get pseudocompound composition */
			for ( j = 0; j < gv.len_ox; j++){				
				SS_ref_db[ph_id].comp_Ppc[m_pc][j] = SS_ref_db[ph_id].ss_comp[j]*SS_ref_db[ph_id].factor;	/** composition */
			}
			for ( j = 0; j < SS_ref_db[ph_id].n_em; j++){												/** save coordinates */
				SS_ref_db[ph_id].p_Ppc[m_pc][j]    = SS_ref_db[ph_id].p[j];										
				SS_ref_db[ph_id].mu_Ppc[m_pc][j]   = SS_ref_db[ph_id].mu[j];										
			}
			/* save xeos */
			for ( j = 0; j < SS_ref_db[ph_id].n_xeos; j++){		
				SS_ref_db[ph_id].xeos_Ppc[m_pc][j] = SS_ref_db[ph_id].iguess[j];							/** compositional variables */
			}	
			SS_ref_db[ph_id].G_Ppc[m_pc] = SS_ref_db[ph_id].df;
			
			/* add increment to the number of considered phases */
			SS_ref_db[ph_id].tot_Ppc += 1;
			SS_ref_db[ph_id].id_Ppc  += 1;


		}
		
		/* solution phase */
		if (d->ph_id_A[i][0] == 3){
			pc_id 					= d->ph_id_A[i][3];
			phase_on[ph_id] 		= 1;
			
			strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */	
			
			cp[id_cp].split 		= 0;							
			cp[id_cp].id 			= ph_id;						/* get phase id */
			cp[id_cp].n_xeos		= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
			cp[id_cp].n_em			= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
			cp[id_cp].n_sf			= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
			
			cp[id_cp].df			= SS_ref_db[ph_id].DF_pc[pc_id];
			cp[id_cp].factor		= SS_ref_db[ph_id].factor_pc[pc_id];	
			
			cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
			cp[id_cp].ss_flags[1] 	= 1;
			cp[id_cp].ss_flags[2] 	= 0;
			
			cp[id_cp].ss_n          = d->n_vec[i];			/* get initial phase fraction */
			
			for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
				cp[id_cp].p_em[ii]  = SS_ref_db[ph_id].p_pc[pc_id][ii];
			}
			for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
				cp[id_cp].dguess[ii]  = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
				cp[id_cp].xeos[ii]    = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
			}

			gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
			gv.n_solvi[ph_id] 	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;


			/** 
			   add PC to Ppc list 
			*/
			if (SS_ref_db[ph_id].id_Ppc >= SS_ref_db[ph_id].n_Ppc){ SS_ref_db[ph_id].id_Ppc = 0; printf("MAXIMUM STORAGE SPACE FOR PC IS REACHED, INCREASED #PC_MAX\n");}
			m_pc = SS_ref_db[ph_id].id_Ppc;

			SS_ref_db[ph_id].info_Ppc[m_pc]   = SS_ref_db[ph_id].info[pc_id];
			SS_ref_db[ph_id].factor_Ppc[m_pc] = SS_ref_db[ph_id].factor_pc[pc_id];
			SS_ref_db[ph_id].DF_Ppc[m_pc]     = SS_ref_db[ph_id].DF_pc[pc_id];
			
			/* get pseudocompound composition */
			for ( j = 0; j < gv.len_ox; j++){				
				SS_ref_db[ph_id].comp_Ppc[m_pc][j] = SS_ref_db[ph_id].comp_pc[pc_id][j];	/** composition */
			}
			for ( j = 0; j < SS_ref_db[ph_id].n_em; j++){												/** save coordinates */
				SS_ref_db[ph_id].p_Ppc[m_pc][j]  = SS_ref_db[ph_id].p_pc[pc_id][j];										
				// SS_ref_db[ph_id].mu_Ppc[m_pc][j] = SS_ref_db[ph_id].mu_pc[pc_id][j];										
			}
			/* save xeos */
			for ( j = 0; j < SS_ref_db[ph_id].n_xeos; j++){		
				SS_ref_db[ph_id].xeos_Ppc[m_pc][j] = SS_ref_db[ph_id].xeos_pc[pc_id][j];							/** compositional variables */
			}	
			SS_ref_db[ph_id].G_Ppc[m_pc] = SS_ref_db[ph_id].G_pc[pc_id];
			
			/* add increment to the number of considered phases */
			SS_ref_db[ph_id].tot_Ppc += 1;
			SS_ref_db[ph_id].id_Ppc  += 1;

		}
	}

	/** 
		get initial conditions for phases not active but close to the G-hyperplane (not predicted to be stable during levelling)
	*/
	int max_n_pc;
	double pc_df;
	for (int i = 0; i < gv.len_ss; i++){
		ph_id = i;
		if (phase_on[i] == 0 && SS_ref_db[i].ss_flags[0] == 1){
			pc_id	 = -1;
			pc_df 	 = 1e6;
			max_n_pc = SS_ref_db[i].tot_pc;
			
			for (int l = 0; l < max_n_pc; l++){
				if (SS_ref_db[i].DF_pc[l]*SS_ref_db[i].factor_pc[l] < pc_df){
					pc_df   = SS_ref_db[i].DF_pc[l]*SS_ref_db[i].factor_pc[l];
					pc_id	= l;
				}	
			}
			
			if (pc_id != -1){
				strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */	
				
				cp[id_cp].id 			= ph_id;						/* get phaseid */
				cp[id_cp].n_xeos		= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
				cp[id_cp].n_em			= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
				cp[id_cp].n_sf			= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
					
				cp[id_cp].df			= SS_ref_db[ph_id].DF_pc[pc_id];
				cp[id_cp].factor		= SS_ref_db[ph_id].factor_pc[pc_id];	
				
				cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
				cp[id_cp].ss_flags[1] 	= 0;
				cp[id_cp].ss_flags[2] 	= 1;
				
				cp[id_cp].ss_n          = 0.0;							/* get initial phase fraction */
				
				for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
					cp[id_cp].p_em[ii]    = SS_ref_db[ph_id].p_pc[pc_id][ii];
				}
				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					cp[id_cp].dguess[ii]  = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
					cp[id_cp].xeos[ii]    = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
				}
				
				gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
				gv.n_solvi[ph_id] 	+= 1;
				gv.len_cp 			+= 1;
				id_cp 				+= 1;
				
				if (gv.len_cp == gv.max_n_cp){
					printf(" !! Maxmimum number of allowed phases under consideration reached !!\n    -> check your problem and potentially increase gv.max_n_cp\n");
					break;
				}
			}
		}
	}

	/**
		Routine to deactivate the liquid endmembers after levelling 
	*/
	char liq_tail[] = "L";
	for (int i = 0; i < gv.len_pp; i++){
		if ( EndsWithTail(gv.PP_list[i], liq_tail) == 1 ) {
			if (gv.pp_flags[i][0] == 1){
				if (gv.pp_flags[i][1] == 1){
					gv.pp_flags[i][0] = 0;
					gv.pp_flags[i][1] = 0;
					gv.pp_flags[i][2] = 0;
					gv.pp_flags[i][3] = 1;
					gv.n_phase       -= 1;
					gv.n_pp_phase    -= 1;
					gv.pp_n[i]        = 0.0;
				}
				else{
					gv.pp_flags[i][0] = 0;
					gv.pp_flags[i][1] = 0;
					gv.pp_flags[i][2] = 0;
					gv.pp_flags[i][3] = 1;
				}
			}
		}
	}
	
	if (gv.verbose == 1){
		printf("\n Initial guesses for compositional variables:\n");
		printf("═════════════════════════════════════════════\n");
		
		for (int id_cp = 0; id_cp < gv.len_cp; id_cp++){
			if (cp[id_cp].ss_flags[0] == 1){
				printf(" %5s [%+10f]->  ",cp[id_cp].name,cp[id_cp].df*cp[id_cp].factor);
				for (int i = 0; i < cp[id_cp].n_xeos; i++){
					printf(" %+10f",cp[id_cp].dguess[i]);
				}
				for (int k = cp[id_cp].n_xeos; k < 11; k++){
					printf(" %10s","-");
				}
				printf("\n");
			}
		}
		printf("\n");	
	}
	
	return gv;
}


/**
  function to run simplex linear programming with pseudocompounds only
*/	
void run_simplex_pseudocompounds(		bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k = 0;

	d->swp = 1;
	while (d->swp == 1){					/** as long as a phase can be added to the guessed assemblage, go on */
		k 		   += 1;
		d->swp     = 0;

		swap_pure_endmembers(				z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db	);	

		swap_pure_phases(					z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db	);	
							
		swap_pseudocompounds(				z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db	);				
	}
	if (gv.verbose == 1){
		printf("    (# iterations %d)",k);	
	}

}

/**
  function to run simplex linear programming with pseudocompounds
*/	
void run_simplex_levelling(				bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										obj_type			*SS_objective
){
	simplex_data *d  = (simplex_data *) splx_data;

	int i, k, iss;

	swap_pure_phases(					z_b,
										splx_data,
										gv,
										PP_ref_db,
										SS_ref_db				);	

	swap_pure_endmembers(				z_b,
										splx_data,
										gv,
										PP_ref_db,
										SS_ref_db				);	


	update_local_gamma(					d->A1,
										d->g0_A,
										d->gamma_ps,
										d->n_Ox					);


	/** update gam_tot using pure species levelling gamma */
	for (i = 0; i < d->n_Ox; i++){
		d->gamma_tot[z_b.nzEl_array[i]] = d->gamma_ps[i];
	}

	/** 
		Pseudocompound levelling with refinement 
	*/
	clock_t t; 
	double time_taken;
	t = clock();

	/** generate the pseudocompounds -> stored in the SS_ref_db structure */
	if (gv.verbose == 1){ printf(" Generate pseudocompounds:\n"); }
	
	PC_ref 			SS_pc_xeos[gv.len_ss];


	if (gv.EM_database == 0){
		for (iss = 0; iss < gv.len_ss; iss++){
			SS_mp_pc_init_function(			SS_pc_xeos, 
											iss,
											gv.SS_list[iss]				);
		}
	}
	else if (gv.EM_database == 2){
		for (iss = 0; iss < gv.len_ss; iss++){
			SS_ig_pc_init_function(			SS_pc_xeos, 
											iss,
											gv.SS_list[iss]				);
		}
	}

	for (iss = 0; iss < gv.len_ss; iss++){
		if (SS_ref_db[iss].ss_flags[0] == 1){

			generate_pseudocompounds(	iss,
										z_b,
										gv,
										SS_ref_db,
										SS_pc_xeos,
										SS_objective				);

			if (gv.verbose == 1){
				printf(" %4s -> %05d active PCs\n",gv.SS_list[iss],SS_ref_db[iss].tot_pc);
			}
		}

	}

	t = clock() - t; 
	time_taken  = ((double)t)/CLOCKS_PER_SEC; 
	if (gv.verbose == 1){ printf("\n [time to generate PC time (ms) %.8f]\n",time_taken*1000);	}
	t = clock();
	
	/** run linear programming with simplex approach */
	run_simplex_pseudocompounds(		z_b,
										splx_data,
										gv,
										PP_ref_db,
										SS_ref_db				);
				
	/* update gamma of SS */
	update_local_gamma(					d->A1,
										d->g0_A,
										d->gamma_ss,
										d->n_Ox			);
										
	t = clock() - t; 
	time_taken  = ((double)t)/CLOCKS_PER_SEC;
	if (gv.verbose == 1){ printf("\n [time to swap SS time (ms) %.8f]\n",time_taken*1000);	}
	
};


void run_localMinimization(				bulk_info 	 		 z_b,
										simplex_data 		*splx_data,
										
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										obj_type			*SS_objective
){
	simplex_data *d  = (simplex_data *) splx_data;

	int i, j, k, ss, l, p, swp;
										
	clock_t t,u; 
	double time_taken;
	t = clock();

	/** generate the pseudocompounds -> stored in the SS_ref_db structure */
	if (gv.verbose == 1){ printf(" Generate pseudocompounds:\n"); }
	
	PC_ref 			SS_pc_xeos[gv.len_ss];

	if (gv.EM_database == 0){
		for (ss = 0; ss < gv.len_ss; ss++){
			SS_mp_pc_init_function(			SS_pc_xeos, 
											ss,
											gv.SS_list[ss]				);
		}
	}
	else if (gv.EM_database == 2){
		for (ss = 0; ss < gv.len_ss; ss++){
			SS_ig_pc_init_function(			SS_pc_xeos, 
											ss,
											gv.SS_list[ss]				);
		}
	}

	ss = 6; // hb index of the solution phase to fully minimize
	// ss = 0; // spn index of the solution phase to fully minimize
	// ss = 3; // cpx index of the solution phase to fully minimize

	struct ss_pc get_ss_pv;		

	gv.gam_tot[0]  = -960.9655;	
	gv.gam_tot[1]  = -1768.2476;	
	gv.gam_tot[2]  = -788.4474;	
	gv.gam_tot[3]  = -678.9683;	
	gv.gam_tot[4]  = -355.2975;	
	gv.gam_tot[5]  = -914.9708;	
	gv.gam_tot[6]  = -839.9561;
	gv.gam_tot[7]  = -1008.3630;
	gv.gam_tot[8]  = -263.7269;
	gv.gam_tot[9]  = -1262.6087;
	gv.gam_tot[10] = -368.4674;

	SS_ref_db[ss].gbase[0]  = -13012.62073;	
	SS_ref_db[ss].gbase[1]  = -13235.27114;	
	SS_ref_db[ss].gbase[2]  = -13472.30496;	
	SS_ref_db[ss].gbase[3]  = -12644.70794;	
	SS_ref_db[ss].gbase[4]  = -12762.02635;	
	SS_ref_db[ss].gbase[5]  = -10496.70590;	
	SS_ref_db[ss].gbase[6]  = -11477.04324;
	SS_ref_db[ss].gbase[7]  = -11155.59746;
	SS_ref_db[ss].gbase[8]  = -11828.15800;
	SS_ref_db[ss].gbase[9]  = -13495.08535;
	SS_ref_db[ss].gbase[10] = -13063.17373;

	// SS_ref_db[ss].gbase[0]  = -2515.94540;	
	// SS_ref_db[ss].gbase[1]  = -2500.25887;	
	// SS_ref_db[ss].gbase[2]  = -2217.27620;	
	// SS_ref_db[ss].gbase[3]  = -2201.58966;	
	// SS_ref_db[ss].gbase[4]  = -1452.03460;	
	// SS_ref_db[ss].gbase[5]  = -1459.64806;	
	// SS_ref_db[ss].gbase[6]  = -2033.47165;
	// SS_ref_db[ss].gbase[7]  = -2445.48343;

	// SS_ref_db[ss].gbase[0]  = -3532.74915;	
	// SS_ref_db[ss].gbase[1]  = -2793.12846;	
	// SS_ref_db[ss].gbase[2]  = -3635.49886;	
	// SS_ref_db[ss].gbase[3]  = -3384.95041;	
	// SS_ref_db[ss].gbase[4]  = -3250.67812;	
	// SS_ref_db[ss].gbase[5]  = -3606.43710;	
	// SS_ref_db[ss].gbase[6]  = -3345.42582;
	// SS_ref_db[ss].gbase[7]  = -3408.36774;
	// SS_ref_db[ss].gbase[8]  = -3105.14810;
	// SS_ref_db[ss].gbase[9]  = -3360.74459;


	// gv.gam_tot[0]  = -1011.909631;	
	// gv.gam_tot[1]  = -1829.092564;	
	// gv.gam_tot[2]  = -819.264126;	
	// gv.gam_tot[3]  = -695.467358;	
	// gv.gam_tot[4]  = -412.948568;	
	// gv.gam_tot[5]  = -971.890270;	
	// gv.gam_tot[6]  = -876.544354;
	// gv.gam_tot[7]  = -1073.640927;
	// gv.gam_tot[8]  = -276.590707;
	// gv.gam_tot[9]  = -1380.299631;
	// gv.gam_tot[10] = 0.0;


	/** rotate gbase with respect to the G-hyperplane (change of base) */
	for (int k = 0; k < SS_ref_db[ss].n_em; k++) {
		SS_ref_db[ss].gb_lvl[k] = SS_ref_db[ss].gbase[k];
		for (int j = 0; j < gv.len_ox; j++) {
			SS_ref_db[ss].gb_lvl[k] -= SS_ref_db[ss].Comp[k][j]*gv.gam_tot[j];
		}
	}	

	printf("minG = [");				
	for (int k = 0; k < gv.n_SS_PC[ss]; k++){
		u = clock();
		get_ss_pv = SS_pc_xeos[ss].ss_pc_xeos[k]; 
		
		for (int i = 0; i < SS_ref_db[ss].n_xeos; i++){
			SS_ref_db[ss].iguess[i] = get_ss_pv.xeos_pc[i];
		}

		SS_ref_db[ss] = 	NLopt_opt_function(		gv, 
													SS_ref_db[ss], 
													ss					);

		u = clock() - u; 
		time_taken  = ((double)u)/CLOCKS_PER_SEC; 
		printf(" %.14f", SS_ref_db[ss].df);
	}
	printf("]\n");	

	printf("tms = [");				
	for (int k = 0; k < gv.n_SS_PC[ss]; k++){
		u = clock();
		get_ss_pv = SS_pc_xeos[ss].ss_pc_xeos[k]; 
		
		for (int i = 0; i < SS_ref_db[ss].n_xeos; i++){
			SS_ref_db[ss].iguess[i] = get_ss_pv.xeos_pc[i];
		}

		SS_ref_db[ss] = 	NLopt_opt_function(		gv, 
													SS_ref_db[ss], 
													ss					);

		u = clock() - u; 
		time_taken  = ((double)u)/CLOCKS_PER_SEC; 
		printf(" %.8f", time_taken );
	}
	printf("]\n");	

	t = clock() - t; 
	time_taken  = ((double)t)/CLOCKS_PER_SEC; 
	if (gv.verbose == 1){ printf("\n [time to local minimization PC time (ms) %.8f]\n",time_taken*1000);	}
	
};


/**
  function to deallocte memory of simplex linear programming (A)
*/	
void destroy_simplex_A(		simplex_data *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

    for (int i = 0; i < d->n_Ox; i++){
		free(d->ph_id_A[i]);
	}
	free(d->ph_id_A);

	free(d->A);
	free(d->A1);
	free(d->Alu);

	free(d->pivot);

	free(d->g0_A);
	free(d->dG_A);
	
	free(d->n_vec);

	free(d->gamma_ps);
	free(d->gamma_ss);
	free(d->gamma_tot);
	free(d->gamma_delta);

	free(d->stage);
};



/**
  function to deallocte memory of simplex linear programming (B)
*/	
void destroy_simplex_B(
	simplex_data *splx_data
){
	simplex_data *d  = (simplex_data *) splx_data;

	free(d->ph_id_B);
	free(d->B);
	free(d->B1);
};

/**
  Levelling function.
  Memory allocation and initialization is divided in part A and B.
  Part A: Reference assemblage against which candidates are tested
  Part B: Candidate phases
*/	
global_variable run_levelling_function(		bulk_info 	 z_b,
											global_variable 	 gv,

											obj_type			*SS_objective,	
											simplex_data		*splx_data,
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp
){
	simplex_data *d  = (simplex_data *) splx_data;

	clock_t t; 
	double time_taken;
	t = clock();

	/** fill matrices */
	fill_simplex_arrays_A(					z_b,
										    splx_data,
											gv,
											PP_ref_db,
											SS_ref_db		);

	/** run linear programming with simplex approach */
	run_simplex_levelling(					z_b,
										    splx_data,
											gv,
											PP_ref_db,
											SS_ref_db,
											SS_objective	);	

	/** run local minimization tests */
	// run_localMinimization(					z_b,
	// 									    splx_data,
	// 										gv,
	// 										PP_ref_db,
	// 										SS_ref_db,
	// 										SS_objective	);	
			
	/* update global variable gamma */
	update_global_gamma_LU(					z_b,
											splx_data		);	
	
										
	/* remove solution from consideration when min driving force is > gv.bnd_filter_pc */
	reduce_ss_list( 						SS_ref_db, 
											gv 				);
	
	
	/* function to send back the updated initial guess, and phases flags */
	gv = update_global_info(				z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db,
											cp,
											SS_objective	);

	if (gv.verbose == 1){
		printf("\nGet initial guess (Gamma and phase fractions) \n");
		printf("══════════════════════════════════════════════\n\n");
		printf("    P: %+10f T: %+10f\n",z_b.P,z_b.T);
		printf(" [----------------------------------------]\n");
		printf(" [  Ph  |   Ph PROP  |   g0_Ph    |  ix   ]\n");
		printf(" [----------------------------------------]\n");

		for (int i = 0; i < d->n_Ox; i++){
			if (d->ph_id_A[i][0] == 0){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]", "F.OX", d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 1){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]", gv.PP_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 2){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]\n", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
			}
			if (d->ph_id_A[i][0] == 3){
				printf(" ['%5s' %+10f  %+12.4f  %5d ]", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
					printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_pc[d->ph_id_A[i][3]][ii] );
				}
				printf("\n");
			}
		}
		printf("\n");
		for (int i = 0; i < d->n_Ox; i++){
			printf(" %g", d->gamma_tot[z_b.nzEl_array[i]]);
		}
		printf("\n");
		printf(" [----------------------------------------]\n");
		printf(" [  OXIDE      GAMMA_EM        GAMMA_PC   ]\n");
		printf(" [----------------------------------------]\n");
		for (int i = 0; i < d->n_Ox; i++){
			printf(" [ %5s %+15f %+15f  ]\n", gv.ox[z_b.nzEl_array[i]], d->gamma_ps[i], d->gamma_tot[z_b.nzEl_array[i]]);
		}
		printf(" [----------------------------------------]\n");
		printf(" [            %4d swaps                  ]\n", d->n_swp);
		printf(" [----------------------------------------]\n");
		
		printf("\n [----------------------------------------]\n");
		printf(" [           ACTIVE PHASES                ]\n");
		printf(" [----------------------------------------]\n");
		for (int i = 0; i < gv.len_ss; i++){
			if (SS_ref_db[i].ss_flags[0] == 1){
				printf(" [                 %5s                  ]\n",gv.SS_list[i]);
			}
		}
		printf(" [----------------------------------------]\n");
		printf(" [           UNACTIVE PHASES              ]\n");
		printf(" [----------------------------------------]\n");
		for (int i = 0; i < gv.len_ss; i++){
			if (SS_ref_db[i].ss_flags[0] == 0){
				printf(" [                 %5s                  ]\n",gv.SS_list[i]);
			}
		}
	}

	t 			= clock() - t; 
	time_taken  = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	gv.LVL_time = time_taken*1000;	
	
	return gv;
};		

/**
  main levelling routine
*/ 
global_variable Levelling(	bulk_info 	z_b,
							global_variable 	gv,

							obj_type 		   *SS_objective,
							simplex_data	   *splx_data,
							PP_ref 			   *PP_ref_db,
							SS_ref 			   *SS_ref_db,
							csd_phase_set  	   *cp
){

	if (gv.verbose == 1){
		printf("\nLevelling (endmembers & solution phase)\n");
		printf("════════════════════════════════════════\n");
	}

	/* pseudosection function to get starting guess */
	gv = run_levelling_function(	z_b,												/** bulk rock informations    */
									gv,													/** global variables (e.g. Gamma) */

									SS_objective,
								    splx_data,
									PP_ref_db,											/** pure phase database */
									SS_ref_db,											/** solution phase database */
									cp				);
	if (gv.verbose == 1){
		printf(" [    Levelling time  %+12f ms     ]\n",gv.LVL_time);
		printf(" [----------------------------------------]\n\n\n");
	}

	/* copy gamma total to the global variables */
	// for (int i = 0; i < gv.len_ox; i++){
	// 	gv.gam_tot[i] -= 1.0;
	// }

	return gv;
};