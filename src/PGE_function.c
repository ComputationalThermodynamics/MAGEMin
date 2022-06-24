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
#include <lapacke.h> 

#include "MAGEMin.h"
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
	printf("ON | phase |  Fraction |  delta_G   |  factor   |   sum_xi   |    Pi - Xi...\n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1 ){

			printf(" %d | %4s | %+10f | %+10f | %+10f | %+10f | ",cp[i].ss_flags[1],cp[i].name,cp[i].ss_n,cp[i].df,cp[i].factor,cp[i].sum_xi);

			for (int k = 0; k < cp[i].n_em; k++) {
				printf(" %+12f",(cp[i].p_em[k]-cp[i].xi_em[k]*cp[i].p_em[k])*SS_ref_db[cp[i].id].z_em[k]);
			}
			for (int k = cp[i].n_em; k < 12; k++){
				printf(" %12s","-");
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
				printf(" %+10f",cp[i].xeos[k]);
			}
			for (int k = cp[i].n_em; k < 12; k++){
				printf(" %10s","-");
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
	printf("OFF| phase |  Fraction |  delta_G   |  factor   |   sum_xi   |    Pi - Xi...\n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].ss_flags[2] == 1){

			printf(" %d | %4s | %+10f | %+10f | %+10f | %+10f | ",cp[i].ss_flags[1],cp[i].name,cp[i].ss_n,cp[i].df,cp[i].factor,cp[i].sum_xi);

			for (int k = 0; k < cp[i].n_em; k++) {
				printf(" %+10f",(cp[i].p_em[k]-cp[i].xi_em[k]*cp[i].p_em[k])*SS_ref_db[cp[i].id].z_em[k]);
			}
			for (int k = cp[i].n_em; k < 12; k++){
				printf(" %10s","-");
			}
			printf("\n");
						
		}
	}
	printf("\n");
	printf("OFF| P. phase |  Fraction  |  delta_G  (< 5.0) | \n");
	printf("═══════════════════════════════════════════════════════════════════════════════════════\n");
	for (int i = 0; i < gv.len_pp; i++){ 
		if (gv.pp_flags[i][2] == 1 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor < 50.0){
			printf(" %d | %4s     | %+10f | %+10f | \n",0,gv.PP_list[i],gv.pp_n[i],PP_ref_db[i].gb_lvl*PP_ref_db[i].factor);
		}
	}
	printf("\n");
	printf("  Oxide   | Chemical potential\n");
	printf("═══════════════════════════════\n");
	for (int i = 0; i < z_b.nzEl_val; i++){		
		printf(" %6s   |  %+12.4f \n",gv.ox[z_b.nzEl_array[i]],gv.gam_tot[z_b.nzEl_array[i]]);
	}

	printf("\n");
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

			/** calculate residual as function xi fraction and not endmember fractions from x-eos */
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
	else if(gv.LP == 0 && gv.PGE == 1){
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
	else {
		printf(" Wrong combination of LP and PGE, check called sub-routines...\n");
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
		if (cp[i].ss_flags[0] == 1 && (cp[i].ss_flags[1] == 1 || cp[i].ss_flags[2] == 1)){
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
global_variable PGE_update_pi(		bulk_info 	z_b,
									global_variable  	gv,

									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db,
									csd_phase_set  		*cp
){
	int ss,ph,i,j,k,x,v;	
	
	for (ph = 0; ph < gv.len_cp; ph++){

		if (cp[ph].ss_flags[1] == 1 && SS_ref_db[cp[ph].id].CstFactor == 0){
			ss = cp[ph].id;

			double dp[cp[ph].n_em];
			for (k = 0; k < cp[ph].n_em; k++){
				dp[k] = (cp[ph].p_em[k]-cp[ph].xi_em[k]*cp[ph].p_em[k])*SS_ref_db[ss].z_em[k];
			}

			// if (norm_vector(dp ,cp[ph].n_em) > 1e-4){
				for (k = 0; k < cp[ph].n_em; k++){
					SS_ref_db[ss].p[k] = cp[ph].p_em[k]*cp[ph].xi_em[k];
				}
				
				SS_ref_db[ss] = P2X(			gv,
												SS_ref_db[ss],
												z_b,
												gv.SS_list[ss]				);
					
				for (j = 0; j < cp[ph].n_xeos; j++){
					SS_ref_db[ss].iguess[j] = cp[ph].xeos[j]*gv.xi_em_cor + SS_ref_db[ss].iguess[j]*(1.0 - gv.xi_em_cor); 
				}												
																		
				SS_ref_db[ss] = PC_function(	gv,
												SS_ref_db[ss], 
												z_b,
												gv.SS_list[ss] 				);
														
				if (SS_ref_db[ss].sf_ok == 1){						
					for (j = 0; j < cp[ph].n_xeos; j++){
						cp[ph].xeos[j] =  SS_ref_db[ss].iguess[j]; 
					}

				}
			// }
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
		if (cp[i].ss_flags[0] == 1 && (cp[i].ss_flags[1] == 1 || cp[i].ss_flags[2] == 1)){
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
	check PC driving force and add phase if below hyperplane
*/
global_variable check_EM(					bulk_info 	 z_b,
											global_variable  	 gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){
	double df, factor;
	for (int i = 0; i < gv.len_ss; i++){												/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){														/** if SS is not filtered out then continue */

			for (int l = 0; l < SS_ref_db[i].n_em; l++){	
				/** if bulk-rock satisfy the compositions of endmembers, retrieve their informations */
				if (SS_ref_db[i].z_em[l] == 1.0){
					
					/* update normalizing factor for solution models than need it */
					factor 	= z_b.fbc/SS_ref_db[i].ape[l];	

					df = SS_ref_db[i].gbase[l];
					for (int j = 0; j < gv.len_ox; j++) {
						df -= SS_ref_db[i].Comp[l][j]*gv.gam_tot[j];
					}
					
					if (df*factor < 0.0){
						printf("WARN: %4s %d %+10f\n",gv.SS_list[i],l,df*factor);
					}
				}
			}
		}
	}
	
	return gv;
}

/**
	check PC driving force and add phase if below hyperplane
*/
global_variable check_PC(					bulk_info 	 z_b,
											global_variable  	 gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){
	double 	min_df, xeos_dist, norm;
	int 	max_n_pc, phase_add, id_cp, dist, ph;
	int 	i,j,k,l,c,m;

	int 	n_candidate = 8;
	int     pc_candidate[n_candidate];
	int     pc_added[n_candidate];
	double  df_candidate[n_candidate];
	int     id_c;



	for (i = 0; i < gv.len_ss; i++){
		min_df    =  1e6;					// high starting value as it is expected to go down
		phase_add =  0;
		id_c 	  =  0;

		for (k = 0; k < n_candidate; k++){
			pc_candidate[k] = -1;
			pc_added[k] 	= -1;
			df_candidate[k] = 0.0;
		}

		if (SS_ref_db[i].ss_flags[0] == 1  && gv.verifyPC[i] == 1){
			
			max_n_pc  = ((SS_ref_db[i].tot_pc >= SS_ref_db[i].n_pc) ? (SS_ref_db[i].n_pc) : (SS_ref_db[i].tot_pc));
			
			for (l = 0; l < max_n_pc; l++){

				dist =  1;
				if (gv.n_solvi[i] > 0){

					for (k = 0; k < gv.n_solvi[i]; k++){  /* go through the upper triangle of the matrix (avoiding diagonal)*/
						ph = SS_ref_db[i].solvus_id[k];

						xeos_dist = euclidean_distance(cp[ph].xeos, SS_ref_db[i].xeos_pc[l], SS_ref_db[i].n_xeos);
						if (xeos_dist < gv.PC_min_dist*gv.SS_PC_stp[i]*sqrt((double)SS_ref_db[i].n_xeos) ){
							 dist = 0;
						} 
					}

				}
				if (dist == 1){
					SS_ref_db[i].DF_pc[l] = SS_ref_db[i].G_pc[l];
					for (j = 0; j < gv.len_ox; j++) {
						SS_ref_db[i].DF_pc[l] -= SS_ref_db[i].comp_pc[l][j]*gv.gam_tot[j];
					}

					if (SS_ref_db[i].DF_pc[l] < min_df){

						if (id_c == n_candidate){ id_c = 0;}
						pc_candidate[id_c] = l;
						df_candidate[id_c] = SS_ref_db[i].DF_pc[l];
						id_c 			  += 1;

						min_df 		= SS_ref_db[i].DF_pc[l];
					}
				}
			}
			id_c -= 1;
			if (id_c == -1){ id_c = n_candidate-1;}

			for (int c = 0; c < n_candidate; c++){
				if (id_c == n_candidate){ id_c = 0;}

				if (df_candidate[id_c] < gv.PC_df_add && pc_candidate[id_c] != -1){

					if(phase_add == 0){

						if (gv.verbose == 1){
							printf("  - %4s %5d added [PC DF check]\n",gv.SS_list[i],pc_candidate[id_c]);
							
							for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
								SS_ref_db[i].iguess[k] = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
							}								
						}
											
						/**
							copy the minimized phase informations to cp structure
						*/
						gv.len_cp				   += 1;
						id_cp 		 				= gv.len_cp-1;
						strcpy(cp[id_cp].name,gv.SS_list[i]);				/* get phase name */				
						cp[id_cp].in_iter			= gv.global_ite;
						cp[id_cp].ss_flags[0] 		= 1;						/* set flags */
						cp[id_cp].ss_flags[1] 		= 0;
						cp[id_cp].ss_flags[2] 		= 1;
						cp[id_cp].split 			= 0;							
						cp[id_cp].id 				= i;						/* get phase id */
						cp[id_cp].n_xeos			= SS_ref_db[i].n_xeos;		/* get number of compositional variables */
						cp[id_cp].n_em				= SS_ref_db[i].n_em;		/* get number of endmembers */
						cp[id_cp].n_sf				= SS_ref_db[i].n_sf;		/* get number of site fractions */
						
						for (k = 0; k < SS_ref_db[i].n_xeos; k++){
							cp[id_cp].dguess[k]     = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
							cp[id_cp].xeos[k]       = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
						}
						for (k = 0; k < SS_ref_db[i].n_xeos; k++){
							cp[id_cp].mu[k]    		=  0.0;
						}

						gv.n_solvi[i] 	       	   += 1;
						gv.id_solvi[i][gv.n_solvi[i]] = id_cp;
						pc_added[phase_add] 		= pc_candidate[id_c];

						phase_add				   += 1;	

						id_c += 1;
					}
					else{
						dist = 1;

						for (m = 0; m < phase_add; m++){
							xeos_dist = euclidean_distance(SS_ref_db[i].xeos_pc[pc_candidate[id_c]], SS_ref_db[i].xeos_pc[pc_added[m]], SS_ref_db[i].n_xeos);
							if (xeos_dist < gv.PC_min_dist*gv.SS_PC_stp[i]*sqrt((double)SS_ref_db[i].n_xeos) ){
								dist = 0;
							} 
						}

						if (dist == 1){

							if (gv.verbose == 1){
								printf("  - %4s %5d added [PC DF check]\n",gv.SS_list[i],pc_candidate[id_c]);
								
								for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
									SS_ref_db[i].iguess[k] = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
								}								
							}
						
							/**
								copy the minimized phase informations to cp structure
							*/
							gv.len_cp				   += 1;
							id_cp 		 				= gv.len_cp-1;
							strcpy(cp[id_cp].name,gv.SS_list[i]);				/* get phase name */				
							cp[id_cp].in_iter			= gv.global_ite;
							cp[id_cp].ss_flags[0] 		= 1;						/* set flags */
							cp[id_cp].ss_flags[1] 		= 0;
							cp[id_cp].ss_flags[2] 		= 1;
							cp[id_cp].split 			= 0;							
							cp[id_cp].id 				= i;						/* get phase id */
							cp[id_cp].n_xeos			= SS_ref_db[i].n_xeos;		/* get number of compositional variables */
							cp[id_cp].n_em				= SS_ref_db[i].n_em;		/* get number of endmembers */
							cp[id_cp].n_sf				= SS_ref_db[i].n_sf;		/* get number of site fractions */
							
							for (k = 0; k < SS_ref_db[i].n_xeos; k++){
								cp[id_cp].dguess[k]     = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
								cp[id_cp].xeos[k]       = SS_ref_db[i].xeos_pc[pc_candidate[id_c]][k];
							}
							for (k = 0; k < SS_ref_db[i].n_xeos; k++){
								cp[id_cp].mu[k]    		=  0.0;
							}

							gv.n_solvi[i] 	       	   += 1;
							gv.id_solvi[i][gv.n_solvi[i]] = id_cp;
							pc_added[phase_add] 		= pc_candidate[id_c];

							phase_add				   += 1;	

							id_c += 1;

						}
						// printf("dist: %d\n",dist);
						
					}

				}


			}
			
		}
	}
	
	// for (i = 0; i < gv.len_ss; i++){
	// 	min_df_id = -1;						// unreallistic index to start with
	// 	min_df    =  1e6;					// high starting value as it is expected to go down
	// 	phase_add =  0;
		
	// 	if (SS_ref_db[i].ss_flags[0] == 1  && gv.verifyPC[i] == 1){
			
	// 		max_n_pc  = ((SS_ref_db[i].tot_pc >= SS_ref_db[i].n_pc) ? (SS_ref_db[i].n_pc) : (SS_ref_db[i].tot_pc));
			
	// 		for (l = 0; l < max_n_pc; l++){

	// 			dist =  1;
	// 			if (gv.n_solvi[i] > 0){

	// 				for (k = 0; k < gv.n_solvi[i]; k++){  /* go through the upper triangle of the matrix (avoiding diagonal)*/
	// 					ph = SS_ref_db[i].solvus_id[k];

	// 					xeos_dist = euclidean_distance(cp[ph].xeos, SS_ref_db[i].xeos_pc[l], SS_ref_db[i].n_xeos);
	// 					if (xeos_dist < gv.PC_min_dist*gv.SS_PC_stp[i]*sqrt((double)SS_ref_db[i].n_xeos) ){
	// 						 dist = 0;
	// 					} 
	// 				}	
	// 			}
	// 			if (dist == 1){
	// 				SS_ref_db[i].DF_pc[l] = SS_ref_db[i].G_pc[l];
	// 				for (j = 0; j < gv.len_ox; j++) {
	// 					SS_ref_db[i].DF_pc[l] -= SS_ref_db[i].comp_pc[l][j]*gv.gam_tot[j];
	// 				}

	// 				if (SS_ref_db[i].DF_pc[l] < min_df){	
	// 					min_df 		= SS_ref_db[i].DF_pc[l];
	// 					min_df_id 	= l;
	// 				}
	// 			}
	// 		}
			
	// 		/* if there is a possible solvus */
	// 		if (min_df < gv.PC_df_add && min_df_id != -1 && phase_add < 2){
	// 			if (gv.verbose == 1){
	// 				printf("  - %4s %5d added [PC DF check]\n",gv.SS_list[i],min_df_id);
					
	// 				for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
	// 					SS_ref_db[i].iguess[k] = SS_ref_db[i].xeos_pc[min_df_id][k];
	// 				}								
	// 			}

	// 			for (k = 0; k < SS_ref_db[i].n_xeos; k++) {
	// 				SS_ref_db[i].iguess[k] = SS_ref_db[i].xeos_pc[min_df_id][k];
	// 			}
				 												 
	// 			/**
	// 				copy the minimized phase informations to cp structure
	// 			*/
	// 			gv.len_cp				   += 1;
	// 			id_cp 		 				= gv.len_cp-1;
	// 			strcpy(cp[id_cp].name,gv.SS_list[i]);				/* get phase name */				
	// 			cp[id_cp].in_iter			= gv.global_ite;
	// 			cp[id_cp].ss_flags[0] 		= 1;						/* set flags */
	// 			cp[id_cp].ss_flags[1] 		= 0;
	// 			cp[id_cp].ss_flags[2] 		= 1;
	// 			cp[id_cp].split 			= 0;							
	// 			cp[id_cp].id 				= i;						/* get phase id */
	// 			cp[id_cp].n_xeos			= SS_ref_db[i].n_xeos;		/* get number of compositional variables */
	// 			cp[id_cp].n_em				= SS_ref_db[i].n_em;		/* get number of endmembers */
	// 			cp[id_cp].n_sf				= SS_ref_db[i].n_sf;		/* get number of site fractions */
				
	// 			for (k = 0; k < SS_ref_db[i].n_xeos; k++){
	// 				cp[id_cp].dguess[k]     = SS_ref_db[i].xeos_pc[min_df_id][k];
	// 				cp[id_cp].xeos[k]       = SS_ref_db[i].xeos_pc[min_df_id][k];
	// 			}
	// 			for (k = 0; k < SS_ref_db[i].n_xeos; k++){
	// 				cp[id_cp].mu[k]    		=  0.0;
	// 			}

	// 			gv.n_solvi[i] 	       	   += 1;
	// 			gv.id_solvi[i][gv.n_solvi[i]] = id_cp;
	// 			phase_add				   += 1;
	// 		}
			
	// 	}
	// }

	return gv;
};



/**
	checks if the pseudocompounds generated during the levelling stage yield a negative driving force
*/
global_variable check_PC_driving_force(		bulk_info 	 z_b,
											global_variable  	 gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){

	int max_n_pc, n_em;
	printf("\n");
	for (int i = 0; i < gv.len_ss; i++){
		if (SS_ref_db[i].ss_flags[0] == 1){
				
			n_em 	 = SS_ref_db[i].n_em;
			max_n_pc = ((SS_ref_db[i].tot_pc >= SS_ref_db[i].n_pc) ? (SS_ref_db[i].n_pc) : (SS_ref_db[i].tot_pc));
			
			for (int l = 0; l < max_n_pc; l++){
				SS_ref_db[i].DF_pc[l] = SS_ref_db[i].G_pc[l];
				for (int j = 0; j < gv.len_ox; j++) {
					SS_ref_db[i].DF_pc[l] -= SS_ref_db[i].comp_pc[l][j]*gv.gam_tot[j];
				}
				
				if (SS_ref_db[i].DF_pc[l] < -1e-10){
					printf("%4s #%4d | %+10f | ",gv.SS_list[i],l,SS_ref_db[i].DF_pc[l]);
					for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
						printf(" %+10f",SS_ref_db[i].xeos_pc[l][k]);
					}
					for (int k = SS_ref_db[i].n_xeos; k < 11; k++){
						printf(" %10s","-");
					}

					printf("\n");
				}
			}	
		}
	}

	return gv;
};


/**
	save initial conditions for the extended Newton method
*/
global_variable save_ini_conditions(		global_variable  	 gv,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){
	for (int ii = 0; ii < gv.len_ox; ii++){
		gv.gam_tot_0[ii] = gv.gam_tot[ii];
	}

	/* Update pure phase (PP) fractions */
	if (gv.n_pp_phase > 0){
		for (int ii = 0; ii < gv.n_pp_phase; ii++){
			gv.pp_n_0[gv.pp_id[ii]]  = gv.pp_n[gv.pp_id[ii]];
		}
	}

	int ph_id;
	for (int iss = 0; iss < gv.len_cp; iss++){ 
		ph_id = cp[iss].id;
		if (cp[iss].ss_flags[0] == 1){
			cp[iss].ss_n_0 = cp[iss].ss_n;

			for (int ii = 0; ii < cp[iss].n_xeos; ii++){
				cp[iss].xeos_0[ii]	= cp[iss].xeos[ii];
			}
		}
	}

	return gv;
}


/**
	save initial conditions for the extended Newton method
*/
global_variable load_ini_conditions(		global_variable  	 gv,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){

	int ph_id;
	for (int iss = 0; iss < gv.len_cp; iss++){ 
		ph_id = cp[iss].id;
		if (cp[iss].ss_flags[0] == 1){
			cp[iss].ss_n_0 = cp[iss].ss_n;

			for (int ii = 0; ii < cp[iss].n_xeos; ii++){
				cp[iss].xeos[ii] = cp[iss].xeos_0[ii];
			}
		}
	}

	return gv;
}




/** 
	function to fill LHS (J)
*/
void PGE_build_Jacobian( 	double 			    *A,
							bulk_info 	 z_b,
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
				A[ix] +=  SS_ref_db[ss].Comp[i][z_b.nzEl_array[j]] * cp[ph].factor * (cp[ph].p_em[i]*cp[ph].xi_em[i])  * SS_ref_db[ss].z_em[i];		
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
										bulk_info 	 z_b,
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
	max_dG_ss   = gv.relax_PGE_val*exp(-8.0*pow(gv.BR_norm,0.29))+1.0;

	g_fac       = (gv.max_g_phase/max_dG_ss)/max_dG;
	n_fac   	= (gv.max_n_phase/max_dG_ss)/max_dn;
	alpha 		= ((n_fac < g_fac) 	 	? 	(n_fac) 		: (g_fac)	);
	alpha 		= ((alpha > gv.max_fac) ? 	(gv.max_fac) 	: (alpha)	);

	gv.alpha	= alpha;
	
	/* Update Gamma */
	for (i = 0; i < z_b.nzEl_val; i++){
		gv.delta_gam_tot[z_b.nzEl_array[i]]  = gv.dGamma[i]*alpha;	
		gv.gam_tot[z_b.nzEl_array[i]] 		+= gv.dGamma[i]*alpha;		
	}
	
	gv.gamma_norm[gv.global_ite] = norm_vector(gv.dGamma, z_b.nzEl_val);

	/* Update solusion phase (SS) fractions */
	for (i = 0; i < gv.n_cp_phase; i++){
		 cp[gv.cp_id[i]].delta_ss_n  = gv.dn_cp[i]*alpha;
		 cp[gv.cp_id[i]].ss_n 		+= gv.dn_cp[i]*alpha;
	}
	
	/* Update pure phase (PP) fractions */
	if (gv.n_pp_phase > 0){
		for (i = 0; i < gv.n_pp_phase; i++){
			 gv.pp_n[gv.pp_id[i]] 		+= gv.dn_pp[i]*alpha;
			 gv.delta_pp_n[gv.pp_id[i]]  = gv.dn_pp[i]*alpha;
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
	int 	ipiv[nEntry];												/** pivot indices*/
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

	/**
		save RHS vector 
	*/
	gv.fc_norm_t1 = norm_vector(	gv.b_PGE,
									nEntry		);
		

	/**
		call lapacke to solve system of linear equation using LU 
	*/
	info = LAPACKE_dgesv(		LAPACK_ROW_MAJOR, 
								nEntry, 
								nrhs, 
								gv.A_PGE, 
								lda, 
								ipiv, 
								gv.b_PGE, 
								ldb					);

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

		/* Update mu and xi of solution phase  */

		gv =	PGE_update_pi(				z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp						);  

		gv =	PGE_update_xi(				z_b,								/** bulk rock constraint 				*/ 
											gv,									/** global variables (e.g. Gamma) 		*/

											PP_ref_db,							/** pure phase database 				*/ 
											SS_ref_db,							/** solution phase database 			*/
											cp						);  

		// add_PGE_pseudocompounds(			z_b,
		// 									splx_data,
		// 									gv,
												
		// 									PP_ref_db,
		// 									SS_ref_db				);

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
  Main PGE routine
*/ 
global_variable LP(		bulk_info 			z_b,
						global_variable 	gv,

						obj_type 			*SS_objective,
						simplex_data	    *splx_data,
						PP_ref 				*PP_ref_db,
						SS_ref 				*SS_ref_db,
						csd_phase_set  		*cp					){
		
	clock_t t; 	

	gv.LP 	 = 1;	gv.PGE 	 = 0;

	int mode = 1;

	// for (int gi = 0; gi < 32; gi++){				
	while ( gv.global_ite < 32 ){

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

		if (gv.global_ite == 12 || gv.global_ite == 24){
			gv = check_PC( 				z_b,						/** bulk rock constraint 				*/ 
										gv,							/** global variables (e.g. Gamma) 		*/

										PP_ref_db,					/** pure phase database 				*/ 
										SS_ref_db,
										cp					); 
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
		ss_min_LP(						mode,
										gv, 							/** global variables (e.g. Gamma) 		*/

										SS_objective,							
										z_b,							/** bulk-rock, pressure and temperature conditions */
										SS_ref_db,						/** solution phase database 			*/	
										cp 					);

		/**
		   Here the linear programming method is used after the PGE step to get a new Gibbs hyper-plane
		*/
		gv = run_LP_with_PGE_phase(		z_b,
										splx_data,
										gv,
												
										PP_ref_db,
										SS_ref_db			);


		gv = init_PGE_using_LP(			z_b,
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
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[0] == 1){
				if (isnan(cp[i].df) == 1 || isinf(cp[i].df) == 1){
					gv.div = 1;	
				}
			}
		}
		if (gv.div == 1){ gv.status = 4; break; }

		t = clock() - t; 
		if (gv.verbose == 1){
			printf("\n __ iteration duration: %+4f ms __\n\n\n",((double)t)/CLOCKS_PER_SEC*1000);
		}
		gv.ite_time[gv.global_ite] = ((double)t)/CLOCKS_PER_SEC*1000;
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
		
	clock_t t; 	

	gv.LP_PGE_switch = gv.global_ite;
	gv.LP 	= 0;	gv.PGE 	= 1;

	int mode = 1;
	// for (int gi = 0; gi < 1; gi++){	
	while (gv.BR_norm > gv.br_max_tol || gv.outter_PGE_ite > (gv.global_ite - gv.LP_PGE_switch)){

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
		if (gv.BR_norm < gv.PC_check_val1 && gv.check_PC1 == 0){
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
		}
		/**
			check driving force of PC when getting close to convergence
		*/
		if (gv.BR_norm < gv.PC_check_val2 && gv.check_PC2 == 0){
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
		gv = split_cp(						gv, 						/** global variables (e.g. Gamma) 		*/
											SS_ref_db,					/** solution phase database 			*/	
											cp 				);		
		/**
			Local minimization of the solution phases
		*/
		ss_min_PGE(							mode,
											gv, 						/** global variables (e.g. Gamma) 		*/

											SS_objective,							
											z_b,						/** bulk-rock, pressure and temperature conditions */
											SS_ref_db,					/** solution phase database 			*/	
											cp 				);	

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

		/**
		 * Solver status
		 * 0: success
		 * 1: under-relaxed
		 * 2: more under-relaxed
		 * 3: reached max iterations (failed)
		 * 4: terminated due to slow convergence or a very large residual (failed)
		**/
		if (gv.global_ite > gv.it_1 && gv.BR_norm < gv.br_max_tol*gv.ur_1){		if (gv.verbose != -1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_1, gv.ur_1);}; gv.status = 1; break;}
		if (gv.global_ite > gv.it_2 && gv.BR_norm < gv.br_max_tol*gv.ur_2){		if (gv.verbose != -1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_2, gv.ur_2);}; gv.status = 2; break;}
		if (gv.global_ite > gv.it_3 && gv.BR_norm < gv.br_max_tol*gv.ur_3){		if (gv.verbose != -1){printf(" >%d iterations, under-relax mass constraint norm (*%.1f)\n\n", gv.it_3, gv.ur_3);}; gv.status = 2; break;}
		if (gv.global_ite > gv.it_f){											if (gv.verbose != -1){printf(" >%d iterations, did not converge  !!!\n\n", gv.it_f);}; gv.status = 3; break;}

		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[0] == 1){
				if (isnan(cp[i].df) == 1 || isinf(cp[i].df) == 1){
					gv.div = 1;	
				}
			}
		}
		if (gv.div == 1){ gv.status = 4; break; }

		t = clock() - t; 
		if (gv.verbose == 1){
			printf("\n __ iteration duration: %+4f ms __\n\n\n",((double)t)/CLOCKS_PER_SEC*1000);
		}

		/* check evolution of mass constraint residual */
		gv.PGE_mass_norm[gv.global_ite]  = gv.BR_norm;				/** save norm for the current global iteration */
		gv.Alg[gv.global_ite] 			 = 1;
		gv.ite_time[gv.global_ite] 		 = ((double)t)/CLOCKS_PER_SEC*1000;
		gv.global_ite 			  		+= 1;

	}

	return gv;
};		

