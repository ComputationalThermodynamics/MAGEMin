/**

The objectives of this function are multiple:  
            
- remove phase when fraction < 0.0                        
- add phase when df < 0.0 and sum_xi > 1.0 (for solution phase only)                               
- swap phase when n_phase = n_oxides and df < 0.0    
- make sure pure phase polymorph are correctly accounted for   
- update phase fraction when adding/removing phase using least square optimization                                  
                                                           
The core of the function revolves around book-keeping and  
the informations stored in PP_flags and SS_flags arrays 
   
Those arrays store the state of the phases:                

PP & SS_flags
-------------

+-------+-------+------+------+-------+-------+------+
| SS/PP |  IN	| CSD  | HLD  |  RMV  |  CYC  | REIN |
+=======+=======+======+======+=======+=======+======+
| [0]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [1]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [2]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| .     | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [m/n]	| 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
										 
- IN:   allowed phase (satisfying bulk rock constraints)
- CSD:  considered phase (part of the active set of phases)
- HLD:  on hold (not in the active set but still scanned at every iteration)
- RMV:  removed (not considered anymore)
- CYC:  number of cycle 
- REIN: phase reintroduced

            
NOTES:
22/06/2021 -> This function should be reworked to make it more efficient/readable
            
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "toolkit.h"

/**
  structure to be passed to compare function 
*/
struct str
{
	double value;
	int    index;
};

/**
  compare double function
*/
int cmp_dbl(const void *a, const void *b)
{
	struct str *a1 = (struct str *)a;
	struct str *a2 = (struct str *)b;
	
	if ((*a1).value < (*a2).value)
		return -1;
	else if ((*a1).value > (*a2).value)
		return 1;
	else
		return 0;
}
	
/**
  compare int function 
*/
int cmp_int(const void *a, const void *b)
{
	int a1 = *((int*)a);
	int a2 = *((int*)b);
	
	if (a1 < a2)
		return -1;
	else if (a1 > a2)
		return 1;
	else
		return 0;
}


/**
	check PC driving force and add phase if below hyperplane
*/
global_variable check_PC(					bulk_info 	 		 z_b,
											global_variable  	 gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
){
	double 	min_df, xeos_dist, norm;
	int 	max_n_pc, phase_add, id_cp, dist, ph,ph_id;
	int 	i,j,k,l,c,m;

	int 	n_candidate = 8;
	int     pc_candidate[n_candidate];
	int     pc_added[n_candidate];
	double  df_candidate[n_candidate];
	int     id_c;

	/* section to add anti-ordering counterpart of active solution phases */
	int add_max = 0;
	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			ph_id = cp[i].id;

			if (SS_ref_db[ph_id].orderVar == 1 && add_max < 2){
				gv.len_cp				   += 1;
				id_cp 		 				= gv.len_cp-1;
				strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */				
				cp[id_cp].in_iter			= gv.global_ite;
				cp[id_cp].ss_flags[0] 		= 1;						/* set flags */
				cp[id_cp].ss_flags[1] 		= 0;
				cp[id_cp].ss_flags[2] 		= 1;
				cp[id_cp].split 			= 0;							
				cp[id_cp].id 				= ph_id;						/* get phase id */
				cp[id_cp].n_xeos			= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
				cp[id_cp].n_em				= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
				cp[id_cp].n_sf				= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
				
				for (k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
					cp[id_cp].dguess[k]     = cp[i].xeos[k]*SS_ref_db[ph_id].idOrderVar[k];
					cp[id_cp].xeos[k]       = cp[i].xeos[k]*SS_ref_db[ph_id].idOrderVar[k];
				}
				for (k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
					cp[id_cp].mu[k]    		=  0.0;
				}

				if (gv.verbose == 1){
						printf("  - %4s    anti-ordering counter-part [act phase]",gv.SS_list[ph_id]);
					
					for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++) {
						printf(" %+8f",cp[id_cp].dguess[k]);
					}	
					printf("\n");							
				}

				add_max += 1;

			}


		}
	}


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
			for (l = 0; l < SS_ref_db[i].tot_pc; l++){
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
						pc_added[phase_add] 		= pc_candidate[id_c];

						phase_add				   += 1;	

						if (gv.verbose == 1){
							printf("  - %4s %5d, DF: %+10f added [PC DF check]",gv.SS_list[i],pc_candidate[id_c],df_candidate[id_c]);
							
							for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
								printf(" %+8f",cp[id_cp].dguess[k]);
							}	
							printf("\n");							
						}

						id_c += 1;





						if (SS_ref_db[i].orderVar == 1){
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
								cp[id_cp].dguess[k]     = cp[id_cp-1].dguess[k]*SS_ref_db[i].idOrderVar[k];
								cp[id_cp].xeos[k]       = cp[id_cp-1].xeos[k]*SS_ref_db[i].idOrderVar[k];
							}
							for (k = 0; k < SS_ref_db[i].n_xeos; k++){
								cp[id_cp].mu[k]    		=  0.0;
							}


							if (gv.verbose == 1){
									printf("     anti-ordering counterpart:");
								
								for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
									printf(" %+8f",cp[id_cp].dguess[k]);
								}	
								printf("\n");							
							}

						}



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
							pc_added[phase_add] 		= pc_candidate[id_c];

							phase_add				   += 1;	


							if (gv.verbose == 1){
								printf("  - %4s %5d, DF: %+10f added [PC DF check]",gv.SS_list[i],pc_candidate[id_c],df_candidate[id_c]);
								
								for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
									printf(" %+8f",cp[id_cp].dguess[k]);
								}	
								printf("\n");							
							}
							id_c += 1;


							if (SS_ref_db[i].orderVar == 1){
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
									cp[id_cp].dguess[k]     = cp[id_cp-1].dguess[k]*SS_ref_db[i].idOrderVar[k];
									cp[id_cp].xeos[k]       = cp[id_cp-1].xeos[k]*SS_ref_db[i].idOrderVar[k];
								}
								for (k = 0; k < SS_ref_db[i].n_xeos; k++){
									cp[id_cp].mu[k]    		=  0.0;
								}

								if (gv.verbose == 1){
									printf("     anti-ordering counterpart:");
									
									for (int k = 0; k < SS_ref_db[i].n_xeos; k++) {
										printf(" %+8f",cp[id_cp].dguess[k]);
									}	
									printf("\n");							
								}

							}
						}

					}

				}


			}
			
		}
	}

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

			for (int l = 0; l < SS_ref_db[i].tot_pc; l++){
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
  Merge solution phase routine 
*/			
global_variable phase_merge_function(		bulk_info 			z_b,
											global_variable 	gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp 
){
	
	if (gv.verbose == 1){
		printf("\nMerge Compositionally close solution phases\n");
		printf("════════════════════════════════════════════\n");
		printf(" phase |  #cp > #cp | Euclidian distance\n");
	}
	
	int ph_id, i, j, k, l, iss, phA, phB;
	double distance;

	/* reinitialize the number of SS instances */
	for (iss = 0; iss < gv.len_ss; iss++){
		gv.n_solvi[iss] = 0;
	}

	/* get number of duplicated phases and their cp id */
	for (i = 0; i < gv.len_cp; i++){
		ph_id = cp[i].id;
		if (cp[i].ss_flags[0] == 1 ){
			SS_ref_db[ph_id].solvus_id[gv.n_solvi[ph_id]] = i;
			gv.n_solvi[ph_id] += 1;
		}
	}

	/* check and merge phases in the active assemblage */	
	for (iss = 0; iss < gv.len_ss; iss++){
		
		/* if there is a possible solvus */
		if (gv.n_solvi[iss] > 1){
			
			/* go through the upper triange of the matrix (avoiding diagonal)*/
			for (k = 0; k < gv.n_solvi[iss]; k++){ 
				for (l = k+1; l < gv.n_solvi[iss]; l++){
					phA = SS_ref_db[iss].solvus_id[k];
					phB = SS_ref_db[iss].solvus_id[l];
					
					/* if the two instances of the same phase can still be considered (not removed from previous iteration)*/
					if (phA != -1 && phB != -1){
						
						/* get the norm of the Euclidian distance */
						distance = euclidean_distance( cp[phA].p_em, cp[phB].p_em, SS_ref_db[iss].n_em);

						/* if the distance is lower than the one given in initiialize.h then proceed to merge */
						if (distance < gv.merge_value){

							/** case 1: both instances are either on ACTIVE or on HOLD */ 
							if (cp[phA].ss_flags[1] + cp[phB].ss_flags[1] != 1){
								if (gv.verbose ==1){
									printf(" %5s | %1d.%1d > %1d.%1d  | %+10f\n",gv.SS_list[iss],l,cp[phB].ss_flags[1],k,cp[phA].ss_flags[1],distance);
								}
								/* in case both phases are active, add phase B fraction on phase A */
								if(cp[phA].ss_flags[1] == 1 && cp[phB].ss_flags[1] == 1){
									cp[phA].ss_n 		  	+=  cp[phB].ss_n;
									/* merge compositional variables at mid point */
									for (i = 0; i < cp[phA].n_xeos; i++){
										cp[phA].xeos[i] = (cp[phA].xeos[i] + cp[phB].xeos[i])/2.0;
									}

									gv.n_cp_phase 		  	-=  1;
									gv.n_phase            	-=  1;
								}
								cp[phB].ss_flags[0]    		 =  0;
								cp[phB].ss_flags[1]    		 =  0;
								cp[phB].ss_flags[2]    		 =  0;
								cp[phB].ss_n		   		 =  0.0;
								SS_ref_db[iss].solvus_id[l]  = -1;	
							}
							else{
								/* Phase A is ACTIVE */
								if (cp[phA].ss_flags[1] == 1){
									if (gv.verbose ==1){
										printf(" %5s | %1d.%1d > %1d.%1d  | %+10f\n",gv.SS_list[iss],l,cp[phB].ss_flags[1],k,cp[phA].ss_flags[1],distance);
									}
									cp[phB].ss_flags[0]    		 =  0;
									cp[phB].ss_flags[1]    		 =  0;
									cp[phB].ss_flags[2]    		 =  0;
									cp[phB].ss_n		   		 =  0.0;
									SS_ref_db[iss].solvus_id[l]  = -1;	
								}
								else{
									/* Phase B is ACTIVE */
									if (gv.verbose ==1){
										printf(" %5s | %1d.%1d > %1d.%1d  | %+10f\n",gv.SS_list[iss],l,cp[phA].ss_flags[1],k,cp[phB].ss_flags[1],distance);
									}
									cp[phA].ss_flags[0]    		 =  0;
									cp[phA].ss_flags[1]    		 =  0;
									cp[phA].ss_flags[2]    		 =  0;
									cp[phA].ss_n		   		 =  0.0;
									SS_ref_db[iss].solvus_id[k]  = -1;		
								}
							}
						}
						
					}
					
				}
			}
		}
	}
	
	/* reinitialize the number of SS instances */
	for (iss = 0; iss < gv.len_ss; iss++){
		gv.n_solvi[iss] = 0;
	}

	/* get number of duplicated phases and their cp id */
	for (i = 0; i < gv.len_cp; i++){
		ph_id = cp[i].id;
		if (cp[i].ss_flags[0] == 1 ){
			SS_ref_db[ph_id].solvus_id[gv.n_solvi[ph_id]] = i;
			gv.n_solvi[ph_id] += 1;
		}
	}
	
   return gv;
};

/**
	from active to hold function
*/
global_variable phase_act2hold(			bulk_info 	z_b,
										global_variable 	gv,

										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp 
){

	/** REMOVE AND DELETE PHASES FROM CONSIDERATION IF FRACTION < 0.0**/
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1 && gv.ph_change == 0){
			if (gv.pp_n[i] < 0.0 ){				
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 1;
				gv.pp_n[i]        = 0.0;
				gv.n_pp_phase    -= 1;
				gv.n_phase       -= 1;
				gv.ph_change      = 1;																/** put to 0 if you want to allow multiple phase removal on top of 1 phase addition */
			}
		}
	}
	
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1  && gv.ph_change == 0){
			if (cp[i].ss_n <= 0.0){
				cp[i].ss_flags[1] = 0;
				cp[i].ss_flags[2] = 1;
				cp[i].ss_n        = 0.0;
				gv.n_cp_phase    -= 1;
				gv.n_phase       -= 1;		
				gv.ph_change      = 1;																	/** phase has been removed, then do not add phase during this iteration */
			}
		}
	}

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1  && gv.ph_change == 0){

			if (cp[i].ss_n < 1e-3 && cp[i].df > 1e-3 && cp[i].sum_xi < 1.0){
				cp[i].ss_flags[1] = 0;
				cp[i].ss_flags[2] = 1;
				cp[i].ss_n        = 0.0;
				gv.n_cp_phase    -= 1;
				gv.n_phase       -= 1;		
				gv.ph_change      = 1;																	/** phase has been removed, then do not add phase during this iteration */
			}
		}
	}

	return gv;
}

/**
	from active to hold function
*/
global_variable phase_hold2rmv(			bulk_info 	z_b,
										global_variable 	gv,

										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp 
){

	/** REMOVE AND DELETE PHASES FROM CONSIDERATION IF DELTA G IS GETTING TOO LARGE**/
	// for (int i = 0; i < gv.len_pp; i++){
	// 	if (gv.pp_flags[i][2] == 1){
	// 		if (PP_ref_db[i].gb_lvl*PP_ref_db[i].factor > gv.bnd_filter_pc){				
	// 			gv.pp_flags[i][0] = 0;
	// 			gv.pp_flags[i][1] = 0;
	// 			gv.pp_flags[i][2] = 0;
	// 			gv.pp_flags[i][3] = 1;
	// 			gv.pp_n[i]        = 0.0;															/** put to 0 if you want to allow multiple phase removal on top of 1 phase addition */
	// 		}
	// 	}
	// }
	
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[2] == 1){
			if (cp[i].df*cp[i].factor > gv.bnd_filter_pge){
				cp[i].ss_flags[0] = 0;
				cp[i].ss_flags[1] = 0;
				cp[i].ss_flags[2] = 0;
				cp[i].ss_flags[3] = 1;
				cp[i].ss_n        = 0.0;																/** phase has been removed, then do not add phase during this iteration */
			}
		}
	}

	return gv;
}

/**
	from active to hold function
*/
global_variable phase_hold2act(		bulk_info 				z_b,
									global_variable 		gv,

									PP_ref 					*PP_ref_db,
									SS_ref 					*SS_ref_db,
									csd_phase_set  			*cp 
){
	double 	dnorm     	= 0.0;			/** max driving force under which a phase can be considered to be added 				*/
	double 	min_sumxi 	= 1.0;				/** min sum of xi fractions over which a solution phase can be considered to be added 	*/
	
	/* get number of hold phase for pure phases */
	int n_pp_hld = 0;
	int pp_act[gv.len_ox];
	int inc = 0;
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			pp_act[inc] = i;
			inc        +=1;
		}
		if (gv.pp_flags[i][2] == 1 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor < gv.min_df){
			n_pp_hld += 1;
		}
	}
	
	/* get number of hold phase for solution phases */
	int n_cp_hld = 0;
	int cp_act[gv.len_ox];
	for (int i = 0; i < gv.len_ox; i++){ cp_act[i] = 0; }
	
	inc = 0;
	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[1] == 1){
			cp_act[inc]            = i;
			inc                   += 1;
		}
		if (cp[i].ss_flags[2] == 1 && cp[i].df*cp[i].factor < gv.min_df && cp[i].sf_ok == 1){
			n_cp_hld              += 1;
		}		
	}

	/** -----------------------------------SORTING PURE AND SOLUTION PHASES BY DRIVING FORCES------------------------------------------------------------------------- **/		
	/* create the structures that will hold the phase array sorted by driving force */
	struct str hld_cp_sort[n_cp_hld];
	struct str hld_pp_sort[n_pp_hld];
	
	inc = 0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[2] == 1 && cp[i].df*cp[i].factor < gv.min_df && cp[i].sf_ok == 1){	
			hld_cp_sort[inc].value  = cp[i].df*cp[i].factor;
			hld_cp_sort[inc].index  = i;
			inc += 1;
		}
	}
	
	/* sort ss array using str structure containing, DF values and phase indices */
	qsort(hld_cp_sort, n_cp_hld, sizeof(hld_cp_sort[0]), cmp_dbl);
	qsort(cp_act, sizeof(cp_act)/sizeof(*cp_act), sizeof(*cp_act), cmp_int);
	
	inc = 0;
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][2] == 1 && PP_ref_db[i].gb_lvl*PP_ref_db[i].factor < gv.min_df){	
			hld_pp_sort[inc].value = PP_ref_db[i].gb_lvl*PP_ref_db[i].factor;
			hld_pp_sort[inc].index = i;
			inc += 1;
		}
	}
	/* sort pp array using str structure containing, DF values and phase indices */
	qsort(hld_pp_sort, n_pp_hld, sizeof(hld_pp_sort[0]), cmp_dbl);
	qsort(pp_act, sizeof(pp_act)/sizeof(*pp_act), sizeof(*pp_act), cmp_int);

	/** 
		ADD NEW SOLUTION PHASE TO THE SYSTEM 
	*/
	for (int i = 0; i < n_cp_hld; i++){
		int 	ixs = hld_cp_sort[i].index;
		double 	df  = hld_cp_sort[i].value;
																			/** loop through sorted SS in hold 														*/
		/** if driving force is negative, phase can be potentially added to the system and decrease Gibbs */
		if (df < gv.min_df){
			/** if phase was never previously added to the system */
			if ( gv.ph_change == 0 ){	

				cp[ixs].ss_flags[1] 	  		  = 1;						/** set to active 																		*/
				cp[ixs].ss_flags[2]	  			  = 0;						/** reset hold 																			*/
				cp[ixs].ss_n       	 		 	  = gv.re_in_n;				/** set initial fraction 																*/
				
				gv.n_cp_phase    				 += 1;						/** set new number of active SS 														*/
				gv.n_phase      				 += 1;						/** set new number of total active phases 												*/
				gv.ph_change 					  = 1;						/** a phase change has been achieved during the iteration 								*/
			}
		}
	}

	/* ADD PURE PHASE TO CURRENT ASSEMBLAGE */	
	int is_polymorph =  0;													/** if = 0, pure phase is a polymorph else is not a polymorph 							*/
	int id_polymorph = -1;													/** save index of the polymorph pure phase 												*/
	for (int i = 0; i < n_pp_hld; i++){ 									/** loop through sorted SS in hold 														*/
		int ixp = hld_pp_sort[i].index;
					
		/* list through pure phases ordered by increasing dG */
		if (hld_pp_sort[i].value < 0.0){									/** if driving force is negative, phase can be potentially added to the system and decrease Gibbs */
			if (gv.ph_change == 0 ){										/** if phase can be potentially added to the system 									*/
				if (gv.n_pp_phase > 0){											/** if a pure phase is already in the active set of phases 								*/
					/* check if pure phase to add is a polymorph of one of active pure phase */
					for (int k = 0; k < gv.len_pp; k++){
						if (gv.pp_flags[k][1] == 1 ){							/** compare pure phase to add with pure phase in the active set 						*/
							is_polymorph  = 0;
							for (int l = 0; l < gv.len_ox; l++){
								if (PP_ref_db[k].Comp[l] != PP_ref_db[ixp].Comp[l]){
									is_polymorph = 1;
									break;
								}
							}
							if (is_polymorph == 0){
								id_polymorph  = k;
								break;
							}
						}						
					}
					
					if (is_polymorph != 0){									/** if the pure phase to add is not a polymorph to an active phase, then just add it 	*/
						gv.pp_flags[ixp][1]  = 1;							/** set to active 																		*/
						gv.pp_flags[ixp][2]  = 0;							/** reset hold 																			*/
						gv.pp_n[ixp]         = gv.re_in_n;					/** set initial fraction 																*/
						gv.n_pp_phase    	+= 1;							/** set new number of active PP 														*/
						gv.n_phase      	+= 1;							/** set new number of total active phases 												*/
						gv.ph_change 		 = 1;							/** a phase change has been achieved during the iteration 								*/
					}
					else{
						if (PP_ref_db[ixp].gb_lvl < PP_ref_db[id_polymorph].gb_lvl){
							gv.pp_flags[ixp][1]  = 1;						/** set to active 																		*/
							gv.pp_flags[ixp][2]  = 0;						/** reset hold 																			*/
							gv.pp_n[ixp]         = gv.pp_n[id_polymorph];	/** set initial fraction to initial polymorph 											*/
							gv.pp_flags[id_polymorph][1]          = 0;		/** set initial polymorph to inactive 													*/
							gv.pp_flags[id_polymorph][2]          = 0;		/** reset hold 																			*/
							gv.pp_flags[id_polymorph][3]          = 1;		/** remove initial polymorph 															*/
							gv.pp_n[id_polymorph]                 = 0.0;	/** reset initial polymorph fraction to 0.0 											*/
							gv.ph_change 						  = 1;		/** a phase change has been achieved during the iteration 								*/
						}
						else{
							gv.pp_flags[ixp][1]  = 0;						/** set to inactive 																	*/
							gv.pp_flags[ixp][2]  = 0;						/** reset hold 																			*/
							gv.pp_flags[ixp][3]  = 1;						/** remove initial polymorph 															*/
						}
					}	
				}
				else{														/** if no pure phase are in the active set of phases 									*/
					gv.pp_flags[ixp][1]  = 1;								/** set to active 																		*/
					gv.pp_flags[ixp][2]  = 0;								/** reset hold 																			*/
					gv.pp_n[ixp]         = gv.re_in_n;						/** set initial fraction 																*/
					gv.n_pp_phase    	+= 1;								/** set new number of active PP 														*/
					gv.n_phase      	+= 1;								/** set new number of total active phases 												*/
					gv.ph_change 		 = 1;								/** a phase change has been achieved during the iteration 								*/
				}
			}
		}
	}

	return gv;
}

/**
  main phase update routine 
*/			
global_variable phase_update_function(		bulk_info 	z_b,
											global_variable 	gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp 
){

	/* initial phase change flag */
	gv.ph_change   = 0;

	/* remove phases with negative fraction */
	gv = 		phase_hold2rmv(	
				z_b,							/** bulk rock constraint 				*/
				gv,								/** global variables (e.g. Gamma) 		*/

				PP_ref_db,						/** pure phase database 				*/
				SS_ref_db,						/** solution phase database 			*/ 
				cp
	); 	

	/* remove phases with negative fraction */
	gv = 		phase_act2hold(	
				z_b,							/** bulk rock constraint 				*/
				gv,								/** global variables (e.g. Gamma) 		*/

				PP_ref_db,						/** pure phase database 				*/
				SS_ref_db,						/** solution phase database 			*/ 
				cp
	); 

	/* check if a phase can be added, and add it */
	if (gv.n_phase < z_b.nzEl_val){	
		gv = 	phase_hold2act(	
				z_b,							/** bulk rock constraint 				*/
				gv,								/** global variables (e.g. Gamma) 		*/

				PP_ref_db,						/** pure phase database 				*/
				SS_ref_db,						/** solution phase database 			*/ 
				cp
		); 
	}

	
   return gv;
};
