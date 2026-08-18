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
Function to call solution phase Minimization        
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#ifdef USE_MPI
	#include "mpi.h"
#endif

#include "nlopt.h"
#include "MAGEMin.h"
#include "gem_function.h"
#include "dump_function.h"
#include "toolkit.h"
#include "GH_database/GH_gem_function.h"
#include "phase_update_function.h"
#include "all_solution_phases.h"

#define LIQ_PC_SYNTH_MAX_DIM 16
#define MAX_LIQ_PHASES 8

/** 
Function to update xi and sum_xi during local minimization.
*/
SS_ref SS_UPDATE_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								bulk_info 	 		 z_b,
								char    			*name){

	/* sf_ok?*/
	if (strcmp(gv.research_group, "tc") 	== 0 ){
		SS_ref_db.sf_ok = 1;
		for (int i = 0; i < SS_ref_db.n_sf; i++){
			if (SS_ref_db.sf[i] < 0.0 || isnan(SS_ref_db.sf[i]) == 1|| isinf(SS_ref_db.sf[i]) == 1){
				SS_ref_db.sf_ok = 0;
				SS_ref_db.sf_id = i;
				break;
			}
		}

		/* xi calculation (phase fraction expression for PGE) */
		SS_ref_db.sum_xi 	= 0.0;
		for (int i = 0; i < SS_ref_db.n_em; i++){
			/* Clamp the exponent's upper bound: for an ordinary phase mu[i]/RT never gets
			   close to this, but fl_DEW's hyperplane-relative mu[i] can be astronomically
			   negative for a formally-favored-but-effectively-absent charged species,
			   overflowing exp() to +Inf; Inf*p[i] with p[i]==0 is then Inf*0=NaN, poisoning
			   sum_xi (and everything downstream in PGE_function.c's mass-balance matrix
			   that reads xi_em). No lower-bound clamp needed - underflow to 0 is the
			   correct, harmless result for a genuinely unfavorable species. */
			SS_ref_db.xi_em[i] = exp(fmin(-SS_ref_db.mu[i]/(SS_ref_db.R*SS_ref_db.T), 700.0));
			SS_ref_db.sum_xi  += SS_ref_db.xi_em[i]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
		}

		/* get composition of solution phase */
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ss_comp[j] = 0.0;
			for (int i = 0; i < SS_ref_db.n_em; i++){
				SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
			} 
		}
	}
	else if (strcmp(gv.research_group, "sb") 	== 0 ){
		SS_ref_db.sf_ok = 1;
	}
	else if (strcmp(gv.research_group, "gh") 	== 0 ){
		SS_ref_db.sf_ok = 1;
	}

	return SS_ref_db;
};


/** 
Function to update xi and sum_xi for the considered phases list (during the inner loop of the PGE stage).
NOTE: When the phase is "liq", the normalization factor is also updated as it depends on the endmember fractions
*/
csd_phase_set CP_UPDATE_function(		global_variable 	gv,
										SS_ref 				SS_ref_db,
										csd_phase_set  		cp, 
										bulk_info 	z_b			){

	/* sf_ok?*/
	cp.sf_ok = 1;
	for (int i = 0; i < cp.n_sf; i++){
		if (cp.sf[i] < 0.0 || isnan(cp.sf[i]) == 1|| isinf(cp.sf[i]) == 1){
			cp.sf_ok = 0;	
			break;
		}
	}
	cp.sum_xi 	= 0.0;
	for (int i = 0; i < cp.n_em; i++){
		/* see the matching clamp + comment in SS_UPDATE_function above */
		cp.xi_em[i] = exp(fmin(-cp.mu[i]/(SS_ref_db.R*SS_ref_db.T), 700.0));
		cp.sum_xi  += cp.xi_em[i]*cp.p_em[i]*SS_ref_db.z_em[i];
	}

	/* get composition of solution phase */
	for (int j = 0; j < gv.len_ox; j++){
		cp.ss_comp[j] = 0.0;
		for (int i = 0; i < cp.n_em; i++){
		   cp.ss_comp[j] += SS_ref_db.Comp[i][j]*cp.p_em[i]*SS_ref_db.z_em[i];
	   } 
	}

	return cp;
};

/** 
	This function ensures that if we drift away from the set of x-eos obtained during levelling, a copy will of the phase will be added to the considered set of phase
	- Drifting occurs when tilting of the hyperplane moves the x-eos far away from their initial guess
	- Note that each instance of the phase, initialized during levelling can only be split once
*/
global_variable split_cp(		global_variable 	 gv,
								SS_ref 			    *SS_ref_db,
								csd_phase_set  		*cp
){
	int id_cp;
	int ph_id;
	double distance;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){

			ph_id= cp[i].id;
			
			distance 	= euclidean_distance( cp[i].xeos, cp[i].dguess, SS_ref_db[ph_id].n_xeos);

			if (distance > 2.0*gv.SS_PC_stp[ph_id]*pow((double)SS_ref_db[ph_id].n_xeos,0.5) && cp[i].split == 0){
				id_cp 					= gv.len_cp;
						
				cp[id_cp].split 		= 1;							/* set split number to one */
				cp[i].split 			= 1;							/* set split number to one */

				strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */	
				
				cp[id_cp].id 			= ph_id;						/* get phaseid */
				cp[id_cp].n_xeos		= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
				cp[id_cp].n_em			= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
				cp[id_cp].n_sf			= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
				
				cp[id_cp].df			= 0.0;
				cp[id_cp].factor		= 0.0;	
				
				cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
				cp[id_cp].ss_flags[1] 	= 0;
				cp[id_cp].ss_flags[2] 	= 1;
				
				cp[id_cp].ss_n          = 0.0;							/* get initial phase fraction */
				
				for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
					cp[id_cp].p_em[ii]      = 0.0;
				}
				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					cp[id_cp].dguess[ii]    = cp[i].dguess[ii];
					cp[id_cp].xeos[ii]      = cp[i].dguess[ii];
					cp[i].dguess[ii]    	= cp[i].xeos[ii];
				}

				gv.n_solvi[ph_id] 	+= 1;
				gv.len_cp 			+= 1;
				
				if (gv.verbose == 1){
					printf("\n  {FYI} %4s cp#%d is grazing away from its field, a copy has been added (xeos = dguess)\n",gv.SS_list[ph_id],i);
				}
				
				if (gv.len_cp == gv.max_n_cp){
					printf(" !! Maxmimum number of allowed phases under consideration reached !!\n    -> check your problem and potentially increase gv.max_n_cp\n");
				}
				
			}
		}	
	}
	return gv;


};

/**
	copy the minimized phase informations to cp structure, if the site fractions are respected
*/
void copy_to_cp(		int 				 i, 
						int 				 ph_id,
						global_variable 	 gv,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp					){

	cp[i].min_time			= SS_ref_db[ph_id].LM_time;
	cp[i].df				= SS_ref_db[ph_id].df_raw;
	cp[i].factor			= SS_ref_db[ph_id].factor;
	cp[i].sum_xi			= SS_ref_db[ph_id].sum_xi;

	for (int ii = 0; ii < cp[i].n_xeos; ii++){
		cp[i].xeos_0[ii]	= cp[i].xeos[ii]; 
		cp[i].xeos[ii]		= SS_ref_db[ph_id].iguess[ii]; 
		cp[i].xeos_1[ii]	= SS_ref_db[ph_id].iguess[ii]; 
		cp[i].dfx[ii]		= SS_ref_db[ph_id].dfx[ii]; 
	}
	
	for (int ii = 0; ii < cp[i].n_em; ii++){
		cp[i].p_em[ii]		= SS_ref_db[ph_id].p[ii];
		cp[i].xi_em[ii]		= SS_ref_db[ph_id].xi_em[ii];
		cp[i].mu[ii]		= SS_ref_db[ph_id].mu[ii];
	}
	for (int ii = 0; ii < gv.len_ox; ii++){
		cp[i].ss_comp[ii]	= SS_ref_db[ph_id].ss_comp[ii];
	}
	
	for (int ii = 0; ii < cp[i].n_sf; ii++){
		cp[i].sf[ii]		= SS_ref_db[ph_id].sf[ii];
	}
}



/**
	add minimized phase to LP PGE pseudocompound list 
*/
void copy_to_Ppc(		int 				 pc_check,
						int 				 add,
						int 				 ph_id,
						global_variable 	 gv,

						obj_type 			*SS_objective,
						SS_ref 			    *SS_ref_db					){

		double G;
		int    m_Ppc;
		int    save_mSS;


		if (pc_check != 2 || SS_ref_db[ph_id].df < gv.mSS_df_min_add || SS_ref_db[ph_id].df > gv.mSS_df_max_add){
			save_mSS = 0;
		}
		else{
			save_mSS = 1;
		}

		/* get unrotated gbase */
		SS_ref_db[ph_id] = non_rot_hyperplane(	gv, 
												SS_ref_db[ph_id]			);

		/* get unrotated minimized point informations */
		G 	=		 (*SS_objective[ph_id])(	SS_ref_db[ph_id].n_xeos,
												SS_ref_db[ph_id].iguess,
												NULL,
												&SS_ref_db[ph_id]			);

		/* check where to add the new phase PC */
		if (SS_ref_db[ph_id].id_Ppc >= SS_ref_db[ph_id].n_Ppc){ SS_ref_db[ph_id].id_Ppc = 0; printf("SS_LP, MAXIMUM STORAGE SPACE FOR PC IS REACHED for %4s, INCREASED #PC_MAX\n",gv.SS_list[ph_id]);}
		m_Ppc = SS_ref_db[ph_id].id_Ppc;

		if (save_mSS == 1){
			SS_ref_db[ph_id].info_Ppc[m_Ppc]   = 9;
		}
		else{
			SS_ref_db[ph_id].info_Ppc[m_Ppc]   = 0;
		}

		SS_ref_db[ph_id].DF_Ppc[m_Ppc]     = G;
		
		/* get pseudocompound composition */
		for (int j = 0; j < gv.len_ox; j++){				
			SS_ref_db[ph_id].comp_Ppc[m_Ppc][j] = SS_ref_db[ph_id].ss_comp[j]*SS_ref_db[ph_id].factor;	/** composition */
		}
		for (int j = 0; j < SS_ref_db[ph_id].n_em; j++){												/** save coordinates */
			SS_ref_db[ph_id].p_Ppc[m_Ppc][j]  = SS_ref_db[ph_id].p[j];												
			SS_ref_db[ph_id].mu_Ppc[m_Ppc][j] = SS_ref_db[ph_id].mu[j]*SS_ref_db[ph_id].z_em[j];										
		}
		/* save xeos */
		for (int j = 0; j < SS_ref_db[ph_id].n_xeos; j++){		
			SS_ref_db[ph_id].xeos_Ppc[m_Ppc][j] = SS_ref_db[ph_id].iguess[j];							/** compositional variables */
		}	
		SS_ref_db[ph_id].G_Ppc[m_Ppc] = G;
		
		/* add increment to the number of considered phases */
		if (SS_ref_db[ph_id].tot_Ppc < SS_ref_db[ph_id].n_Ppc){
			SS_ref_db[ph_id].tot_Ppc += 1;
		}
		else{
			SS_ref_db[ph_id].tot_Ppc = SS_ref_db[ph_id].n_Ppc;
		}
		SS_ref_db[ph_id].id_Ppc  += 1;
}

/** 
	Minimization function for PGE 
*/
void ss_min_PGE(		global_variable 	 gv,
						PC_type             *PC_read,

						obj_type 			*SS_objective,
						NLopt_type			*NLopt_opt,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id;
	int 	pc_check;
	clock_t u;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			pc_check = gv.PC_checked;
			ph_id = cp[i].id;
			cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */
			u = clock();
			/**
				set the iguess of the solution phase to the one of the considered phase
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k] = cp[i].xeos[k];
			}

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			/**
				Define a sub-hypervolume for the solution phases bounds
			*/
			SS_ref_db[ph_id] = restrict_SS_HyperVolume(	gv, 
														SS_ref_db[ph_id],
														gv.box_size_mode_PGE	);
			
			/**
				call to NLopt for non-linear + inequality constraints optimization
			*/
			SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
														SS_ref_db[ph_id]		);										
			
			/**
				establish a set of conditions to update initial guess for next round of local minimization 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k]   =  SS_ref_db[ph_id].xeos[k];
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

			u = clock() - u;
			SS_ref_db[ph_id].LM_time = ((double)u)/CLOCKS_PER_SEC*1000.0; 

			;
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

				// here we need to save the pseudocompound to have an estimate of the LP Matrix
				if (pc_check == 1){
					copy_to_Ppc(							pc_check,
															0,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db						);
				}					
			}
			else{
				if (gv.verbose == 1){
					printf(" !> SF [:%d] not respected for %4s (SS not updated)\n",SS_ref_db[ph_id].sf_id,gv.SS_list[ph_id]);
				}											
			}
		}
	}

};


/** 
	Minimization function for PGE 
*/
void init_PGE_from_LP(	global_variable 	 gv,
						PC_type				*PC_read,

						obj_type 			*SS_objective,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			ph_id = cp[i].id;

			/**
				set the iguess of the solution phase to the one of the considered phase 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k] = cp[i].xeos[k];
			}

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			
			SS_ref_db[ph_id] = PC_function(				gv,
														PC_read,
														SS_ref_db[ph_id], 
														z_b,
														ph_id 		);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
														SS_ref_db[ph_id], 
														z_b, 
														gv.SS_list[ph_id]		);

			copy_to_cp(									i, 
														ph_id,
														gv,
														SS_ref_db,
														cp						);	

		}
	}

};


/**
	Minimization function for PGE
*/
void ss_min_LP(			global_variable 	 gv,
						PC_type				*PC_read,

						obj_type 			*SS_objective,
						NLopt_type 			*NLopt_opt,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){

	double r;
	int 	ph_id;
	int     pc_check;
	int 	act;
	int 	candidate_ok, is_liq_synth_candidate;
	clock_t u;
	for (int i = 0; i < gv.len_ss; i++){
		gv.n_min[i] = 0;
	}
	int    is_gh            = (strcmp(gv.research_group, "gh") == 0);
	int    is_tc            = (strcmp(gv.research_group, "tc") == 0);

	int    n_liq_ph          = 0;
	int    ph_id_liq[MAX_LIQ_PHASES];
	int    liq_synth_active[MAX_LIQ_PHASES];
	int    liq_real_min_found[MAX_LIQ_PHASES];
	int    liq_candidate_index[MAX_LIQ_PHASES]; /* index of cp[] entry with the successful minimization */

	if (gv.liq_pc_synth_active && (is_gh || is_tc)){
		for (int iss = 0; iss < gv.len_ss && n_liq_ph < MAX_LIQ_PHASES; iss++){
			if (SS_ref_db[iss].is_liq == 1){ ph_id_liq[n_liq_ph++] = iss; }
		}
	}

	for (int l = 0; l < n_liq_ph; l++){
		int N_liq = 0;
		for (int i = 0; i < gv.len_ox; i++){
			if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id_liq[l]){ N_liq += 1; }
		}
		liq_synth_active[l]    = (N_liq >= gv.gh_liq_pc_synth_threshold);
		liq_real_min_found[l]  = 0;
		liq_candidate_index[l] = -1;
	}

	pc_check = gv.PC_checked;
	for (int i = 0; i < gv.len_cp; i++){

		if (cp[i].ss_flags[0] == 1){
			candidate_ok 	= 1;
			ph_id 			= cp[i].id;

			int liq_l = -1;
			for (int l = 0; l < n_liq_ph; l++){ if (ph_id_liq[l] == ph_id){ liq_l = l; break; } }

			/* NR-17/07/26
				Here I added a target rule for rMELTS. The liquid model including H2O and CO2 leads to a large number of local minimum. Brute force exploration is needed to lead the minimization toward the global minimum.
				This is something that is not needed for the other databases, as the Gibbs surface of other liquid models have smoother landscape.
			*/
			if (strcmp(gv.db, "rMELTS") == 0){
				if (gv.global_ite < gv.act_rMELTS_liq_pc_synth){
					act = 1;
				}
				else{
					is_liq_synth_candidate = (liq_l >= 0 && liq_synth_active[liq_l] && !liq_real_min_found[liq_l]);
					if (liq_l >= 0 && liq_synth_active[liq_l]){
						act = is_liq_synth_candidate ? 1 : 0;
					}
					else if ( SS_ref_db[ph_id].is_liq == 1 && gv.n_min[ph_id] > gv.n_max_val){
						act = 0;
					}
					else{
						act = 1;
					}
				}

			}
			else{
				// deactivating the next part helps for IGAD database at VHT
				is_liq_synth_candidate = (liq_l >= 0 && liq_synth_active[liq_l] && !liq_real_min_found[liq_l]);
				if (liq_l >= 0 && liq_synth_active[liq_l]){
					act = is_liq_synth_candidate ? 1 : 0;
				}
				else if ( SS_ref_db[ph_id].is_liq == 1 && gv.n_min[ph_id] > gv.n_max_val){
					act = 0;
				}
				else{
					act = 1;
				}
			}

			gv.n_min[ph_id] += 1;
			if (act == 1){
				cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */
				u = clock(); 
				/**
					set the iguess of the solution phase to the one of the considered phase 
				*/
				for (int k = 0; k < cp[i].n_xeos; k++) {
					SS_ref_db[ph_id].iguess[k] 	= cp[i].xeos[k];
					cp[i].xeos_0[k] 			= cp[i].xeos[k];
				}

				/**
					Rotate G-base hyperplane
				*/
				SS_ref_db[ph_id] = rotate_hyperplane(		gv, 
															SS_ref_db[ph_id]		);

				/**
					Define a sub-hypervolume for the solution phases bounds
				*/
				SS_ref_db[ph_id] = restrict_SS_HyperVolume(	gv, 
															SS_ref_db[ph_id],
															gv.box_size_mode_LP		);

				/**
					call to NLopt for non-linear + inequality constraints optimization
				*/
				SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
															SS_ref_db[ph_id]		);


				u = clock() - u;
				SS_ref_db[ph_id].LM_time = ((double)u)/CLOCKS_PER_SEC*1000.0; 

				if (gv.verbose == 1){

					printf(" NLopt status: %d\n", SS_ref_db[ph_id].status);
					SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																SS_ref_db[ph_id], 
																z_b, 
																gv.SS_list[ph_id]		);

					print_SS_informations(  				gv,
															SS_ref_db[ph_id],
															ph_id					);
				}
				if (strcmp(gv.research_group, "gh") == 0){
					if (SS_ref_db[ph_id].status != 3){
						double sum_x = 0.0;
						for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
							sum_x += SS_ref_db[ph_id].xeos[k];
						}
						if (fabs(sum_x - 1.0) > 1.0e-4 && sum_x > 0.0){
							for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
								SS_ref_db[ph_id].iguess[k] = SS_ref_db[ph_id].xeos[k] / sum_x;
							}
							SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
																		SS_ref_db[ph_id]		);
						}
					}
				}
				if (SS_ref_db[ph_id].status != 3){
					candidate_ok = 0;
					if (gv.verbose == 1){
						printf(" Failed minimization of %4s -> code %d\n",gv.SS_list[ph_id], SS_ref_db[ph_id].status);
					}
				}

				/**
					This next part is important to check if the minimized phase is valid and can be added to the PC list
				*/
				for (int k = 0; k < cp[i].n_xeos; k++) {
					cp[i].xeos_1[k] 			 =  SS_ref_db[ph_id].xeos[k];
				}
				for (int k = 0; k < cp[i].n_xeos; k++) {
					SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_1[k];
				}
				SS_ref_db[ph_id] = PC_function(				gv,
															PC_read,
															SS_ref_db[ph_id], 
															z_b,
															ph_id 		);
														
				SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
															SS_ref_db[ph_id],
															z_b,
															gv.SS_list[ph_id]		);
				/**
					add minimized phase to LP PGE pseudocompound list
				*/
				if (SS_ref_db[ph_id].sf_ok == 1){
					for (int k = 0; k < cp[i].n_xeos; k++) {
						SS_ref_db[ph_id].iguess[k]   =  SS_ref_db[ph_id].xeos[k];
					}
					SS_ref_db[ph_id] = PC_function(				gv,
																PC_read,
																SS_ref_db[ph_id], 
																z_b,
																ph_id 		);
															
					SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
																SS_ref_db[ph_id],
																z_b,
																gv.SS_list[ph_id]		);
					copy_to_Ppc(							pc_check,
															1,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db					);	
				}


				if (is_liq_synth_candidate && candidate_ok == 1){
					liq_real_min_found[liq_l] = 1;
					liq_candidate_index[liq_l] = i;
				}

			}

		}	
	}

	/**
	  	For each liquid phase present that has been minimized successfully,
		we can try to synthesize new PCs by linear interpolation between the real minimum and the other PCs of the same phase
	*/
	for (int l = 0; l < n_liq_ph; l++){
		if (!(liq_synth_active[l] && liq_real_min_found[l])){ continue; }

		int phid   = ph_id_liq[l];
		int n_xeos = SS_ref_db[phid].n_xeos;
		if (n_xeos > LIQ_PC_SYNTH_MAX_DIM){ continue; }

		double x_star[LIQ_PC_SYNTH_MAX_DIM];
		for (int k = 0; k < n_xeos; k++){ x_star[k] = SS_ref_db[phid].xeos[k]; }

		double steps[3] = {0.3, 0.6, 0.9};

		for (int ic = 0; ic < gv.len_ox; ic++){
			if (cp[ic].ss_flags[0] == 1 && cp[ic].id == phid && ic != liq_candidate_index[l]){
				for (int si = 0; si < 3; si++){
					double s = steps[si];
					int ok = 1;
					for (int k = 0; k < n_xeos; k++){
						double x_syn = x_star[k] + s * (cp[ic].xeos[k] - x_star[k]);
						if (x_syn < SS_ref_db[phid].bounds[k][0] || x_syn > SS_ref_db[phid].bounds[k][1]){ ok = 0; break; }
						SS_ref_db[phid].iguess[k] = x_syn;
					}
					if (!ok) { continue; }

					SS_ref_db[phid] = PC_function(gv, PC_read, SS_ref_db[phid], z_b, phid);
					SS_ref_db[phid] = SS_UPDATE_function(gv, SS_ref_db[phid], z_b, gv.SS_list[phid]);

					if (SS_ref_db[phid].sf_ok == 1){
						copy_to_Ppc(pc_check, 1, phid, gv, SS_objective, SS_ref_db);
						if (gv.verbose == 1){ printf(" [liq pc synth-linear] added PC for cp#%d step %g\n", ic, s); }
					}
					else{
						if (gv.verbose == 1){ printf(" [liq pc synth-linear] discarded PC (sf fail) for cp#%d step %g\n", ic, s); }
					}
				}
			}
		}
	}

};

/**
  initialize solution phase database
**/
global_variable init_ss_db(		int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){


	if (EM_database == 0){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mp_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 1){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mb_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 11){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mb_ext_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 2){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_ig_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 22){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_igd_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 3){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_igad_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 4 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_um_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 5 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_um_ext_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 6 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mtl_EM_function(		gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 7 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mpe_EM_function(		gv,
														SS_ref_db[i],
														gv.EM_dataset,
														z_b,
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 8 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;
			SS_ref_db[i].R  = 0.0083144;

			SS_ref_db[i]    = G_SS_all_EM_function(		gv,
														SS_ref_db[i],
														gv.EM_dataset,
														z_b,
														gv.SS_list[i]		);
		}
	}

	return gv;
};



/**
  initialize solution phase database
**/
global_variable init_ss_db_sb(	int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){

	if (EM_database == 0){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb11_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
											
										/** can become a global variable instead */
		}
	}
	else if (EM_database == 1){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb21_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
											
										/** can become a global variable instead */
		}
	}
	else if (EM_database == 2){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb24_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset,
														z_b,
														gv.SS_list[i]		);

										/** can become a global variable instead */
		}
	}
	return gv;
};

/**
  initialize solution phase database for the "gh" (Ghiorso/MELTS) research group
**/
global_variable init_ss_db_gh(	int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){
	/* see init_em_db_gh's own comment on why this is also set here, not
	   just in GH_SS_objective_init_function. */
	GH_actual_EM_database = gv.EM_database;
	for (int i = 0; i < gv.len_ss; i++){
		SS_ref_db[i].P  = z_b.P;
		SS_ref_db[i].T  = z_b.T;
		SS_ref_db[i].R  = 0.0083144;

		SS_ref_db[i]    = G_SS_gh_EM_function(	gv,
												SS_ref_db[i],
												gv.EM_dataset,
												z_b,
												gv.SS_list[i]		);
	}
	return gv;
};
