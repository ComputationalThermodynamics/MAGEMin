/**
Function to call solution phase Minimization        
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "mpi.h"
#include "nlopt.h"
#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "NLopt_opt_function.h"
#include "dump_function.h"
#include "toolkit.h"
#include "phase_update_function.h"
#include "objective_functions.h"

/** 
Function to update xi and sum_xi during local minimization.
*/
SS_ref SS_UPDATE_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								bulk_info 	 		 z_b,
								char    			*name){

	/* sf_ok?*/
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
		SS_ref_db.xi_em[i] = exp(-SS_ref_db.mu[i]/(SS_ref_db.R*SS_ref_db.T));
		SS_ref_db.sum_xi  += SS_ref_db.xi_em[i]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
	}

	/* get composition of solution phase */
	for (int j = 0; j < gv.len_ox; j++){
		SS_ref_db.ss_comp[j] = 0.0;
		for (int i = 0; i < SS_ref_db.n_em; i++){
		   SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
	   } 
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
		cp.xi_em[i] = exp(-cp.mu[i]/(SS_ref_db.R*SS_ref_db.T));
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
int copy_to_Ppc_composite(		int 				 ph_id,
								global_variable 	 gv,

								obj_type 			*SS_objective,
								SS_ref 			    *SS_ref_db				){

		double G;
		int    m_Ppc;

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

		SS_ref_db[ph_id].info_Ppc[m_Ppc]   = 0;

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
		SS_ref_db[ph_id].tot_Ppc += 1;
		SS_ref_db[ph_id].id_Ppc  += 1;

		return m_Ppc;
}



/**
	add minimized phase to LP PGE pseudocompound list 
*/
void copy_to_Ppc(		int 				 pc_check,
						int 				 add,
						int 				 ph_id,
						global_variable 	 gv,

						obj_type 			*SS_objective,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp					){

		double G;
		int    m_Ppc;

		if (add != 0 || SS_ref_db[ph_id].df_raw < 1e-3 || SS_ref_db[ph_id].df_raw > 0.25){
			pc_check = 0;
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

		if (pc_check == 1){
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
		SS_ref_db[ph_id].tot_Ppc += 1;
		SS_ref_db[ph_id].id_Ppc  += 1;
}

/** 
	Minimization function for PGE 
*/
void ss_min_PGE(		global_variable 	 gv,

						obj_type 			*SS_objective,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id;
	int 	pc_check;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			pc_check = gv.PC_checked;
			ph_id = cp[i].id;
			cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */

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
			SS_ref_db[ph_id] = NLopt_opt_function(		gv, 
														SS_ref_db[ph_id], 
														ph_id					);
			
			/**
				establish a set of conditions to update initial guess for next round of local minimization 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k]   =  SS_ref_db[ph_id].xeos[k];
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

				// here we need to save the pseudocompound to have an estimate of the LP Matrix
				if (pc_check == 1){
					copy_to_Ppc(							pc_check,
															0,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db,
															cp						);
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
														SS_ref_db[ph_id], 
														z_b,
														gv.SS_list[ph_id] 		);
													
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


void compute_cst_dG_Ppc(	global_variable 	 gv,
							obj_type 			*SS_objective,
							bulk_info 	 		 z_b,
							SS_ref 			    *SS_ref_db,
							csd_phase_set  		*cp,

							int					 ph_id,
							int					 cp_id
){
	double 	delta_G, a, b, c;
	int		n, conv, sign_a, sign_c;
	int    	i, j, k;

	double 	tol        	= 1e-6;						// tolerance on delta_G
	double 	target_dg   = 1e-4;						// delta_G for generated set of Ppc
	double 	n_max       = 8;						// maximum number of iterations
	double  ref_df 		= SS_ref_db[ph_id].df;

	int    	n_xeos 		= SS_ref_db[ph_id].n_xeos;
	int    	n_em 		= SS_ref_db[ph_id].n_em;

	for (i = 0; i < n_em; i++){
		for (k = 0; k < cp[cp_id].n_xeos; k++) {
			cp[i].xeos_r[k] = (rnd(1.0) -0.5) / 100.0;
		}

        delta_G     = 1.0;                                                  // initialize missfit
        a           = 0.0;
        b           = 1.0;
        conv        = 0;
        n           = 0;
        sign_a      = -1;

		while (n < n_max && conv == 0){

			c = (a+b)/2.0;

			for (k = 0; k < cp[cp_id].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k]   =  cp[cp_id].xeos_1[k] + cp[cp_id].xeos_r[k]*c;
			}

			SS_ref_db[ph_id] = PC_function(				gv,
														SS_ref_db[ph_id], 
														z_b,
														gv.SS_list[ph_id] 		);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
														SS_ref_db[ph_id], 
														z_b, 
														gv.SS_list[ph_id]		);

			delta_G = SS_ref_db[ph_id].df - ref_df - target_dg;

			// printf(" df: %+10f delta_G: %+10f\n",SS_ref_db[ph_id].df,delta_G);
			if (fabs(delta_G) < tol){
				conv = 1;
			}
			else{
				sign_c = delta_G/fabs(delta_G);
				if (sign_c == sign_a){
					a = c;
					sign_a = sign_c;
				}
				else{
					b = c;
				}
			}

			n += 1;
		}

		if (SS_ref_db[ph_id].sf_ok == 1){
			copy_to_Ppc(								0,
														1,
														ph_id,
														gv,

														SS_objective,
														SS_ref_db,
														cp						);	
		}

	}

}


/** 
	Minimization function for PGE 
*/
void ss_min_LP(			global_variable 	 gv,

						obj_type 			*SS_objective,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){

	double r;
	int 	ph_id;
	int     pc_check;
	int 	act;

	for (int i = 0; i < gv.len_ss; i++){ 
		gv.n_min[i] = 0;
	}

	for (int i = 0; i < gv.len_cp; i++){ 
		pc_check = gv.PC_checked;

		if (cp[i].ss_flags[0] == 1){
			ph_id = cp[i].id;

			if ( strcmp( gv.SS_list[ph_id], "liq") == 0 && gv.n_min[ph_id] > 2){
				act = 0;
			}
			else{
				act = 1;
			}
			gv.n_min[ph_id] += 1;

			if (act == 1){
				cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */

				/**
					set the iguess of the solution phase to the one of the considered phase 
				*/
				for (int k = 0; k < cp[i].n_xeos; k++) {
					SS_ref_db[ph_id].iguess[k] 	= cp[i].xeos[k];
					cp[i].xeos_0[k] 			= cp[i].xeos[k];;
					// SS_ref_db[ph_id].dguess[k] = cp[i].xeos[k];			//dguess can be used of LP, it is used for PGE to check for drifting
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
				SS_ref_db[ph_id] = NLopt_opt_function(		gv, 
															SS_ref_db[ph_id], 
															ph_id					);

				/** 
					print solution phase informations (print has to occur before saving PC)
				*/
				if (gv.verbose == 1){
					SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																SS_ref_db[ph_id], 
																z_b, 
																gv.SS_list[ph_id]		);

					print_SS_informations(  				gv,
															SS_ref_db[ph_id],
															ph_id					);
				}


				for (int k = 0; k < cp[i].n_xeos; k++) {
					cp[i].xeos_1[k] 			 =  SS_ref_db[ph_id].xeos[k];
				}
				

				// compute_cst_dG_Ppc(	gv,
				// 					SS_objective,
				// 					z_b,
				// 					SS_ref_db,
				// 					cp,

				// 					ph_id,
				// 					i
				// );

				double shift = 0.0;
				double sh_array[] = {0.0,-0.0001,0.0001,0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.75};

				int add_def = 0;
				for (int add = 0; add < 11; add++){
					
					shift = sh_array[add];
					for (int k = 0; k < cp[i].n_xeos; k++) {
						SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_1[k] * (1.0-shift) + cp[i].xeos_0[k] * (shift);
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
						add minimized phase to LP PGE pseudocompound list 
					*/
					if (SS_ref_db[ph_id].sf_ok == 1){
						copy_to_Ppc(							pc_check,
																add,
																ph_id,
																gv,

																SS_objective,
																SS_ref_db,
																cp						);	
					}
					else{
						if (add_def == 0){
							for (int k = 0; k < cp[i].n_xeos; k++) {
								SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_0[k];
							}
							
							SS_ref_db[ph_id] = PC_function(				gv,
																		SS_ref_db[ph_id], 
																		z_b,
																		gv.SS_list[ph_id] 		);
																	
							SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																		SS_ref_db[ph_id], 
																		z_b, 
																		gv.SS_list[ph_id]		);

							copy_to_Ppc(								0,
																		1,
																		ph_id,
																		gv,

																		SS_objective,
																		SS_ref_db,
																		cp						);	
							add_def = 1;
						}

						// if (gv.verbose == 1){
						// 	printf(" !> SF [:%d] not respected for %4s (SS not updated)\n",SS_ref_db[ph_id].sf_id,gv.SS_list[ph_id]);
						// }											
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
													EM_database, 
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

			SS_ref_db[i]    = G_SS_mb_EM_function(	gv, 
													SS_ref_db[i], 
													EM_database, 
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

			SS_ref_db[i]    = G_SS_ig_EM_function(	gv, 
													SS_ref_db[i], 
													EM_database, 
													z_b, 
													gv.SS_list[i]		);
											
										/** can become a global variable instead */
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
													EM_database, 
													z_b, 
													gv.SS_list[i]		);
											
										/** can become a global variable instead */
		}
	}
	return gv;
};