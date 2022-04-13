/**
Function to call solution phase Minimization        
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include <lapacke.h> 

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
SS_ref SS_UPDATE_function(		global_variable gv,
								SS_ref SS_ref_db, 
								struct bulk_info z_b,
								char    *name){

	/* sf_ok?*/
	SS_ref_db.sf_ok = 1;
	for (int i = 0; i < SS_ref_db.n_sf; i++){
		if (SS_ref_db.sf[i] <= 0.0 || isnan(SS_ref_db.sf[i]) == 1|| isinf(SS_ref_db.sf[i]) == 1){
			SS_ref_db.sf_ok = 0;	
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
	for (int j = 0; j < nEl; j++){
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
										struct bulk_info 	z_b			){

	/* sf_ok?*/
	cp.sf_ok = 1;
	for (int i = 0; i < cp.n_sf; i++){
		if (cp.sf[i] <= 0.0 || isnan(cp.sf[i]) == 1|| isinf(cp.sf[i]) == 1){
			cp.sf_ok = 0;	
			break;
		}
	}

	/* xi calculation (phase fraction expression for PGE) */
	cp.sum_xi 	= 0.0;	
	for (int i = 0; i < cp.n_em; i++){ 
		cp.xi_em[i] = exp(-cp.mu[i]/(SS_ref_db.R*SS_ref_db.T));
		cp.sum_xi  += cp.xi_em[i]*cp.p_em[i]*SS_ref_db.z_em[i];
	}

	/* get composition of solution phase */
	for (int j = 0; j < nEl; j++){
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
global_variable split_cp(		int 				 i, 
								global_variable 	 gv,
								SS_ref 			    *SS_ref_db,
								csd_phase_set  		*cp
){
	int id_cp;
	int ph_id = cp[i].id;
	double distance;

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
		
		gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
		gv.n_solvi[ph_id] 	+= 1;
		gv.len_cp 			+= 1;
		
		if (gv.verbose == 1){
			printf("\n  {FYI} %4s cp#%d is grazing away from its field, a copy has been added (xeos = dguess)\n",gv.SS_list[ph_id],i);
		}
		
		if (gv.len_cp == gv.max_n_cp){
			printf(" !! Maxmimum number of allowed phases under consideration reached !!\n    -> check your problem and potentially increase gv.max_n_cp\n");
		}
		
	}
	
	return gv;
};

/** 
	Minimization function for PGE 
*/
void ss_min_PGE(		int 				mode, 
						int 				i, 
						global_variable 	gv,
						struct 	bulk_info 	z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id = cp[i].id;
	
	cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */
	SS_ref_db[ph_id].min_mode 	= mode;								/** send the right mode to the local minimizer */
	gv.maxeval   		  		= gv.maxeval_mode_1;

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

	double relax, norm;
	if (cp[i].in_iter > 0 && abs(cp[i].in_iter - gv.global_ite) <= 8){
		relax = 1.0;
	}
	else{
		relax = 64.0;
	}

	norm = 1.0;
	if (gv.BR_norm < gv.relax_PGE){
		norm = norm_vector(cp[i].mu,cp[i].n_em)/relax;
		if (norm > 1.0){ 
			norm = 1.0;
		}
	}

	/**
		Define a sub-hypervolume for the solution phases bounds
	*/
	SS_ref_db[ph_id] = restrict_SS_HyperVolume(	gv, 
												SS_ref_db[ph_id],
												gv.box_size_mode_1*norm	);
	
	/**
		call to NLopt for non-linear + inequality constraints optimization
	*/
	SS_ref_db[ph_id] = NLopt_opt_function(	gv, 
											SS_ref_db[ph_id], 
											ph_id						);
	
	/**
		establish a set of conditions to update initial guess for next round of local minimization 
	*/
	for (int k = 0; k < cp[i].n_xeos; k++) {
		SS_ref_db[ph_id].iguess[k]  = SS_ref_db[ph_id].xeos[k];
	}	

	SS_ref_db[ph_id] = PC_function(			gv,
											SS_ref_db[ph_id], 
											z_b,
											gv.SS_list[ph_id] 			);
											
	SS_ref_db[ph_id] = SS_UPDATE_function(	gv, 
											SS_ref_db[ph_id], 
											z_b, 
											gv.SS_list[ph_id]			);
	
	/**
		copy the minimized phase informations to cp structure, if the site fractions are respected
	*/
	if (SS_ref_db[ph_id].sf_ok == 1){
		cp[i].min_time			= SS_ref_db[ph_id].LM_time;
		cp[i].df				= SS_ref_db[ph_id].df_raw;
		cp[i].factor			= SS_ref_db[ph_id].factor;
		cp[i].sum_xi			= SS_ref_db[ph_id].sum_xi;

		for (int ii = 0; ii < cp[i].n_xeos; ii++){
			cp[i].xeos[ii]		= SS_ref_db[ph_id].iguess[ii]; 
			cp[i].dfx[ii]		= SS_ref_db[ph_id].dfx[ii]; 
		}
		
		for (int ii = 0; ii < cp[i].n_em; ii++){
			cp[i].p_em[ii]		= SS_ref_db[ph_id].p[ii];
			cp[i].xi_em[ii]		= SS_ref_db[ph_id].xi_em[ii];
			cp[i].mu[ii]		= SS_ref_db[ph_id].mu[ii];
		}
		for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
			for (int jj = 0; jj < SS_ref_db[ph_id].n_xeos; jj++){
				cp[i].dpdx[ii][jj] = SS_ref_db[ph_id].dp_dx[ii][jj];
			}
		}
		for (int ii = 0; ii < gv.len_ox; ii++){
			cp[i].ss_comp[ii]	= SS_ref_db[ph_id].ss_comp[ii];
		}
		
		for (int ii = 0; ii < cp[i].n_sf; ii++){
			cp[i].sf[ii]		= SS_ref_db[ph_id].sf[ii];
		}	
	}
	else{
		if (gv.verbose == 1){
			printf(" !> SF not respected for %4s (SS not updated)\n",gv.SS_list[ph_id]);
		}	
	}

	/** 
		print solution phase informations 
	*/
	if (gv.verbose == 1){
		print_SS_informations(  gv,
								SS_ref_db[ph_id],
								ph_id									);
	}
};


/**
  reset global variable for parallel calculations 
*/
global_variable reset_gv(					global_variable 	 gv,
											struct bulk_info 	 z_b,
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db
){
	int i,k;
	for (k = 0; k < gv.n_flags; k++){
		for (i = 0; i < gv.len_pp; i++){
			gv.pp_flags[i][k]   = 0;
		}
		for (int i = 0; i < gv.len_ss; i++){ 
			SS_ref_db[i].ss_flags[k]   = 0;
		}
	}
	
	/* reset pure phases fractions and xi */
	for (int i = 0; i < gv.len_pp; i++){		
		gv.pp_n[i] 		  = 0.0;
		gv.delta_pp_n[i]  = 0.0;
		gv.pp_xi[i] 	  = 0.0;
		gv.delta_pp_xi[i] = 0.0;
	}
	
	/* reset pure phases */
	char liq_tail[] = "L";
	for (int i = 0; i < gv.len_pp; i++){
		if ( EndsWithTail(gv.PP_list[i], liq_tail) == 1 ){
			if (z_b.T < 773.0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}
			else{
				gv.pp_flags[i][0] = 1;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 0;
			}
		}
		else{
			gv.pp_flags[i][0] = 1;
			gv.pp_flags[i][1] = 0;
			gv.pp_flags[i][2] = 1;
			gv.pp_flags[i][3] = 0;
		}
	}
	/* reset solution phases frations */
	for (int i = 0; i < gv.len_ss; i++){ 	
		SS_ref_db[i].ss_n = 0.0;
	}
	
	gv.system_density     = 0.;
	gv.system_bulkModulus = 0.;
	gv.system_shearModulus= 0.;
	gv.system_Vp 		  = 0.;
	gv.system_Vs 		  = 0.;

	gv.check_PC_ite		  = 0;
	gv.check_PC			  = 0;
	gv.maxeval		      = gv.maxeval_mode_1;
	gv.len_cp 		  	  = 0;
	gv.div				  = 0;
	gv.status			  = 0;
	gv.n_phase            = 0;                  /** reset the number of phases to start with */
	gv.global_ite		  = 0;					/** reset the number of global iteration to zero */
	gv.G_system           = 0.0;
	gv.n_cp_phase         = 0;					/** reset the number of ss phases to start with */
	gv.n_pp_phase         = 0;					/** reset the number of pp phases to start with */
	gv.alpha          	  = gv.max_fac;

    for (i = 0; i < gv.it_f; i++){
        gv.PGE_mass_norm[i] = 0.0;
    }

    for (i = 0; i < gv.len_ox; i++){	
        gv.mass_residual[i] = 0.0;
        gv.gam_tot[i]     	= 0.0;
        gv.del_gam_tot[i]   = 0.0;
        gv.delta_gam_tot[i] = 0.0;
    }

    for (i = 0; i < gv.len_ss; i++){	
        gv.n_solvi[i] = 0;
		for (k = 0; k < gv.max_n_cp; k++){	
			gv.id_solvi[i][k] = 0;
		} 
    }

	return gv;
}

/**
  initilize solution phase database
**/
global_variable init_ss_db(		int 				 EM_database,
								struct bulk_info 	 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){
	double R  = 0.0083144; 
	for (int i = 0; i < gv.len_ss; i++){
	
		SS_ref_db[i]    = G_SS_EM_function(		gv, 
												SS_ref_db[i], 
												EM_database, 
												z_b, 
												gv.SS_list[i]		);
										 
		SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
		SS_ref_db[i].T  = z_b.T;		
		SS_ref_db[i].R  = R;										/** can become a global variable instead */

	}

	return gv;
};

/**
  reset stable phases entries
*/
void reset_sp(						global_variable 	 gv,
									stb_system  		*sp
){
	/* reset system */
	for (int i = 0; i < gv.len_ox; i++){
		strcpy(sp[0].ph[i],"");	
		sp[0].bulk[i] 					= 0.0;
		sp[0].gamma[i] 					= 0.0;
		sp[0].bulk_S[i] 				= 0.0;
		sp[0].bulk_M[i] 				= 0.0;
		sp[0].bulk_F[i] 				= 0.0;

		sp[0].ph_type[i] 				= -1;
		sp[0].ph_id[i] 					=  0;
		sp[0].ph_frac[i] 				=  0.0;
	}

	/* reset phases */
	for (int n = 0; n < gv.len_ox; n++){
		for (int i = 0; i < gv.len_ox; i++){
			sp[0].PP[n].Comp[i] 			= 0.0;
			sp[0].SS[n].Comp[i] 			= 0.0;
			sp[0].SS[n].compVariables[i] 	= 0.0;
		}
		for (int i = 0; i < gv.len_ox+1; i++){
			strcpy(sp[0].SS[n].emNames[i],"");	
			
			sp[0].SS[n].emFrac[i] 			= 0.0;
			sp[0].SS[n].emChemPot[i] 		= 0.0;

			for (int j = 0; j < gv.len_ox; j++){
				sp[0].SS[n].emComp[i][j]	= 0.0;
			}
		}
	}

}

/**
  reset considered phases entries
*/
void reset_cp(						global_variable 	 gv,
									struct bulk_info 	 z_b,
									csd_phase_set  		*cp
){
	
	for (int i = 0; i < gv.max_n_cp; i++){		
		strcpy(cp[i].name,"");					/* get phase name */	
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
		cp[i].delta_ss_n    	= 0.0;				/* get initial phase fraction */
		
		for (int ii = 0; ii < gv.len_ox + 1; ii++){
			cp[i].p_em[ii]      = 0.0;
			cp[i].xi_em[ii]     = 0.0;
			cp[i].dguess[ii]    = 0.0;
			cp[i].xeos[ii]      = 0.0;
			cp[i].dfx[ii]       = 0.0;
			cp[i].mu[ii]        = 0.0;
			cp[i].mu0[ii]       = 0.0;
			cp[i].ss_comp[ii]   = 0.0;
		}
		 
		for (int ii = 0; ii < (gv.len_ox + 1)*2; ii++){
			cp[i].sf[ii]    	= 0.0;
		}
		cp[i].mass 				= 0.0;
		cp[i].volume 			= 0.0;
		cp[i].phase_density 	= 0.0;
		cp[i].phase_cp 			= 0.0;
	}

};

/**
  reset compositional variables (xeos) when something goes wrong during minimization 
*/
void reset_SS(						global_variable 	 gv,
									struct bulk_info 	 z_b,
									SS_ref 				*SS_ref_db
){
	/* reset solution phases */
	for (int iss = 0; iss < gv.len_ss; iss++){
		SS_ref_db[iss].min_mode	= 1;
		SS_ref_db[iss].tot_pc 	= 0;
		SS_ref_db[iss].id_pc  	= 0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db[iss].solvus_id[j] = -1;	
		}
		for (int i = 0; i < (SS_ref_db[iss].n_pc); i++){
			SS_ref_db[iss].factor_pc[i] = 0.0;
			SS_ref_db[iss].n_swap[i] = 0;
			SS_ref_db[iss].G_pc[i]   = 0.0;
			SS_ref_db[iss].DF_pc[i]  = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				SS_ref_db[iss].comp_pc[i][j]  = 0.0;
			}
			for (int j = 0; j < SS_ref_db[iss].n_em; j++){
				SS_ref_db[iss].p_pc[i][j]  = 0;	
				SS_ref_db[iss].mu_pc[i][j] = 0;	
			}
			for (int j = 0; j < (SS_ref_db[iss].n_xeos); j++){
				SS_ref_db[iss].xeos_pc[i][j]  = 0.0;
			}
		}
		for (int j = 0; j < SS_ref_db[iss].n_em; j++){
			SS_ref_db[iss].xi_em[j]     = 0.0;
			SS_ref_db[iss].z_em[j]      = 1.0;
			SS_ref_db[iss].mu[j] 	    = 0.0;
		}
		SS_ref_db[iss].sum_xi		    = 0.0;
		SS_ref_db[iss].df			    = 0.0;
		SS_ref_db[iss].df_raw		    = 0.0;

		for (int k = 0; k < SS_ref_db[iss].n_xeos; k++) {					/** initialize initial guess using default thermocalc guess	*/ 
			SS_ref_db[iss].iguess[k]  = 0.0;
			SS_ref_db[iss].dguess[k]  = 0.0;
			SS_ref_db[iss].xeos[k]    = 0.0;
			SS_ref_db[iss].bounds[k][0] = SS_ref_db[iss].bounds_ref[k][0];
			SS_ref_db[iss].bounds[k][1] = SS_ref_db[iss].bounds_ref[k][1];
		}		
	}	

};


/**
Function to destroy allocated memory of solution phase structure
*/
void SS_ref_destroy(	global_variable gv, 
						SS_ref *SS_ref_db		){	
											  
	for (int i = 0; i < gv.len_ss; i++){
		
		free(SS_ref_db[i].ss_flags);
		for (int j = 0; j < SS_ref_db[i].n_em; j++) {
			free(SS_ref_db[i].dp_dx[j]);
			free(SS_ref_db[i].Comp[j]);
		}
		free(SS_ref_db[i].dp_dx);
		free(SS_ref_db[i].Comp);
		
		free(SS_ref_db[i].gbase);
		free(SS_ref_db[i].gb_lvl);
		free(SS_ref_db[i].z_em);
		free(SS_ref_db[i].density);
		free(SS_ref_db[i].dguess);
		free(SS_ref_db[i].iguess);
		free(SS_ref_db[i].p);
		free(SS_ref_db[i].mat_phi);
		free(SS_ref_db[i].mu_Gex);
		free(SS_ref_db[i].sf);
		free(SS_ref_db[i].mu);
		free(SS_ref_db[i].dfx);
		free(SS_ref_db[i].ss_comp);
		free(SS_ref_db[i].xi_em);	
		free(SS_ref_db[i].xeos);
		free(SS_ref_db[i].dsf);

		/** destroy box bounds */
		for (int j = 0; j< SS_ref_db[i].n_xeos; j++) {
			free(SS_ref_db[i].bounds[j]);
			free(SS_ref_db[i].bounds_ref[j]);
		}
		free(SS_ref_db[i].bounds_ref);			
		free(SS_ref_db[i].bounds);
		
		/** free pseudocompound related memory */
		for (int j = 0; j< SS_ref_db[i].n_pc; j++) {
			free(SS_ref_db[i].comp_pc[j]);
			free(SS_ref_db[i].p_pc[j]);
			free(SS_ref_db[i].xeos_pc[j]);
		}
		free(SS_ref_db[i].comp_pc);
		free(SS_ref_db[i].n_swap);
		free(SS_ref_db[i].info);
		free(SS_ref_db[i].xeos_pc);
		free(SS_ref_db[i].p_pc);
		free(SS_ref_db[i].G_pc);
		free(SS_ref_db[i].factor_pc);
		free(SS_ref_db[i].DF_pc);
		free(SS_ref_db[i].ub_pc);
		free(SS_ref_db[i].lb_pc);
		free(SS_ref_db[i].xeos_sf_ok);
	}
};

/**
	Function to destroy allocated memory of considered phase structure
*/
void CP_destroy(		global_variable gv, 
						csd_phase_set  *cp		){	
											  
	for (int i = 0; i < gv.max_n_cp; i++){
		free(cp[i].name);
		free(cp[i].p_em);
		free(cp[i].xi_em);
		free(cp[i].dguess);
		free(cp[i].xeos);
		free(cp[i].ss_flags);
		free(cp[i].ss_comp);
		free(cp[i].dfx);
		free(cp[i].sf);
		free(cp[i].mu);
	}
};
