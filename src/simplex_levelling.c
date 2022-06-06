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
#include <lapacke.h> 
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
#include "SS_xeos_PC.h"


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
	associate the array of pointer with the right solution phase
*/
void SS_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv				){	
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "bi")  == 0 ){
			SS_objective[iss]  = obj_bi; 		}
		else if (strcmp( gv.SS_list[iss], "cd")  == 0){
			SS_objective[iss]  = obj_cd; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_objective[iss]  = obj_cpx; 		}
		else if (strcmp( gv.SS_list[iss], "ep")  == 0){
			SS_objective[iss]  = obj_ep; 		}
		else if (strcmp( gv.SS_list[iss], "fl")  == 0){
			SS_objective[iss]  = obj_fl; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_objective[iss]  = obj_g; 		}
		else if (strcmp( gv.SS_list[iss], "hb")  == 0){
			SS_objective[iss]  = obj_hb; 		}
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			SS_objective[iss]  = obj_ilm; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			SS_objective[iss]  = obj_liq; 		}
		else if (strcmp( gv.SS_list[iss], "mu")  == 0){
			SS_objective[iss]  = obj_mu; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_objective[iss]  = obj_ol; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_objective[iss]  = obj_opx; 		}
		else if (strcmp( gv.SS_list[iss], "pl4T") == 0){
			SS_objective[iss]  = obj_pl4T; 		}
		else if (strcmp( gv.SS_list[iss], "spn") == 0){
			SS_objective[iss]  = obj_spn; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};			
}

/**
  update global Gamma 
*/	
void update_global_gamma( 				struct bulk_info 	z_b,
										simplex_data 	   *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

	/** update gam_tot using solution phase levelling gamma */
	for (int i = 0; i < d->n_Ox; i++){
		d->gamma_delta[z_b.nzEl_array[i]] 	= d->gamma_ss[i] - d->gamma_tot[z_b.nzEl_array[i]];
		d->gamma_tot[z_b.nzEl_array[i]]     	= d->gamma_ss[i];
	}	
	
};

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
	if (d->dG_B < -1e-6){	
		d->min_F  =  1e6;												/** max value for F, tentative here because F can tend to +inf */
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
  function to run simplex linear programming 
*/	
void fill_simplex_arrays_A(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){

	simplex_data *d  = (simplex_data *) splx_data;
	/* fill reference assemblage */

	for (int k = 0; k < z_b.nzEl_val; k++) {
		d->g0_A[k]    		   = 1e10;								/** penalty G */
		d->ph_id_A[k][0]  	   = 0;									/** phase_id for penalty phase */
		d->A[k+k*z_b.nzEl_val] = 1.0;								/** eye matrix for stoechiometry */
		d->n_vec[k] 		   = z_b.bulk_rock[z_b.nzEl_array[k]];	/** initial phase fraction simply corresponds to the bulk rock composition Ax = br */
	}

}

/**
	print levelling informations 
*/
void print_levelling(		struct bulk_info 	 z_b,
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
		max_n_pc = get_max_n_pc(SS_ref_db[i].tot_pc, SS_ref_db[i].n_pc);	

		for (int l = 0; l < max_n_pc; l++){
			SS_ref_db[i].DF_pc[l] = SS_ref_db[i].G_pc[l];
			for (int j = 0; j < gv.len_ox; j++) {
				SS_ref_db[i].DF_pc[l] -= SS_ref_db[i].comp_pc[l][j]*gv.gam_tot[j];
			}
			
			if (SS_ref_db[i].DF_pc[l] < 1.0){
				printf(" %4s %04d  #swap: %04d #stage %04d | ",gv.SS_list[i],l,SS_ref_db[i].n_swap[l],SS_ref_db[i].info[l]);
				printf("DF: %+4f | ", SS_ref_db[i].DF_pc[l]);
				//for (int j = 0; j < SS_ref_db[i].n_xeos; j++){
					//printf(" %+4f",SS_ref_db[i].xeos_pc[l][j]);
				//}
				//for (int k = SS_ref_db[i].n_xeos; k < 10; k++){
					//printf(" %4s","-");
				//}
				//printf(" | ");

				for (int j = 0; j < SS_ref_db[i].n_em; j++){
					//for (int k = 0; k < gv.len_ox; k++) {
						//SS_ref_db[i].mu_pc[l][j] -= SS_ref_db[i].comp_pc[l][k]*gv.gam_tot[k];
					//}
					printf(" %+4f",SS_ref_db[i].mu_pc[l][j]);
				}
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

			}
		}
	}
}

/**
	Generate pseudocompounds
*/
void generate_pseudocompounds(	int 		 		 ss,
								struct bulk_info 	 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db,
								PC_ref 				*SS_PC_xeos,
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
		get_ss_pv = SS_PC_xeos[ss].ss_pc_xeos[k]; 
		
		/* TMP, not so elegant way to deal with cases were an oxide of the bulk rock composition = 0.0 */	
		for (int i = 0; i < SS_ref_db[ss].n_xeos; i++){
			if (get_ss_pv.xeos_pc[i] > SS_ref_db[ss].bounds_ref[i][1]){
				get_ss_pv.xeos_pc[i] = SS_ref_db[ss].bounds_ref[i][1];
			}
		}
		G 	= (*SS_objective[ss])(SS_ref_db[ss].n_xeos, get_ss_pv.xeos_pc, 	NULL, &SS_ref_db[ss]);

		/** store pseudocompound */
		if ( G < gv.max_G_pc ){

			swp     = 0;
			while (swp == 0){

				/* get composition of solution phase */
				for (j = 0; j < nEl; j++){
					SS_ref_db[ss].ss_comp[j] = 0.0;
					for (i = 0; i < SS_ref_db[ss].n_em; i++){
					SS_ref_db[ss].ss_comp[j] += SS_ref_db[ss].Comp[i][j]*SS_ref_db[ss].p[i]*SS_ref_db[ss].z_em[i];
				} 
				}

				if (SS_ref_db[ss].id_pc >= SS_ref_db[ss].n_pc){ SS_ref_db[ss].id_pc = 0; printf("MAXIMUM STORAGE SPACE FOR PC IS REACHED, INCREASED #PC_MAX\n");}
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
					SS_ref_db[ss].mu_pc[m_pc][j] = SS_ref_db[ss].mu[j]*SS_ref_db[ss].z_em[j];										
				}
				/* save xeos */
				for ( j = 0; j < SS_ref_db[ss].n_xeos; j++){		
					SS_ref_db[ss].xeos_pc[m_pc][j] = get_ss_pv.xeos_pc[j];							/** compositional variables */
				}	
				
				SS_ref_db[ss].G_pc[m_pc] = G;
				
				/* add increment to the number of considered phases */
				SS_ref_db[ss].tot_pc += 1;
				SS_ref_db[ss].id_pc  += 1;
				swp = 1;	
			}	
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
			
			max_n_pc = get_max_n_pc(SS_ref_db[iss].tot_pc, SS_ref_db[iss].n_pc);
			
			for (int l = 0; l < max_n_pc; l++){
				/* if the driving force of the pseudocompound is lower to the filter then consider it */
				if (SS_ref_db[iss].DF_pc[l]*SS_ref_db[iss].factor_pc[l] < gv.bnd_filter_pc){
					phase_on = 1;
				}	
			}
			/* if no generated pseudocompound are close to the hyperplane then turn off the solution phase */
			if (phase_on == 0){
				if (gv.verbose != -1){
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

global_variable update_global_info(		struct bulk_info 	 z_b,
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
	int add_phase;
	int phase_on[gv.len_ss];
	
	/* initialiaze phase active in the considered assemblage */
	for (int i = 0; i < gv.len_ss; i++){
		phase_on[i] = 0;
	}

	/** 
		get initial conditions for active phases
	*/
	for (int i = 0; i < d->n_Ox; i++){
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
		/* pure endmembers as solution phase */
		if (d->ph_id_A[i][0] == 2){
			phase_on[ph_id] 		= 1;
			em_id 					= d->ph_id_A[i][3];

			for (int j = 0; j < SS_ref_db[ph_id].n_em; j++) {	
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
			
			for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
				cp[id_cp].p_em[ii]  = SS_ref_db[ph_id].p[ii];
			}
			for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
				cp[id_cp].dguess[ii]  = SS_ref_db[ph_id].iguess[ii];
				cp[id_cp].lvlxeos[ii] = SS_ref_db[ph_id].iguess[ii];
				cp[id_cp].xeos[ii]    = SS_ref_db[ph_id].iguess[ii];
			}
	
			gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
			gv.n_solvi[ph_id] 	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;
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
				cp[id_cp].lvlxeos[ii] = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
			}

			gv.id_solvi[ph_id][gv.n_solvi[ph_id]] = id_cp;
			gv.n_solvi[ph_id] 	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;
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
			max_n_pc = get_max_n_pc(SS_ref_db[i].tot_pc, SS_ref_db[i].n_pc);
			
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
					cp[id_cp].p_em[ii]      = SS_ref_db[ph_id].p_pc[pc_id][ii];
				}
				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					cp[id_cp].dguess[ii]    = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
					cp[id_cp].xeos[ii]      = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
					cp[id_cp].lvlxeos[ii]   = SS_ref_db[ph_id].xeos_pc[pc_id][ii];
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
  function to swp pure phases
*/	
void swap_pure_phases(				struct bulk_info 	 z_b,
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
				d->g0_A[d->ph2swp] 	   = d->g0_B;
				
				for (int j = 0; j < d->n_Ox; j++){				
					k = d->ph2swp + j*d->n_Ox;
					d->A[k] = d->B[j];
				}
				
				for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

				/** inverse guessed assemblage stoechiometry matrix */
				inverseMatrix(	d->A1,
								d->n_Ox		);
				
				/** update phase fractions */
				MatVecMul(		d->A1,
								// br,
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
void swap_pure_endmembers(				struct bulk_info 	 z_b,
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


					d->g0_B		 = SS_ref_db[i].gbase[l]*factor;
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
						d->g0_A[d->ph2swp] 	   = d->g0_B;
						
						for (int j = 0; j < d->n_Ox; j++){				
							int k = d->ph2swp + j*d->n_Ox;
							d->A[k] = d->B[j];
						}
						for (int k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

						/** inverse guessed assemblage stoechiometry matrix */
						inverseMatrix(	d->A1,
										d->n_Ox		);
						
						/** update phase fractions */
						MatVecMul(		d->A1,
										// br,
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
  function to swp pure phases
*/	
void swap_pseudocompounds(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k, max_n_pc;
	
	for (int i = 0; i < gv.len_ss; i++){										/**loop to pass informations from active endmembers */
		if (SS_ref_db[i].ss_flags[0] == 1){

			max_n_pc = get_max_n_pc(SS_ref_db[i].tot_pc, SS_ref_db[i].n_pc);

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
					SS_ref_db[i].n_swap[l]   = d->n_swp;
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
					inverseMatrix(	d->A1,
									d->n_Ox		);
					
					/** update phase fractions */
					MatVecMul(		d->A1,
									// br,
									z_b.bulk_rock_cat,
									d->n_vec,
									d->n_Ox		);
				}
			}
		}
	}

}

/**
  function to run simplex linear programming with pseudocompounds only
*/	
void run_simplex_vPC_only(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db
){
	simplex_data *d  = (simplex_data *) splx_data;

	int     k = 0;

	d->swp = 1;
	while (d->swp == 1){																			/** as long as a phase can be added to the guessed assemblage, go on */
		k 				 += 1;
		d->swp     = 0;
		
		swap_pure_phases(					z_b,
											splx_data,
											gv,
											PP_ref_db,
											SS_ref_db	);	
											
		swap_pure_endmembers(				z_b,
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
void run_simplex_vPC_stage(				struct bulk_info 	 z_b,
										simplex_data 		*splx_data,
										
										global_variable 	 gv,
										
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db,
										obj_type			*SS_objective
){
	simplex_data *d  = (simplex_data *) splx_data;

	int i, k, iss;

	/** copy A onto A1 in order to inverse it using LAPACKE */
	for (k = 0; k < d->n_Ox*d->n_Ox; k++){ d->A1[k] = d->A[k];}

	/** inverse guessed assemblage stoechiometry matrix */
	inverseMatrix(						d->A1, 
										d->n_Ox					);
	
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
	
	PC_ref 			SS_PC_xeos[gv.len_ss];
	
	for (iss = 0; iss < gv.len_ss; iss++){
		SS_PC_init_function(			SS_PC_xeos, 
										iss,
										gv.SS_list[iss]				);
	}
	
	for (iss = 0; iss < gv.len_ss; iss++){
		if (SS_ref_db[iss].ss_flags[0] == 1){

			generate_pseudocompounds(	iss,
										z_b,
										gv,
										SS_ref_db,
										SS_PC_xeos,
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
	run_simplex_vPC_only(				z_b,
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

/**
  function to deallocte memory of simplex linear programming (A)
*/	
void destroy_simplex_A(		simplex_data *splx_data
){

	simplex_data *d  = (simplex_data *) splx_data;

    for (int i = 0; i < d->n_Ox; i++){
		free(d->ph_id_A[i]);
	}
	free(d->A);
	free(d->pivot);
	free(d->A1);
	free(d->ph_id_A);
	free(d->g0_A);
	free(d->dG_A);
	free(d->n_vec);
	free(d->gamma_ps);
	free(d->gamma_ss);
	free(d->gamma_tot);
	free(d->gamma_delta);
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
global_variable run_levelling_function(		struct bulk_info 	 z_b,
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
	run_simplex_vPC_stage(					z_b,
										    splx_data,
											gv,
											PP_ref_db,
											SS_ref_db,
											SS_objective	);	
			
	/* update global variable gamma */
	update_global_gamma(					z_b,
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
		printf("══════════════════════════════════════════════\n");
		printf("   STEP 1: Pure species guess\n");
		printf("   ═══════════════════════════\n\n");
		printf("    P: %+10f T: %+10f\n",z_b.P,z_b.T);
		printf("\t[---------------------------------------]\n");
		printf("\t[  EM  |   EM PROP  |   g0_EM    |  ix  ]\n");
		printf("\t[---------------------------------------]\n");
		

		for (int i = 0; i < d->n_Ox; i++){
			if (d->ph_id_A[i][0] == 1){
				printf("\t['%5s' %+10f  %+10f  %5d ]", gv.PP_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				printf("\n");
			}
			if (d->ph_id_A[i][0] == 2){
				printf("\t['%5s' %+10f  %+10f  %5d ]\n", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
			}
			if (d->ph_id_A[i][0] == 3){
				printf("\t['%5s' %+10f  %+10f  %5d ]", gv.SS_list[d->ph_id_A[i][1]], d->n_vec[i], d->g0_A[i], d->ph_id_A[i][0]);
				for (int ii = 0; ii < SS_ref_db[d->ph_id_A[i][1]].n_xeos; ii++){
					printf(" %+10f", SS_ref_db[d->ph_id_A[i][1]].xeos_pc[d->ph_id_A[i][3]][ii] );
				}
				printf("\n");
			}
		}
		printf("\t[---------------------------------------]\n");
		printf("\t[  OXIDE      GAMMA_EM        GAMMA_PC  ]\n");
		printf("\t[---------------------------------------]\n");
		for (int i = 0; i < d->n_Ox; i++){
			printf("\t[ %5s %+15f %+15f ]\n", gv.ox[z_b.nzEl_array[i]], d->gamma_ps[i], d->gamma_tot[z_b.nzEl_array[i]]);
		}
		printf("\t[---------------------------------------]\n");
		printf("\t[            %4d swaps                 ]\n", d->n_swp);
		printf("\t[---------------------------------------]\n");
		
		printf("\n\t[---------------------------------------]\n");
		printf("\t[           ACTIVE PHASES               ]\n");
		printf("\t[---------------------------------------]\n");
		for (int i = 0; i < gv.len_ss; i++){
			if (SS_ref_db[i].ss_flags[0] == 1){
				printf("\t[                %5s                  ]\n",gv.SS_list[i]);
			}
		}
		printf("\t[---------------------------------------]\n");
		printf("\t[           UNACTIVE PHASES             ]\n");
		printf("\t[---------------------------------------]\n");
		for (int i = 0; i < gv.len_ss; i++){
			if (SS_ref_db[i].ss_flags[0] == 0){
				printf("\t[                %5s                  ]\n",gv.SS_list[i]);
			}
		}
		printf("\t[---------------------------------------]\n");
	}

	t 			= clock() - t; 
	time_taken  = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	gv.LVL_time = time_taken*1000;	
	
	return gv;
};		

/**
  main levelling routine
*/ 
global_variable Levelling(	struct bulk_info 	z_b,
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

	/** declare structure to store simplex arrays */
	// obj_type 		SS_objective[gv.len_ss];	
	
	// SS_objective_init_function(				SS_objective,
	// 										gv				);

	/* pseudosection function to get starting guess */
	gv = run_levelling_function(	z_b,												/** bulk rock informations    */
									gv,													/** global variables (e.g. Gamma) */

									SS_objective,
								    splx_data,
									PP_ref_db,											/** pure phase database */
									SS_ref_db,											/** solution phase database */
									cp				);
	if (gv.verbose == 1){
		printf("    [    Levelling time  %+12f ms    ]\n",gv.LVL_time);
		printf("    [---------------------------------------]\n\n\n");
	}


	return gv;
};		
