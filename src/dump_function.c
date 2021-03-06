/**
 The goal is to retrieve logged data and dump them into      
 files to track the behaviour of the solver                  
 This should be used in production runs as dumping files to  
 disk will drastically slow down the performances            
                                                             
 TRACKED DOWN DATA:                                          
  - Execution time per function                      
  - Residual                                                                                                     
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "mpi.h"
#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "toolkit.h"

/**
  Initialize dumping function by creating needed files
*/
void dump_init(global_variable gv){
	FILE *loc_min;
	char 	out_lm[255];
	struct 	stat st = {0};
	int 	rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (stat(gv.outpath, &st) == -1) {
    	mkdir(gv.outpath, 0700);
	}

	/** ----------------------------------------------------------------------------------------------- **/
	/** THERMOCALC LIKE FINAL OUTPUT **/
	if (gv.verbose == 1){
		sprintf(out_lm,	"%s_thermocalc_style_output.txt"		,gv.outpath); 
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "\n");
		fclose(loc_min);	
	}
	/** ----------------------------------------------------------------------------------------------- **/
	if (gv.verbose == 0){
		/** MATLAB GRID OUTPUT **/
		if (numprocs==1){	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath); 		}
		else 			{	sprintf(out_lm,	"%s_pseudosection_output.%i.txt"	,gv.outpath, rank); }
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "// NUMBER\tSTATUS[S,R1,R2,F]\tP[kbar]\tT[C]\tG_sys[G]\tBR_norm[wt]\tVp[km/s]\tVs[km/s]\tGAMMA[G] PHASE[name]\tMODE[wt]\tRHO[kg.m-3]\tX-EOS\n");
		fclose(loc_min);	


		if (gv.save_residual_evolution == 1){
			/** OUTPUT RESIDUAL EVOLUTION **/
			if (numprocs==1){	sprintf(out_lm,	"%s_residual_norm.txt"		,gv.outpath); 		}
			else 			{	sprintf(out_lm,	"%s_residual_norm.%i.txt"	,gv.outpath, rank); }
			loc_min 	= fopen(out_lm, 	"w"); 
			fclose(loc_min);	
		}
			
		/** MODE 2 - LOCAL MINIMA **/
		if (gv.Mode == 2){
			if (numprocs==1){	 sprintf(out_lm,	"%s__LOCAL_MINIMA.txt"		,gv.outpath); 	   }
			else 			{	 sprintf(out_lm,	"%s__LOCAL_MINIMA.%i.txt"	,gv.outpath, rank);}
			loc_min 	= fopen(out_lm, 	"w"); 
			fprintf(loc_min, "// PHASE_NAME[char]\tN_x-eos[n]\tN_POINTS\tGAMMA[G]\n");
			fprintf(loc_min, "// NUMBER\t INITIAL ENDMEMBER PROPORTIONS[n+1]\tINITIAL_GUESS_x_eos[n]\tFINAL_x-eos[n]\tFINAL ENDMEMBER PROPORTIONS[n+1]\tDRIVING_FORCE[dG]\n");
		
			fclose(loc_min);	
		}
		/** MODE 2 - LEVELLING_GAMMA **/
		if (gv.Mode == 3){
			if (numprocs==1){	 sprintf(out_lm,	"%s__LEVELLING_GAMMA.txt"		,gv.outpath); 	   }
			else 			{	 sprintf(out_lm,	"%s__LEVELLING_GAMMA.%i.txt"	,gv.outpath, rank);}
			loc_min 	= fopen(out_lm, 	"w"); 
			fprintf(loc_min, "// BULK-ROCK[len_ox]\tP[kbar]\tT[??C]\tGAMMA[G]\n");

			fclose(loc_min);	
		}
	}
}

/**
  Save final result of minimization
*/
void fill_output_struct(		global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	double G = 0.0;
	double sum;
	double sum_wt;
	double sum_em_wt;
	double sum_ph_mass;

	double sum_Molar_mass_bulk;
	double atp2wt;

	int nox  = gv.len_ox;
	int i, j, k, m, n;
	
	strcpy(sp[0].MAGEMin_ver,gv.version);	
	
	sp[0].bulk_res_norm 		 = gv.BR_norm;
	sp[0].n_iterations 		     = gv.global_ite;
	sp[0].status 		         = gv.div;
	
	sp[0].nOx					 = gv.len_ox;
	sp[0].rho					 = gv.system_density;
	sp[0].bulkMod				 = gv.system_bulkModulus;
	sp[0].shearMod				 = gv.system_shearModulus;
	sp[0].Vp					 = gv.system_Vp;
	sp[0].Vs					 = gv.system_Vs;

	sp[0].bulkModulus_M  		 = gv.melt_bulkModulus,
	sp[0].bulkModulus_S  		 = gv.solid_bulkModulus,
	sp[0].shearModulus_S  		 = gv.solid_shearModulus,
	sp[0].Vp_S  		 		 = gv.solid_Vp,
	sp[0].Vs_S  		 		 = gv.solid_Vs,

	sp[0].G = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		strcpy(sp[0].oxides[j],gv.ox[j]);	
		sp[0].G 				+= z_b.bulk_rock[j]*gv.gam_tot[j];
	}

	sp[0].P 			 		 = z_b.P;
	sp[0].T 			 		 = z_b.T;
	
	sum_Molar_mass_bulk = 0.0;
	for (i = 0; i < nox; i++){
		sp[0].bulk[i] 	 		 = z_b.bulk_rock[i];
		sp[0].gamma[i] 	 		 = gv.gam_tot[i];
		sp[0].bulk_wt[i] 	 	 = z_b.bulk_rock[i]*z_b.masspo[i];
		sum_Molar_mass_bulk     += sp[0].bulk_wt[i];
	}
	for (i = 0; i < nox; i++){
		sp[0].bulk_wt[i] 	 	/= sum_Molar_mass_bulk;
	}

	sp[0].n_ph			 		 = gv.n_phase;
	sp[0].n_PP			 		 = 0;
	sp[0].n_SS			 		 = 0;
	sp[0].frac_S				 = 0.0;
	sp[0].frac_M				 = 0.0;
	sp[0].frac_F				 = 0.0;
	sp[0].rho_S				 	 = 0.0;
	sp[0].rho_M				 	 = 0.0;
	sp[0].rho_F				 	 = 0.0;

	/* copy data from solution phases */
	n = 0;
	m = 0;
	for (int i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			strcpy(sp[0].ph[n],cp[i].name);	

			atp2wt = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				atp2wt	+= cp[i].ss_comp[j]*cp[i].factor*z_b.masspo[j];
			}
			atp2wt		/= sum_Molar_mass_bulk;

			sp[0].ph_frac[n]  	 = cp[i].ss_n;
			sp[0].ph_frac_wt[n]  = cp[i].ss_n*atp2wt;

			sp[0].ph_type[n]  	 = 1;
			sp[0].ph_id[n] 		 = m;
			sp[0].n_SS 			+= 1;
			
			G = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}
			
			sp[0].SS[m].nOx 	 = gv.len_ox;
			sp[0].SS[m].f 		 = cp[i].factor;
			sp[0].SS[m].G 		 = G;
			sp[0].SS[m].deltaG	 = cp[i].df;
			sp[0].SS[m].V 		 = cp[i].volume;
			sp[0].SS[m].cp 		 = cp[i].phase_cp;
			sp[0].SS[m].rho 	 = cp[i].phase_density;
			sp[0].SS[m].alpha 	 = cp[i].phase_expansivity;
			sp[0].SS[m].bulkMod  = cp[i].phase_bulkModulus/10.;
			sp[0].SS[m].shearMod = cp[i].phase_shearModulus/10.;
			sp[0].SS[m].Vp 		 = sqrt((cp[i].phase_bulkModulus/10. + 4.0/3.0*cp[i].phase_shearModulus/10.)/(cp[i].phase_density/1e3));
			sp[0].SS[m].Vs 		 = sqrt(cp[i].phase_shearModulus/10.0/(cp[i].phase_density/1e3));	

			sp[0].SS[m].n_xeos   = cp[i].n_xeos;
			sp[0].SS[m].n_em 	 = cp[i].n_em;

			/* solution phase composition */
			sum_wt = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				sp[0].SS[m].Comp[j]				= cp[i].ss_comp[j]*cp[i].factor;
				sp[0].SS[m].Comp_wt[j]			= sp[0].SS[m].Comp[j]*z_b.masspo[j];
				sum_wt 						   += sp[0].SS[m].Comp_wt[j];
			}
			for (j = 0; j < gv.len_ox; j++){
				sp[0].SS[m].Comp_wt[j]		   /= sum_wt;
			}
	
			for (j = 0; j < cp[i].n_xeos; j++){	
				sp[0].SS[m].compVariables[j] 	= cp[i].xeos[j];
			}
			sum_ph_mass = 0.0;
			for (j = 0; j < cp[i].n_em; j++){
				sum_em_wt = 0.0;
				for (k = 0; k < gv.len_ox; k++){
					sum_em_wt += SS_ref_db[cp[i].id].Comp[j][k]*cp[i].p_em[j]*z_b.masspo[k];
				}
				sp[0].SS[m].emFrac_wt[j] 		= sum_em_wt/sum_wt;
				sum_ph_mass					   += sp[0].SS[m].emFrac_wt[j];
				strcpy(sp[0].SS[m].emNames[j],SS_ref_db[cp[i].id].EM_list[j]);	
				sp[0].SS[m].emFrac[j] 			= cp[i].p_em[j];
				sp[0].SS[m].emChemPot[j] 		= cp[i].mu[j];
				
				sum = 0.0;
				for (k = 0; k < gv.len_ox; k++){
					sp[0].SS[m].emComp[j][k]	= SS_ref_db[cp[i].id].Comp[j][k]*cp[i].factor;
					sp[0].SS[m].emComp_wt[j][k]	= sp[0].SS[m].emComp[j][k]*z_b.masspo[k];
					sum 					   += sp[0].SS[m].emComp_wt[j][k];
				}
				for (k = 0; k < gv.len_ox; k++){
					sp[0].SS[m].emComp_wt[j][k]/= sum;
				}
			}
			for (j = 0; j < cp[i].n_em; j++){
				sp[0].SS[m].emFrac_wt[j] 	   /= sum_ph_mass;

			}

			if (strcmp( cp[i].name, "liq") == 0 || strcmp( cp[i].name, "fl") == 0 ){
				if (strcmp( cp[i].name, "liq") == 0){
					sp[0].frac_M 				= cp[i].ss_n;
					sp[0].rho_M  				= cp[i].phase_density;
					sum = 0.0;
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_M[j]	   		= cp[i].ss_comp[j]*cp[i].factor;
						sp[0].bulk_M_wt[j]	    = sp[0].bulk_M[j]*z_b.masspo[j];
						sum 				   += sp[0].bulk_M_wt[j];
					}
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_M_wt[j]	   /= sum;
					}
					atp2wt = sum/sum_Molar_mass_bulk;
					sp[0].frac_M_wt 		    = sp[0].frac_M*atp2wt;

				}
				else{
					sp[0].frac_F 				= cp[i].ss_n;
					sp[0].rho_F  				= cp[i].phase_density;
					sum = 0.0;
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_F[j]	   		= cp[i].ss_comp[j]*cp[i].factor;
						sp[0].bulk_F_wt[j]	   	= cp[i].ss_comp[j]*cp[i].factor*z_b.masspo[j];
						sum 				   += sp[0].bulk_F_wt[j];
					}
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_F_wt[j]	   /= sum;	
					}
					atp2wt = sum/sum_Molar_mass_bulk;
					sp[0].frac_F_wt 		    = sp[0].frac_F*atp2wt;
				}
			}
			else {
				sp[0].frac_S 				   += cp[i].ss_n;
				sp[0].rho_S  				   += cp[i].ss_n*cp[i].phase_density;
				for (j = 0; j < gv.len_ox; j++){
					sp[0].bulk_S[j]	   		   += cp[i].ss_n*cp[i].ss_comp[j]*cp[i].factor;
				}
			}

			n 					+= 1;
			m 					+= 1;
		}
	}
	/* copy data from pure phases */
	m = 0;
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			strcpy(sp[0].ph[n],gv.PP_list[i]);

			atp2wt = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				atp2wt	+= PP_ref_db[i].Comp[j]*PP_ref_db[i].factor*z_b.masspo[j];
			}
			atp2wt		/= sum_Molar_mass_bulk;

			sp[0].ph_frac[n]  	 = gv.pp_n[i];
			sp[0].ph_frac_wt[n]  = gv.pp_n[i]*atp2wt;
			sp[0].ph_type[n]  	 = 0;
			sp[0].ph_id[n] 		 = m;
			sp[0].n_PP 			+= 1;
					
			sp[0].PP[m].nOx 	 = gv.len_ox;
			sp[0].PP[m].f 		 = PP_ref_db[i].factor;
			sp[0].PP[m].G 		 = PP_ref_db[i].gbase;
			sp[0].PP[m].deltaG	 = PP_ref_db[i].gb_lvl;
			sp[0].PP[m].V 		 = PP_ref_db[i].volume;
			sp[0].PP[m].cp 		 = PP_ref_db[i].phase_cp;
			sp[0].PP[m].rho 	 = PP_ref_db[i].phase_density;
			sp[0].PP[m].alpha 	 = PP_ref_db[i].phase_expansivity;
			sp[0].PP[m].bulkMod  = PP_ref_db[i].phase_bulkModulus/10.;
			sp[0].PP[m].shearMod = PP_ref_db[i].phase_shearModulus/10.;
			sp[0].PP[m].Vp 		 = sqrt((PP_ref_db[i].phase_bulkModulus/10. + 4.0/3.0*PP_ref_db[i].phase_shearModulus/10.)/(PP_ref_db[i].phase_density/1e3));
			sp[0].PP[m].Vs 		 = sqrt(PP_ref_db[i].phase_shearModulus/10.0/(PP_ref_db[i].phase_density/1e3));	

			sum = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				sp[0].PP[m].Comp[j]		 = PP_ref_db[i].Comp[j]*PP_ref_db[i].factor;
				sp[0].PP[m].Comp_wt[j]   = sp[0].PP[m].Comp[j]*z_b.masspo[j];
				sum 					+= sp[0].PP[m].Comp_wt[j];
			}
			for (j = 0; j < gv.len_ox; j++){
				sp[0].PP[m].Comp_wt[j]  /= sum;
			}
		
			sp[0].frac_S 		+= gv.pp_n[i];
			sp[0].rho_S  		+= gv.pp_n[i]*PP_ref_db[i].phase_density;
			for (j = 0; j < gv.len_ox; j++){
				sp[0].bulk_S[j]	+= gv.pp_n[i]*PP_ref_db[i].Comp[j]*PP_ref_db[i].factor;;
			}
			n 			    	+= 1;
			m 					+= 1;
		}
	}
	
	/* normalize rho_S and bulk_S */
	sp[0].rho_S  				/= sp[0].frac_S;
	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_S[j]	   		/= sp[0].frac_S;
	}

	sum = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_S_wt[j]	    = sp[0].bulk_S[j]*z_b.masspo[j];
		sum 				   += sp[0].bulk_S_wt[j];
	}

	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_S_wt[j]	   /= sum;
	}

	atp2wt = sum/sum_Molar_mass_bulk;
	sp[0].frac_S_wt 		    = sp[0].frac_S*atp2wt;



	// debug print
	if (1 == 0){
		for (int m = 0; m < gv.n_phase; m++){
			printf(" %4s %+10f\n",sp[0].ph[m],sp[0].ph_frac_wt[m]);
		}
		printf("\n");
		for (int m = 0; m < gv.n_cp_phase; m++){
			printf(" %4s composition [wt]\n",sp[0].ph[m]);
			for (j = 0; j < gv.len_ox; j++){
				printf(" %+10f", sp[0].SS[m].Comp_wt[j]);
			}
			printf("\n");

			for (int k = 0; k < sp[0].SS[m].n_em; k++){
				printf(" %+10f",sp[0].SS[m].emFrac_wt[k]);
			}
			printf("\n");

		}
		n = 0;
		for (int m = gv.n_cp_phase; m < gv.n_phase; m++){
			printf(" %4s composition [wt]\n",sp[0].ph[m]);
			for (j = 0; j < gv.len_ox; j++){
				printf(" %+10f", sp[0].PP[n].Comp_wt[j]);
			}
			n += 1;
			printf("\n");
		}

		printf("Bulk solid:\n  %+10f |",sp[0].frac_S_wt );
		for (j = 0; j < gv.len_ox; j++){
			printf(" %+10f", sp[0].bulk_S_wt[j]);
		}
		printf("\n");

		printf("Bulk melt:\n  %+10f |",sp[0].frac_M_wt );
		for (j = 0; j < gv.len_ox; j++){
			printf(" %+10f", sp[0].bulk_M_wt[j]);
		}
		printf("\n");

		printf("Bulk fluid:\n  %+10f |",sp[0].frac_F_wt );
		for (j = 0; j < gv.len_ox; j++){
			printf(" %+10f", sp[0].bulk_F_wt[j]);
		}
		printf("\n");
	}
	
}

/**
  Save final result of minimization
*/
void dump_results_function(		global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp
){
	
	FILE *loc_min;
	char out_lm[255];
	
	int i,j, rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/** ----------------------------------------------------------------------------------------------- **/
	/** THERMOCALC LIKE FINAL OUTPUT **/
	if (gv.verbose == 1){
		/* output active phase fraction*/
		if (numprocs==1){	sprintf(out_lm,	"%s_thermocalc_style_output.txt"		,gv.outpath); 	}
		else 			{	sprintf(out_lm,	"%s_thermocalc_style_output.%i.txt"		,gv.outpath, rank); 	}

		loc_min 	= fopen(out_lm, 	"a"); 
		fprintf(loc_min, "============================================================\n");
		
		for (int i = 0; i < gv.len_cp; i++){
			if ( cp[i].ss_flags[1] == 1){
				fprintf(loc_min, 	"%4s ", cp[i].name);

			}
		}
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){
				fprintf(loc_min, 	"%4s ", gv.PP_list[i]);
			}
		}	
		fprintf(loc_min, " {%10.5f, %10.5f} kbar/??C\n\n",z_b.P,z_b.T-273.15);
		
		fprintf(loc_min, "Compositional variables (solution phase):\n");		
		for (i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				fprintf(loc_min, 	" %5s", cp[i].name);
				for (j = 0; j < cp[i].n_xeos; j++){
					fprintf(loc_min, 	"%10.5f ", cp[i].xeos[j]);
				}
				for (int k = j; k < gv.len_ox; k++){
					fprintf(loc_min, 	"%10s ", "-");
				}		
				fprintf(loc_min, "\n");	
			}
		}
		
		fprintf(loc_min, "\nEnd-members fraction (solution phase):\n");		
		for (i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				
				fprintf(loc_min, 	" %5s", "");
				for (j = 0; j < cp[i].n_em; j++){
					fprintf(loc_min, 	"%10s ", SS_ref_db[cp[i].id].EM_list[j]);
				}
				for (int k = j; k < gv.len_ox; k++){
					fprintf(loc_min, 	"%10s ", "-");
				}		
				fprintf(loc_min, "\n");	
				
				fprintf(loc_min, 	" %5s", cp[i].name);
				for (j = 0; j < cp[i].n_em; j++){
					fprintf(loc_min, 	"%10.5f ", cp[i].p_em[j]);
				}
				for (int k = j; k < gv.len_ox; k++){
					fprintf(loc_min, 	"%10s ", "-");
				}		
				fprintf(loc_min, "\n");	
			}
		}
		
		if (gv.Mode == 1){
			fprintf(loc_min, "\nGibbs energy of reference G0 (solution phase):\n");		
			for (i = 0; i < gv.len_cp; i++){
				if (cp[i].ss_flags[1] == 1){
					fprintf(loc_min, 	" %5s", cp[i].name);
					for (j = 0; j < cp[i].n_em; j++){
						fprintf(loc_min, 	"%14.6f ", cp[i].gbase[j]);
					}
					for (int k = j; k < gv.len_ox; k++){
						fprintf(loc_min, 	"%14s ", "-");
					}		
					fprintf(loc_min, "\n");	
				}
			}
					
			fprintf(loc_min, "\nChemical potentials [J] (solution phase):\n");		
			for (i = 0; i < gv.len_cp; i++){
				if (cp[i].ss_flags[1] == 1){
					fprintf(loc_min, 	" %5s", cp[i].name);
					for (j = 0; j < cp[i].n_em; j++){
						fprintf(loc_min, 	"%14.6f ", cp[i].mu[j]);
					}
					for (int k = j; k < gv.len_ox; k++){
						fprintf(loc_min, 	"%14s ", "-");
					}		
					fprintf(loc_min, "\n");	
				}
			}
		}

		fprintf(loc_min, "\nSite fractions:\n");		
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				fprintf(loc_min, 	" %5s", cp[i].name);
				for (j = 0; j < (cp[i].n_sf); j++){
					fprintf(loc_min, 	"%10.5f ", cp[i].sf[j]); // *-1.0 because inequality are given as -x <= 0 in NLopt
				}
				for (int k = j; k < gv.len_ox; k++){
					fprintf(loc_min, 	"%10s ", "-");
				}		
				fprintf(loc_min, "\n");	
			}
		}
		
		fprintf(loc_min, "\nOxide compositions [mol%%] (normalized on 1 atom basis):\n");	
		fprintf(loc_min, "%5s"," ");
		for (i = 0; i < gv.len_ox; i++){
			fprintf(loc_min, " %10s", gv.ox[i]);
		}
		fprintf(loc_min, "\n %5s","SYS");	
		for (int i = 0; i < gv.len_ox; i++){
			fprintf(loc_min, "%10.5f ",z_b.bulk_rock[i]);
		}
		fprintf(loc_min, "\n");	
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				fprintf(loc_min, 	" %5s", cp[i].name);
				for (int j = 0; j < gv.len_ox; j++){
					fprintf(loc_min, "%10.5f ", cp[i].ss_comp[j]*cp[i].factor);
				}
				fprintf(loc_min, "\n");
			}
		}

		for (int i = 0; i < gv.len_pp; i++){ 
			if (gv.pp_flags[i][1] == 1){
				fprintf(loc_min, 	" %5s", gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					fprintf(loc_min, "%10.5f ", PP_ref_db[i].Comp[j]*PP_ref_db[i].factor);
				}
				fprintf(loc_min, "\n");
			}
		}
		
		double G;
		fprintf(loc_min, "\n");	
		fprintf(loc_min, "Stable mineral assemblage:\n");	
		fprintf(loc_min, "%6s%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","phase","mode","f","G" ,"V" ,"Cp","rho","Thermal_Exp","BulkMod[GPa]","ShearMod[GPa]","Vp[km/s]","Vs[km/s]");
					
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				
				if (gv.Mode == 1){
					G = cp[i].df;
				}
				else{
					G = 0.0;
					for (int j = 0; j < gv.len_ox; j++){
						G += cp[i].ss_comp[j]*gv.gam_tot[j];
					}
				}

				fprintf(loc_min, "%6s", cp[i].name);
				fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.2f %+12.2f %+12.2f %+12.2f",cp[i].ss_n,cp[i].factor,G,cp[i].volume,cp[i].phase_cp,cp[i].phase_density,cp[i].phase_expansivity,cp[i].phase_bulkModulus/10.,cp[i].phase_shearModulus/10.,sqrt((cp[i].phase_bulkModulus/10. +4.0/3.0*cp[i].phase_shearModulus/10.)/(cp[i].phase_density/1e3)),sqrt(cp[i].phase_shearModulus/10.0/(cp[i].phase_density/1e3)) );
				fprintf(loc_min, "\n");
			}
		}
		
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){ 
				fprintf(loc_min, "%6s", gv.PP_list[i]);
				fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.2f %+12.2f %+12.2f %+12.2f",gv.pp_n[i],PP_ref_db[i].factor,PP_ref_db[i].gbase,PP_ref_db[i].volume,PP_ref_db[i].phase_cp,PP_ref_db[i].phase_density,PP_ref_db[i].phase_expansivity,PP_ref_db[i].phase_bulkModulus/10.,PP_ref_db[i].phase_shearModulus/10.,sqrt((PP_ref_db[i].phase_bulkModulus/10. +4.0/3.0*PP_ref_db[i].phase_shearModulus/10.)/(PP_ref_db[i].phase_density/1e3)),sqrt(PP_ref_db[i].phase_shearModulus/10.0/(PP_ref_db[i].phase_density/1e3)));
				fprintf(loc_min, "\n");
			}
		}
		
		G = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			G += z_b.bulk_rock[j]*gv.gam_tot[j];
		}
		fprintf(loc_min, "%6s %24s %+12.5f %25s %+12.5f %12s %+12.5f %+12.5f %+12.5f %+12.5f\n","SYS"," ",G," ",gv.system_density," ",gv.system_bulkModulus,gv.system_shearModulus,gv.system_Vp,gv.system_Vs);
		
		/* output solution phase composition */ 	
		fprintf(loc_min, "\nGamma (chemical potential of oxides):\n");
		for (i = 0; i < gv.len_ox; i++){
			fprintf(loc_min, "%6s %+12.5f\n", gv.ox[i], gv.gam_tot[i]);
		}
		fprintf(loc_min, "\ndelta Gibbs energy (G-hyperplane distance):\n");
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				fprintf(loc_min, 	"%5s %+10e\n", cp[i].name,cp[i].df);
			}
		}
		fprintf(loc_min, "\n\n");
		fclose(loc_min);	
	}
	
	/** ----------------------------------------------------------------------------------------------- **/
	/** MATLAB GRID OUTPUT **/
	if (gv.verbose == 0){
		if (numprocs==1){	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath); 		}
		else 			{	sprintf(out_lm,	"%s_pseudosection_output.%i.txt"	,gv.outpath, rank); }
		/* get number of repeated phases for the solvi */
		int n_solvi[gv.len_ss];
		for (int i = 0; i < gv.len_ss; i++){
			n_solvi[i] = 0;
		}
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[1] == 1){
				n_solvi[cp[i].id] += 1;
			}
		}
	
		loc_min 	= fopen(out_lm, 	"a"); 
		fprintf(loc_min, "%i %i %.10f %.10f %.10f %.10f", gv.numPoint+1, gv.status, z_b.P, z_b.T-273.15, gv.G_system,gv.BR_norm);

		for (i = 0; i < gv.len_ox; i++){
			fprintf(loc_min," %0.10f", gv.gam_tot[i]);
		}

		// fprintf(loc_min, " %.10f %.10f",gv.system_Vp,gv.system_Vs);
		fprintf(loc_min, " %.10f %.10f",gv.V_cor[0],gv.V_cor[1]);
		fprintf(loc_min, "\n");
		for (int i = 0; i < gv.len_cp; i++){ 
			if (cp[i].ss_flags[1] == 1){
				
				
				if (n_solvi[cp[i].id] > 1){
					fprintf(loc_min, 	"%s_%d \t %.10f \t %.10f \t", cp[i].name,n_solvi[cp[i].id], cp[i].ss_n, cp[i].phase_density);
				}
				else{
					fprintf(loc_min, 	"%s \t %.10f \t %.10f \t", cp[i].name, cp[i].ss_n, cp[i].phase_density);
				}

				fprintf(loc_min, 	"%d ", cp[i].n_xeos);
				for (int j = 0; j < (cp[i].n_xeos); j++){
					fprintf(loc_min, 	"%.10f ", cp[i].xeos[j]);
				}
				for (int j = 0; j < (cp[i].n_em); j++){
					fprintf(loc_min, 	"%10s ",  SS_ref_db[cp[i].id].EM_list[j]);	
					fprintf(loc_min, 	"%.10f ", cp[i].p_em[j]);
				}
				fprintf(loc_min, "\n");
			}
		}	
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){
				fprintf(loc_min, 	"%s \t %.10f \t %.10f \t", gv.PP_list[i], gv.pp_n[i], PP_ref_db[i].phase_density);
				fprintf(loc_min, "\n");
			}
		}	
		fprintf(loc_min, "\n");
		fclose(loc_min);

		if (gv.save_residual_evolution == 1){
			/** OUTPUT RESIDUAL EVOLUTION **/
			if (numprocs==1){	sprintf(out_lm,	"%s_residual_norm.txt"		,gv.outpath); 		}
			else 			{	sprintf(out_lm,	"%s_residual_norm.%i.txt"	,gv.outpath, rank); }

			loc_min 	= fopen(out_lm, 	"a"); 

			for (int j = 0; j < gv.global_ite; j++){
				fprintf(loc_min, "%.6f ", gv.PGE_mass_norm[j]);
			}
			fprintf(loc_min, "\n");

			fclose(loc_min);
		}

	}
};

/**
  Parallel file dump for phase diagrams
*/
void mergeParallelFiles(global_variable gv){

	int i, rank, numprocs,  MAX_LINE_LENGTH=200;
	char out_lm[255];
	char in_lm[255];
	char c; 
	char buf[MAX_LINE_LENGTH];
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (numprocs == 1){ return; }

	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath);
   	FILE *fp2 = fopen(out_lm, "w"); 

	fprintf(fp2, "// NUMBER\tSTATUS[S,R1,R2,F]\tP[kbar]\tT[C]\tG_sys[G]\tBR_norm[wt]\tVp[km/s]\tVs[km/s]\tGAMMA[G]; PHASE[name]\tMODE[wt]\tRHO[kg.m-3]\tX-EOS\n");

	// Open file to be merged 
	for (i = 0; i < numprocs; i++){
		// open file
		sprintf(in_lm,	"%s_pseudosection_output.%i.txt"		,gv.outpath, i);
		FILE *fp1 = fopen(in_lm, "r"); 
		
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip first line = comment (we don't want to copy that)
	
		// Copy contents of first file to file3.txt 
		while ((c = fgetc(fp1)) != EOF){ 
			fputc(c, fp2); 
		}
		fclose(fp1); 
	}
   fclose(fp2);
}

/**
  Parallel file dump for phase diagrams
*/
void mergeParallel_residual_Files(global_variable gv){

	int i, rank, numprocs,  MAX_LINE_LENGTH=2048;
	char out_lm[255];
	char in_lm[255];
	char c; 
	char buf[MAX_LINE_LENGTH];
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (numprocs == 1){ return; }

	sprintf(out_lm,	"%s_residual_norm.txt"		,gv.outpath);
   	FILE *fp2 = fopen(out_lm, "w"); 

	// Open file to be merged 
	for (i = 0; i < numprocs; i++){
		// open file
		sprintf(in_lm,	"%s_residual_norm.%i.txt"		,gv.outpath, i);
		FILE *fp1 = fopen(in_lm, "r"); 
		
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip first line = comment (we don't want to copy that)
	
		// Copy contents of first file to file3.txt 
		while ((c = fgetc(fp1)) != EOF){ 
			fputc(c, fp2); 
		}
		fclose(fp1); 
	}
   fclose(fp2);
}


/**
  Parallel file dump for local minima search
*/
void mergeParallel_LocalMinima_Files(global_variable gv){

	int i, rank, numprocs,  MAX_LINE_LENGTH=200;
	char out_lm[255];
	char in_lm[255];
	char c; 
	char buf[MAX_LINE_LENGTH];
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (numprocs == 1){ return; }

	sprintf(out_lm,	"%s__LOCAL_MINIMA.txt"		,gv.outpath);
   	FILE *fp2 = fopen(out_lm, "w"); 

	fprintf(fp2, "// PHASE_NAME[char]\tN_x-eos[n]\tN_POINTS\tGAMMA[G]\n");
	fprintf(fp2, "// NUMBER\t INITIAL ENDMEMBER PROPORTIONS[n+1]\tINITIAL_GUESS_x_eos[n]\tFINAL_x-eos[n]\tFINAL ENDMEMBER PROPORTIONS[n+1]\tDRIVING_FORCE[dG]\n");
	

	// Open file to be merged 
	for (i = 0; i < numprocs; i++){
		// open file
		sprintf(in_lm,	"%s__LOCAL_MINIMA.%i.txt"		,gv.outpath, i);
		FILE *fp1 = fopen(in_lm, "r"); 
		
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip 1th line = comment (we don't want to copy that)
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip 2nd line = comment (we don't want to copy that)
		if (i>0){
			fgets(buf, MAX_LINE_LENGTH, fp1);				// skip 3rd line = info about Gamma (only needed once)
		}
	
		// Copy contents of first file to file3.txt 
		while ((c = fgetc(fp1)) != EOF){ 
			fputc(c, fp2); 
		}
		fclose(fp1); 
	}
   fclose(fp2); 
}

/**
  Parallel file dump for first stage of levelling minimization only
*/
void mergeParallel_LevellingGamma_Files(global_variable gv){

	int i, rank, numprocs,  MAX_LINE_LENGTH=200;
	char out_lm[255];
	char in_lm[255];
	char c; 
	char buf[MAX_LINE_LENGTH];
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (numprocs == 1){ return; }

	sprintf(out_lm,	"%s__LEVELLING_GAMMA.txt"		,gv.outpath);
   	FILE *fp2 = fopen(out_lm, "w"); 

	fprintf(fp2, "// BULK-ROCK[len_ox]\tP[kbar]\tT[??C]\tGAMMA[G]\n");

	// Open file to be merged 
	for (i = 0; i < numprocs; i++){
		// open file
		sprintf(in_lm,	"%s__LEVELLING_GAMMA.%i.txt"		,gv.outpath, i);
		FILE *fp1 = fopen(in_lm, "r"); 
		
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip 1th line = comment (we don't want to copy that)
		fgets(buf, MAX_LINE_LENGTH, fp1);					// skip 2nd line = comment (we don't want to copy that)
		if (i>0){
			fgets(buf, MAX_LINE_LENGTH, fp1);					// skip 3rd line = info about Gamma (only needed once)
		}
	
		// Copy contents of first file to file3.txt 
		while ((c = fgetc(fp1)) != EOF){ 
			fputc(c, fp2); 
		}
		fclose(fp1); 
	}
   fclose(fp2); 
}

