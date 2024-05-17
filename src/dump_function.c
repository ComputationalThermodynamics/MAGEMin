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
#include "objective_functions.h"
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


	/** THERMOCALC-LIKE OUTPUT **/
	if (gv.verbose == 1 && gv.output_matlab == 0){
		sprintf(out_lm,	"%s_thermocalc_style_output.txt"		,gv.outpath); 
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "\n");
		fclose(loc_min);	
	}
	/** OUTPUT FOR MATLAB, print out wt and mol fractions instead of TC atom basis **/
	if (gv.output_matlab >= 1){
		if (numprocs==1){	sprintf(out_lm,	"%s_matlab_output.txt"		,gv.outpath); 		}
		else 			{	sprintf(out_lm,	"%s_matlab_output.%i.txt"	,gv.outpath, rank); }
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "\n");
		fclose(loc_min);
	}

	if (gv.verbose == 0){
		/**GUI OUTPUT **/
		if (numprocs==1){	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath); 		}
		else 			{	sprintf(out_lm,	"%s_pseudosection_output.%i.txt"	,gv.outpath, rank); }
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "// {number status[] P[kbar] T[C] G_sys[G] BR_norm[wt] Gamma[G] Vp[km/s] Vs[km/s] entropy[J/K]} nextline {Phase[name] mode[wt] density[kg.m-3] x-eos}\n");
		fclose(loc_min);	
	}
}

/**
  Save final result of minimization
*/
void fill_output_struct(		global_variable 	 gv,
								simplex_data	    *splx_data,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	double G = 0.0;
	double sum;
	double sum_wt;
	double sum_mol;
	double sum_vol;
	double sum_em_wt;
	double sum_ph_mass;

	double sum_Molar_mass_bulk;
	double atp2wt;

	int nox  = gv.len_ox;
	int i, j, k, m, n, em_id, ph_id, pc_id;
	
	strcpy(sp[0].MAGEMin_ver,gv.version);	


	if (gv.EM_database == 0){	
		strcpy(sp[0].dataset,"tc_ds62");	
	}
	else if (gv.EM_database == 1){	
		strcpy(sp[0].dataset,"tc_ds62");	
	}
	else if (gv.EM_database == 2){		
		strcpy(sp[0].dataset,"tc_ds634");	
	}
	else if (gv.EM_database == 4){	
		strcpy(sp[0].dataset,"tc_ds633");	
	}

	sp[0].bulk_res_norm 		 = gv.BR_norm;
	sp[0].n_iterations 		     = gv.global_ite;
	sp[0].status 		         = gv.status;
	
	sp[0].nOx					 = gv.len_ox;
	sp[0].rho					 = gv.system_density;
	sp[0].fO2					 = log10(gv.system_fO2);
	sp[0].dQFM					 = log10(gv.system_deltaQFM);
	sp[0].aH2O					 = gv.system_aH2O;
	sp[0].aSiO2					 = gv.system_aSiO2;
	sp[0].aTiO2					 = gv.system_aTiO2;
	sp[0].aAl2O3				 = gv.system_aAl2O3;
	sp[0].aMgO					 = gv.system_aMgO;
	sp[0].aFeO				 	 = gv.system_aFeO;

	sp[0].alpha				 	 = gv.system_expansivity;
	sp[0].V				 	 	 = gv.system_volume*10.0;	
	sp[0].cp				 	 = gv.system_cp;	
	sp[0].entropy				 = gv.system_entropy;
	sp[0].enthalpy				 = gv.system_enthalpy;

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
	sp[0].X 			 		 = 1.0;
	
	sum_Molar_mass_bulk = 0.0;
	sp[0].M_sys  		= 0.0;
	for (i = 0; i < nox; i++){
		sp[0].M_sys 			+= z_b.bulk_rock[i]*z_b.masspo[i];	
		sp[0].bulk[i] 	 		 = z_b.bulk_rock[i];
		sp[0].gamma[i] 	 		 = gv.gam_tot[i];
		sp[0].bulk_wt[i] 	 	 = z_b.bulk_rock[i]*z_b.masspo[i];
		sum_Molar_mass_bulk     += sp[0].bulk_wt[i];

		sp[0].bulk_S[i] 	 		 = 0.0;
		sp[0].bulk_S_wt[i] 	 		 = 0.0;
		sp[0].bulk_M[i] 	 		 = 0.0;
		sp[0].bulk_M_wt[i] 	 		 = 0.0;
		sp[0].bulk_F[i] 	 		 = 0.0;
		sp[0].bulk_F_wt[i] 	 		 = 0.0;
	}
	for (i = 0; i < nox; i++){
		sp[0].bulk_wt[i] 	 	/= sum_Molar_mass_bulk;
	}

	sp[0].n_ph			 		 = gv.n_phase;
	sp[0].n_PP			 		 = 0;
	sp[0].n_SS			 		 = 0;
	sp[0].n_mSS			 		 = 0;
	sp[0].frac_S				 = 0.0;
	sp[0].frac_M				 = 0.0;
	sp[0].frac_F				 = 0.0;
	sp[0].frac_S_wt				 = 0.0;
	sp[0].frac_M_wt				 = 0.0;
	sp[0].frac_F_wt				 = 0.0;
	sp[0].rho_S				 	 = 0.0;
	sp[0].rho_M				 	 = 0.0;
	sp[0].rho_F				 	 = 0.0;
	sp[0].cp_wt					 = 0.0;
	sp[0].s_cp					 = 0.0;

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
			sp[0].cp_wt 			+= cp[i].phase_cp * cp[i].ss_n*atp2wt * cp[i].factor;
			G = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}
			
			sp[0].SS[m].nOx 	 = gv.len_ox;
			sp[0].SS[m].f 		 = cp[i].factor;
			sp[0].SS[m].G 		 = G;
			sp[0].SS[m].deltaG	 = cp[i].df;
			sp[0].SS[m].V 		 = cp[i].volume*10.;
			sp[0].SS[m].cp 		 = cp[i].phase_cp;
			sp[0].SS[m].rho 	 = cp[i].phase_density;
			sp[0].SS[m].alpha 	 = cp[i].phase_expansivity;
			sp[0].SS[m].entropy  = cp[i].phase_entropy;
			sp[0].SS[m].enthalpy = cp[i].phase_enthalpy;
			sp[0].SS[m].bulkMod  = cp[i].phase_bulkModulus/10.;
			sp[0].SS[m].shearMod = cp[i].phase_shearModulus/10.;
			sp[0].SS[m].Vp 		 = sqrt((cp[i].phase_bulkModulus/10. + 4.0/3.0*cp[i].phase_shearModulus/10.)/(cp[i].phase_density/1e3));
			sp[0].SS[m].Vs 		 = sqrt(cp[i].phase_shearModulus/10.0/(cp[i].phase_density/1e3));	

			sp[0].SS[m].n_xeos   = cp[i].n_xeos;
			sp[0].SS[m].n_em 	 = cp[i].n_em;

			/* solution phase composition */
			sum_wt = 0.0;
			sum_mol = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				sp[0].SS[m].Comp[j]				= cp[i].ss_comp[j]*cp[i].factor;
				sp[0].SS[m].Comp_wt[j]			= sp[0].SS[m].Comp[j]*z_b.masspo[j];
				sum_wt 						   += sp[0].SS[m].Comp_wt[j];
				sum_mol 					   += sp[0].SS[m].Comp[j];
			}
			for (j = 0; j < gv.len_ox; j++){
				sp[0].SS[m].Comp_wt[j]		   /= sum_wt;
				sp[0].SS[m].Comp[j]		   	   /= sum_mol;
			}
	
			for (j = 0; j < cp[i].n_xeos; j++){	
				sp[0].SS[m].compVariables[j] 	= cp[i].xeos[j];
			}

			for (j = 0; j < cp[i].n_xeos; j++){	
				strcpy(sp[0].SS[m].compVariablesNames[j],SS_ref_db[cp[i].id].CV_list[j]);	
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
				
				sum_wt  = 0.0;
				sum_mol = 0.0;
				for (k = 0; k < gv.len_ox; k++){
					sp[0].SS[m].emComp[j][k]	= SS_ref_db[cp[i].id].Comp[j][k]*cp[i].factor;
					sp[0].SS[m].emComp_wt[j][k]	= sp[0].SS[m].emComp[j][k]*z_b.masspo[k];
					sum_wt 					   += sp[0].SS[m].emComp_wt[j][k];
					sum_mol 				   += sp[0].SS[m].emComp[j][k];
				}
				for (k = 0; k < gv.len_ox; k++){
					sp[0].SS[m].emComp_wt[j][k]	/= sum_wt;
					sp[0].SS[m].emComp[j][k]	/= sum_mol;
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
					sum_mol = 0.0;
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_F[j]	   		= cp[i].ss_comp[j]*cp[i].factor;
						sp[0].bulk_F_wt[j]	   	= cp[i].ss_comp[j]*cp[i].factor*z_b.masspo[j];
						sum 				   += sp[0].bulk_F_wt[j];
						sum_mol		   		   += sp[0].bulk_F[j];
					}
					for (j = 0; j < gv.len_ox; j++){
						sp[0].bulk_F_wt[j]	   /= sum;
						sp[0].bulk_F[j]	   	   /= sum_mol;	
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
			sp[0].cp_wt 		+= PP_ref_db[i].phase_cp * gv.pp_n[i]*atp2wt *PP_ref_db[i].factor;

			sp[0].PP[m].nOx 	 = gv.len_ox;
			sp[0].PP[m].f 		 = PP_ref_db[i].factor;
			sp[0].PP[m].G 		 = PP_ref_db[i].gbase;
			sp[0].PP[m].deltaG	 = PP_ref_db[i].gb_lvl;
			sp[0].PP[m].V 		 = PP_ref_db[i].volume*10.;
			sp[0].PP[m].cp 		 = PP_ref_db[i].phase_cp;
			sp[0].PP[m].rho 	 = PP_ref_db[i].phase_density;
			sp[0].PP[m].alpha 	 = PP_ref_db[i].phase_expansivity;
			sp[0].PP[m].entropy  = PP_ref_db[i].phase_entropy;
			sp[0].PP[m].enthalpy = PP_ref_db[i].phase_enthalpy;
			sp[0].PP[m].bulkMod  = PP_ref_db[i].phase_bulkModulus/10.;
			sp[0].PP[m].shearMod = PP_ref_db[i].phase_shearModulus/10.;
			sp[0].PP[m].Vp 		 = sqrt((PP_ref_db[i].phase_bulkModulus/10. + 4.0/3.0*PP_ref_db[i].phase_shearModulus/10.)/(PP_ref_db[i].phase_density/1e3));
			sp[0].PP[m].Vs 		 = sqrt(PP_ref_db[i].phase_shearModulus/10.0/(PP_ref_db[i].phase_density/1e3));	

			sum_wt = 0.0;
			sum_mol = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				sp[0].PP[m].Comp[j]		 = PP_ref_db[i].Comp[j]*PP_ref_db[i].factor;
				sp[0].PP[m].Comp_wt[j]   = sp[0].PP[m].Comp[j]*z_b.masspo[j];
				sum_wt 					+= sp[0].PP[m].Comp_wt[j];
				sum_mol					+= sp[0].PP[m].Comp[j];
			}
			for (j = 0; j < gv.len_ox; j++){
				sp[0].PP[m].Comp_wt[j]  /= sum_wt;
				sp[0].PP[m].Comp[j]  	/= sum_mol;
			}
		
			if  (strcmp( gv.PP_list[i], "H2O") != 0){
				sp[0].frac_S 		+= gv.pp_n[i];
				sp[0].rho_S  		+= gv.pp_n[i]*PP_ref_db[i].phase_density;
				for (j = 0; j < gv.len_ox; j++){
					sp[0].bulk_S[j]	+= gv.pp_n[i]*PP_ref_db[i].Comp[j]*PP_ref_db[i].factor;;
				}
			}
			if  (strcmp( gv.PP_list[i], "H2O") == 0){
				sp[0].frac_F 		= gv.pp_n[i];
				sp[0].rho_F  		= gv.pp_n[i]*PP_ref_db[i].phase_density;
				for (j = 0; j < gv.len_ox; j++){
					sp[0].bulk_F[j]	= gv.pp_n[i]*PP_ref_db[i].Comp[j]*PP_ref_db[i].factor;;
				}
			}
			n 			    	+= 1;
			m 					+= 1;
		}
	}
	

	// compute volume fraction and normalize other fractions
	sum_vol = 0.0;
	sum_mol = 0.0;
	sum_wt  = 0.0;
	n = 0;
	for (int i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			sp[0].ph_frac_vol[n] = sp[0].ph_frac_wt[n] / sp[0].SS[n].rho;
			sum_vol += sp[0].ph_frac_vol[n];
			sum_mol += sp[0].ph_frac[n];
			sum_wt  += sp[0].ph_frac_wt[n];
			n+=1;
		}
	}
	m = 0;
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			sp[0].ph_frac_vol[n] =  sp[0].ph_frac_wt[n] / sp[0].PP[m].rho;
			sum_vol += sp[0].ph_frac_vol[n];
			sum_mol += sp[0].ph_frac[n];
			sum_wt  += sp[0].ph_frac_wt[n];
			m +=1;
			n +=1;
		}
	}
	for (int i = 0; i < gv.n_phase; i++){
		sp[0].ph_frac_vol[i] 	/= sum_vol;
		sp[0].ph_frac[i] 		/= sum_mol;
		sp[0].ph_frac_wt[i] 	/= sum_wt;
	}

	/* The following section normalizes the entries for S (solid), M (melt) and F (fluid) which are entries useful for geodynamic coupling */
	// normalize rho_S and bulk_S
	if (sp[0].frac_S > 0.0){
		sp[0].rho_S  				/= sp[0].frac_S;
		for (j = 0; j < gv.len_ox; j++){
			sp[0].bulk_S[j]	   		/= sp[0].frac_S;
		}
	}

	sum = 0.0;
	sum_mol = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_S_wt[j]	    = sp[0].bulk_S[j]*z_b.masspo[j];
		sum 				   += sp[0].bulk_S_wt[j];
		sum_mol 			   += sp[0].bulk_S[j];
	}

	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_S_wt[j]	   /= sum;
		sp[0].bulk_S[j] 	   /= sum_mol;
	}

	atp2wt = sum/sum_Molar_mass_bulk;
	sp[0].frac_S_wt 		    = sp[0].frac_S*atp2wt;

	// normalize bulk_M
	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_M[j]	   		/= sp[0].frac_M;
	}

	sum 	= 0.0;
	sum_mol = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		sum 				   += sp[0].bulk_M_wt[j];
		sum_mol 			   += sp[0].bulk_M[j];
	}

	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_M_wt[j]	   /= sum;
		sp[0].bulk_M[j] 	   /= sum_mol;
	}

	// normalize rho_F and bulk_F
	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_F[j]	   		/= sp[0].frac_F;
	}

	sum 	= 0.0;
	sum_mol = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		sum 				   += sp[0].bulk_F_wt[j];
		sum_mol 			   += sp[0].bulk_F[j];
	}

	for (j = 0; j < gv.len_ox; j++){
		sp[0].bulk_F_wt[j]	   /= sum;
		sp[0].bulk_F[j] 	   /= sum_mol;
	}

	// Normalize the fraction of S + M + F = 1.0
	sum = sp[0].frac_F + sp[0].frac_M + sp[0].frac_S;
	sp[0].frac_F /= sum;
	sp[0].frac_M /= sum;
	sp[0].frac_S /= sum;


	sum = sp[0].frac_F_wt + sp[0].frac_M_wt + sp[0].frac_S_wt;

	sp[0].frac_F_wt /= sum;
	sp[0].frac_M_wt /= sum;
	sp[0].frac_S_wt /= sum;

	/* compute cp as J/K/kg for given bulk-rock composition */
	double MolarMass_system = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		MolarMass_system += z_b.bulk_rock[i]*(z_b.masspo[i]);
	}
	sp[0].s_cp 					= sp[0].cp_wt/MolarMass_system*1e6;


	/* get LP assemblage - This routine retrieves the information of the solution phases and pure phase as computed at equilibrium -> to be used as initial guess */
	m = 0;
	simplex_data *d  = (simplex_data *) splx_data;

	for (i = 0; i < d->n_Ox; i++){
		ph_id 		= d->ph_id_A[i][1];
			
		if (d->ph_id_A[i][0] == 1){			/* if phase is a pure species */
			for (int j = 0; j < gv.len_ox; j++){
				sp[0].mSS[m].comp_Ppc[j] = PP_ref_db[ph_id].Comp[j]*PP_ref_db[ph_id].factor;
			}

			sp[0].mSS[m].G_Ppc 	= PP_ref_db[ph_id].gbase*PP_ref_db[ph_id].factor;
			sp[0].mSS[m].DF_Ppc = sp[0].mSS[m].G_Ppc;
			for (int j = 0; j < gv.len_ox; j++) {
				sp[0].mSS[m].DF_Ppc -= sp[0].mSS[m].comp_Ppc[j] *gv.gam_tot[j];
			}

			strcpy(sp[0].mSS[m].info,"lpig");
			strcpy(sp[0].mSS[m].ph_type,"pp");
			strcpy(sp[0].mSS[m].ph_name,gv.PP_list[ph_id]);
			sp[0].mSS[m].ph_id 		= ph_id;
			sp[0].mSS[m].nOx 		= gv.len_ox;
			sp[0].mSS[m].n_xeos		= 0;
			sp[0].mSS[m].n_em 		= 0;

			sp[0].n_mSS += 1;
			m += 1;
		}
		else if (d->ph_id_A[i][0] == 2){ 		/* pure endmembers as solution phase */
			em_id 					= d->ph_id_A[i][3];

			for (j = 0; j < SS_ref_db[ph_id].n_em; j++) {	
				SS_ref_db[ph_id].p[j] = gv.em2ss_shift;
			}
			SS_ref_db[ph_id].p[em_id] = 1.0 - gv.em2ss_shift*SS_ref_db[ph_id].n_em;
			
			SS_ref_db[ph_id] = P2X(			gv,
											SS_ref_db[ph_id],
											z_b,
											gv.SS_list[ph_id]					);

			/* get unrotated gbase */
			SS_ref_db[ph_id] = non_rot_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			SS_ref_db[ph_id] = PC_function(				gv,
														SS_ref_db[ph_id], 
														z_b,
														gv.SS_list[ph_id] 		);


			for (j = 0; j < gv.len_ox; j++){
				sp[0].mSS[m].comp_Ppc[j] = SS_ref_db[ph_id].ss_comp[j]*SS_ref_db[ph_id].factor;
			}
											
			sp[0].mSS[m].G_Ppc 	= SS_ref_db[ph_id].df;
			sp[0].mSS[m].DF_Ppc = SS_ref_db[ph_id].df;
			for (j = 0; j < gv.len_ox; j++) {
				sp[0].mSS[m].DF_Ppc -= sp[0].mSS[m].comp_Ppc[j]*gv.gam_tot[j];
			}

			strcpy(sp[0].mSS[m].info,"lpig");
			strcpy(sp[0].mSS[m].ph_type,"ss_em");
			strcpy(sp[0].mSS[m].ph_name,gv.SS_list[ph_id]);
			sp[0].mSS[m].ph_id 		= ph_id;
			sp[0].mSS[m].em_id 		= em_id;
			sp[0].mSS[m].nOx 		= gv.len_ox;
			sp[0].mSS[m].n_xeos		= SS_ref_db[ph_id].n_xeos;
			sp[0].mSS[m].n_em 		= SS_ref_db[ph_id].n_em;

			for (j = 0; j < SS_ref_db[ph_id].n_em; j++){
				sp[0].mSS[m].p_Ppc[j] 	= SS_ref_db[ph_id].p[j];
				sp[0].mSS[m].mu_Ppc[j] 	= SS_ref_db[ph_id].mu[j]*SS_ref_db[ph_id].z_em[j];
			}
			for (j = 0; j < SS_ref_db[ph_id].n_xeos; j++){
				sp[0].mSS[m].xeos_Ppc[j] = SS_ref_db[ph_id].iguess[j];	
			}

			sp[0].n_mSS += 1;
			m += 1;
		}
		else if (d->ph_id_A[i][0] == 3){				/* solution phase */
			pc_id 					= d->ph_id_A[i][3];

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
			
			/* get unrotated gbase */
			SS_ref_db[ph_id] = non_rot_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			SS_ref_db[ph_id] = PC_function(				gv,
														SS_ref_db[ph_id], 
														z_b,
														gv.SS_list[ph_id] 		);
											
			for (j = 0; j < gv.len_ox; j++){
				sp[0].mSS[m].comp_Ppc[j] = SS_ref_db[ph_id].ss_comp[j]*SS_ref_db[ph_id].factor;
			}

			sp[0].mSS[m].G_Ppc 	= SS_ref_db[ph_id].df;
			sp[0].mSS[m].DF_Ppc = SS_ref_db[ph_id].df;
			for (j = 0; j < gv.len_ox; j++) {
				sp[0].mSS[m].DF_Ppc -= sp[0].mSS[m].comp_Ppc[j]*gv.gam_tot[j];
			}

			strcpy(sp[0].mSS[m].info,"lpig");
			strcpy(sp[0].mSS[m].ph_type,"ss");
			strcpy(sp[0].mSS[m].ph_name,gv.SS_list[ph_id]);
			sp[0].mSS[m].ph_id 		= ph_id;
			sp[0].mSS[m].nOx 		= gv.len_ox;
			sp[0].mSS[m].n_xeos		= SS_ref_db[ph_id].n_xeos;
			sp[0].mSS[m].n_em 		= SS_ref_db[ph_id].n_em;

			for (int j = 0; j < SS_ref_db[ph_id].n_em; j++){
				sp[0].mSS[m].p_Ppc[j] 	= SS_ref_db[ph_id].p[j];
				sp[0].mSS[m].mu_Ppc[j] 	= SS_ref_db[ph_id].mu[j]*SS_ref_db[ph_id].z_em[j];
			}
			for (int j = 0; j < SS_ref_db[ph_id].n_xeos; j++){
				sp[0].mSS[m].xeos_Ppc[j] = SS_ref_db[ph_id].iguess[j];	
			}

			sp[0].n_mSS += 1;
			m += 1;
		}
	}

	/* copy metastable phases to sb structure */
	int n_xeos, n_em;
	for (int i = 0; i < gv.len_ss; i++){
		if (SS_ref_db[i].ss_flags[0] == 1){

			n_em 	 = SS_ref_db[i].n_em;
			n_xeos 	 = SS_ref_db[i].n_xeos;
			for (int l = 0; l < SS_ref_db[i].tot_Ppc; l++){
				if (SS_ref_db[i].info_Ppc[l] == 9 && m < gv.max_n_mSS){

					sp[0].n_mSS += 1;
					strcpy(sp[0].mSS[m].info,"ppc");
					strcpy(sp[0].mSS[m].ph_type,"ss");
					strcpy(sp[0].mSS[m].ph_name,gv.SS_list[i]);
					sp[0].mSS[m].ph_id 		= i;
					sp[0].mSS[m].nOx 		= gv.len_ox;
					sp[0].mSS[m].n_xeos		= n_xeos;
					sp[0].mSS[m].n_em 		= n_em;

					sp[0].mSS[m].G_Ppc 		= SS_ref_db[i].G_Ppc[l];

					sp[0].mSS[m].DF_Ppc = SS_ref_db[i].G_Ppc[l];
					for (int j = 0; j < gv.len_ox; j++) {
						sp[0].mSS[m].DF_Ppc -= SS_ref_db[i].comp_Ppc[l][j]*gv.gam_tot[j];
					}

					for (int j = 0; j < gv.len_ox; j++){
						sp[0].mSS[m].comp_Ppc[j] = SS_ref_db[i].comp_Ppc[l][j];
					}
					for (int j = 0; j < n_em; j++){
						sp[0].mSS[m].p_Ppc[j] 	= SS_ref_db[i].p_Ppc[l][j];
						sp[0].mSS[m].mu_Ppc[j] 	= SS_ref_db[i].mu_Ppc[l][j];
					}
					for (int j = 0; j < n_xeos; j++){
						sp[0].mSS[m].xeos_Ppc[j] 	= SS_ref_db[i].xeos_Ppc[l][j];
					}
					
					m += 1;
				}
			}

		}
	}
	if (m >= gv.max_n_mSS){
		printf("WARNING: maximum number of metastable pseudocompounds has been reached, increase the value in gss_init_function.c (SP_INIT_function)\n");
	}


	// debug print
	if (1 == 0){
		printf("Phase vol\n");
		for (int m = 0; m < gv.n_phase; m++){
			printf(" %4s %+10f\n",sp[0].ph[m],sp[0].ph_frac_vol[m]);
		}
		printf("\n");
		printf("Phase wt\n");
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
 	thermocalc like output
*/
void output_thermocalc(			global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	FILE *loc_min;
	char out_lm[255];
	
	int rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int i,j,m,n,k,n_ss;

	/* output active phase fraction*/
	if (numprocs==1){	sprintf(out_lm,	"%s_thermocalc_style_output.txt"		,gv.outpath); 	}
	else 			{	sprintf(out_lm,	"%s_thermocalc_style_output.%i.txt"		,gv.outpath, rank); 	}

	loc_min 	= fopen(out_lm, 	"a"); 
	fprintf(loc_min, "============================================================\n");
	
	n_ss = 0;
	for (i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	"%4s ", cp[i].name);
			n_ss += 1;
		}
	}
	for (i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			fprintf(loc_min, 	"%4s ", gv.PP_list[i]);
		}
	}	
	fprintf(loc_min, " {%10.5f, %10.5f} kbar/°C\n\n",z_b.P,z_b.T-273.15);
	
	fprintf(loc_min, "Compositional variables (solution phase):\n");		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < cp[i].n_xeos; j++){
				fprintf(loc_min, 	"%10.5f ", cp[i].xeos[j]);
			}
			for (k = j; k < gv.len_ox; k++){
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
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
			
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < cp[i].n_em; j++){
				fprintf(loc_min, 	"%10.5f ", cp[i].p_em[j]);
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
		}
	}

	fprintf(loc_min, "\nEnd-members PGE expression [exp(-mu/(RT))]:\n");		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			
			fprintf(loc_min, 	" %5s", "");
			for (j = 0; j < cp[i].n_em; j++){
				fprintf(loc_min, 	"%10s ", SS_ref_db[cp[i].id].EM_list[j]);
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
			
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < cp[i].n_em; j++){
				fprintf(loc_min, 	"%10.5f ", cp[i].xi_em[j]);
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
		}
	}

	fprintf(loc_min, "\nEnd-members delta_mu:\n");		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			
			fprintf(loc_min, 	" %5s", "");
			for (j = 0; j < cp[i].n_em; j++){
				fprintf(loc_min, 	"%10s ", SS_ref_db[cp[i].id].EM_list[j]);
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
			
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < cp[i].n_em; j++){
				fprintf(loc_min, 	"%10.5f ", cp[i].mu[j]);
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
		}
	}

	fprintf(loc_min, "\nSite fractions:\n");		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < (cp[i].n_sf); j++){
				fprintf(loc_min, 	"%10.5f ", cp[i].sf[j]); // *-1.0 because inequality are given as -x <= 0 in NLopt
			}
			for (k = j; k < gv.len_ox; k++){
				fprintf(loc_min, 	"%10s ", "-");
			}		
			fprintf(loc_min, "\n");	
		}
	}
	
	fprintf(loc_min, "\nOxide compositions [mol fraction] (normalized on 1 atom basis):\n");	
	fprintf(loc_min, "%5s"," ");
	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, " %10s", gv.ox[i]);
	}
	fprintf(loc_min, "\n %5s","SYS");	
	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, "%10.5f ",z_b.bulk_rock[i]);
	}
	fprintf(loc_min, "\n");	
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < gv.len_ox; j++){
				fprintf(loc_min, "%10.5f ", cp[i].ss_comp[j]*cp[i].factor);
			}
			fprintf(loc_min, "\n");
		}
	}

	for (i = 0; i < gv.len_pp; i++){ 
		if (gv.pp_flags[i][1] == 1){
			fprintf(loc_min, 	" %5s", gv.PP_list[i]);
			for (j = 0; j < gv.len_ox; j++){
				fprintf(loc_min, "%10.5f ", PP_ref_db[i].Comp[j]*PP_ref_db[i].factor);
			}
			fprintf(loc_min, "\n");
		}
	}
	
	double G;
	fprintf(loc_min, "\n");	
	fprintf(loc_min, "Stable mineral assemblage:\n");	
	fprintf(loc_min, "%6s%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","phase","mod[mol fr]","f","G[J]" ,"V[cm3/mol]" ,"Cp[kJ/K]","Rho[kg/m3]","Alpha[1/K]","Entropy[J/K]","Enthalpy[J]","BulkMod[GPa]","ShearMod[GPa]","Vp[km/s]","Vs[km/s]");
				
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			

			G = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}

			fprintf(loc_min, "%6s", cp[i].name);
			fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.6f %+12.6f %+12.2f %+12.2f %+12.2f %+12.2f",
						cp[i].ss_n,cp[i].factor,
						G,
						cp[i].volume*10.,
						cp[i].phase_cp,
						cp[i].phase_density,
						cp[i].phase_expansivity,
						cp[i].phase_entropy,
						cp[i].phase_enthalpy,
						cp[i].phase_bulkModulus/10.,
						cp[i].phase_shearModulus/10.,
						sqrt((cp[i].phase_bulkModulus/10. +4.0/3.0*cp[i].phase_shearModulus/10.)/(cp[i].phase_density/1e3)),
						sqrt(cp[i].phase_shearModulus/10.0/(cp[i].phase_density/1e3)) 
					);
			fprintf(loc_min, "\n");
		}
	}
	
	for (i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){ 
			fprintf(loc_min, "%6s", gv.PP_list[i]);
			fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.6f %+12.6f %+12.2f %+12.2f %+12.2f %+12.2f",
					gv.pp_n[i],
					PP_ref_db[i].factor,
					PP_ref_db[i].gbase,
					PP_ref_db[i].volume*10.,
					PP_ref_db[i].phase_cp,
					PP_ref_db[i].phase_density,
					PP_ref_db[i].phase_expansivity,
					PP_ref_db[i].phase_entropy,
					PP_ref_db[i].phase_enthalpy,
					PP_ref_db[i].phase_bulkModulus/10.,
					PP_ref_db[i].phase_shearModulus/10.,
					sqrt((PP_ref_db[i].phase_bulkModulus/10. +4.0/3.0*PP_ref_db[i].phase_shearModulus/10.)/(PP_ref_db[i].phase_density/1e3)),
					sqrt(PP_ref_db[i].phase_shearModulus/10.0/(PP_ref_db[i].phase_density/1e3))
				);
			fprintf(loc_min, "\n");
		}
	}
	
	G = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		G += z_b.bulk_rock[j]*gv.gam_tot[j];
	}
	fprintf(loc_min, "%6s %24s %+12.5f %+12.5f %+12.5f %+12.5f %12s %+12.6f %+12.6f %+12.5f %+12.5f %+12.5f %+12.5f\n",
			"SYS",
			" ",
			G,
			gv.system_volume*10.,
			gv.system_cp,
			gv.system_density,
			" ",
			gv.system_entropy,
			gv.system_enthalpy,
			gv.system_bulkModulus,
			gv.system_shearModulus,
			gv.system_Vp,
			gv.system_Vs
	);
	
	/* output solution phase composition */ 	
	fprintf(loc_min, "\nGamma[J] (chemical potential of oxides):\n");
	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, "%6s %+12.5f\n", gv.ox[i], gv.gam_tot[i]);
	}
	fprintf(loc_min, "\nG-hyperplane distance[J]:\n");
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	"%5s %+10e\n", cp[i].name,cp[i].df);
		}
	}
	fprintf(loc_min, "\n\n");

	/* save initial guess for THERMOCALC */
	fprintf(loc_min, "Initial guess for THERMOCALC:\n");
	fprintf(loc_min, "%% ----------------------------------------------------------\n");
	fprintf(loc_min, "%% at P =  %12.8f, T = %12.8f, for: ",z_b.P,z_b.T-273.15);
	for (i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	"%s ", cp[i].name);

		}
	}
	fprintf(loc_min, "\n");
	fprintf(loc_min, "%% ----------------------------------------------------------\n");
	fprintf(loc_min, "ptguess  %12.8f %12.8f\n",z_b.P,z_b.T-273.15);
	fprintf(loc_min, "%% ----------------------------------------------------------\n");
	n = 1;
	for (i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			for (j = 0; j < cp[i].n_xeos; j++){
				if (strlen(cp[i].name) == 1){
					fprintf(loc_min, 	"xyzguess %5s(%1s) %10f\n", SS_ref_db[cp[i].id].CV_list[j], cp[i].name, cp[i].xeos[j]);
				}
				else if (strlen(cp[i].name) == 2){
					fprintf(loc_min, 	"xyzguess %5s(%2s) %10f\n", SS_ref_db[cp[i].id].CV_list[j], cp[i].name, cp[i].xeos[j]);
				} 
				else if (strlen(cp[i].name) == 3){
					fprintf(loc_min, 	"xyzguess %5s(%3s) %10f\n", SS_ref_db[cp[i].id].CV_list[j], cp[i].name, cp[i].xeos[j]);
				} 
				else if (strlen(cp[i].name) == 4){
					fprintf(loc_min, 	"xyzguess %5s(%4s) %10f\n", SS_ref_db[cp[i].id].CV_list[j], cp[i].name, cp[i].xeos[j]);
				} 
				else if (strlen(cp[i].name) == 5){
					fprintf(loc_min, 	"xyzguess %5s(%5s) %10f\n", SS_ref_db[cp[i].id].CV_list[j], cp[i].name, cp[i].xeos[j]);
				} 												
			}	
			if (n < n_ss){
				fprintf(loc_min, 	"%% -----------------------------\n");
			}
			n += 1;
		}
	}
	fprintf(loc_min, 	"%% —————————————————————————————\n");

	fclose(loc_min);
}

/**
 	output used for the graphic user interface
*/
void output_gui(				global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	FILE *loc_min;
	char out_lm[255];
	
	int rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int i,j,m,n,k;

	if (numprocs==1){	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath); 		}
	else 			{	sprintf(out_lm,	"%s_pseudosection_output.%i.txt"	,gv.outpath, rank); }
	/* get number of repeated phases for the solvi */
	int n_solvi[gv.len_ss];
	for (i = 0; i < gv.len_ss; i++){
		n_solvi[i] = 0;
	}
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			n_solvi[cp[i].id] += 1;
		}
	}

	loc_min 	= fopen(out_lm, 	"a"); 
	fprintf(loc_min, "%i %i %.10f %.10f %.10f %.10f", gv.numPoint+1, gv.status, z_b.P, z_b.T-273.15, gv.G_system,gv.BR_norm);

	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min," %0.10f", gv.gam_tot[i]);
	}
	for (i = gv.len_ox; i < 11; i++){
		fprintf(loc_min," %0.10f", 0.0);
	}
	fprintf(loc_min, " %.10f %.10f %.10f",gv.system_Vp,gv.system_Vs,gv.system_entropy);

	fprintf(loc_min, "\n");
	m = 0;
	for (i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[1] == 1){
			
			if (n_solvi[cp[i].id] > 1){
				fprintf(loc_min, 	"%s_%d \t %.10f \t %.10f \t", cp[i].name,n_solvi[cp[i].id], cp[i].ss_n, cp[i].phase_density);
			}
			else{
				fprintf(loc_min, 	"%s \t %.10f \t %.10f \t", cp[i].name, cp[i].ss_n, cp[i].phase_density);
			}

			fprintf(loc_min, 	"%d ", cp[i].n_xeos);
			for (j = 0; j < (cp[i].n_xeos); j++){
				fprintf(loc_min, 	"%.10f ", cp[i].xeos[j]);
			}
			for (j = 0; j < (cp[i].n_em); j++){
				fprintf(loc_min, 	"%10s ",  SS_ref_db[cp[i].id].EM_list[j]);	
				fprintf(loc_min, 	"%.10f ", cp[i].p_em[j]);
			}
			fprintf(loc_min, 	"%d ",   gv.len_ox);
			for (int j = 0; j < gv.len_ox; j++){
				fprintf(loc_min, 	"%10s ",   gv.ox[j]);	
				fprintf(loc_min, 	"%.10f ",  sp[0].SS[m].Comp_wt[j]);	
				// fprintf(loc_min, 	"%.10f ",  cp[i].ss_comp[j]*cp[i].factor);	
			}
			fprintf(loc_min, 	"%.10f ",  sp[0].ph_frac_wt[m]);	
			fprintf(loc_min, "\n");
			m += 1;
		}
	}	
	n = 0;
	for (i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			fprintf(loc_min, 	"%s \t %.10f \t %.10f \t", gv.PP_list[i], gv.pp_n[i], PP_ref_db[i].phase_density);
			fprintf(loc_min, 	"%d ",   0);
			fprintf(loc_min, 	"%d ",   gv.len_ox);

			for (int j = 0; j < gv.len_ox; j++){
				fprintf(loc_min, 	"%10s ",   gv.ox[j]);	
				fprintf(loc_min, 	"%.10f ",  sp[0].PP[n].Comp_wt[j]);	
				// fprintf(loc_min, 	"%.10f ",  PP_ref_db[i].Comp[j]*PP_ref_db[i].factor);	
			}
			fprintf(loc_min, 	"%.10f ",  sp[0].ph_frac_wt[m+n]);	
			fprintf(loc_min, "\n");
			n += 1;
		}
	}	
	fprintf(loc_min, "\n");
	fclose(loc_min);


	/**
	 * @brief output parameters to calculate seismic wave velocity correction
	 * 
	 */
	// FILE *tot_min;
	// char tot_lm[255];
	
	// if (numprocs==1){	sprintf(tot_lm,	"%s_wave_output.txt"		,gv.outpath); 		}
	// else 			{	sprintf(tot_lm,	"%s_wave_output.%i.txt"		,gv.outpath, rank); }

	// tot_min 	= fopen(tot_lm, 	"a"); 
	
	// if (sp[0].frac_M < 1.0){
	// 	if (sp[0].frac_M > 0.0){
	// 	fprintf(tot_min, "%i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n",gv.numPoint+1, z_b.P, z_b.T-273.15,sp[0].Vs_S,sp[0].Vp_S,sp[0].bulkModulus_S,sp[0].shearModulus_S,sp[0].bulkModulus_M,sp[0].rho_M,sp[0].rho_S,sp[0].frac_M,
	// 																										sp[0].bulk_S_wt[5]+sp[0].bulk_S_wt[6],sp[0].bulk_M_wt[5]+sp[0].bulk_M_wt[6], sp[0].bulk_S_wt[0], sp[0].bulk_M_wt[0], sp[0].bulk_M_wt[10]);
	// 	}
	// 	else{
	// 	fprintf(tot_min, "%i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n", gv.numPoint+1, z_b.P, z_b.T-273.15,sp[0].Vs_S,sp[0].Vp_S,sp[0].bulkModulus_S,sp[0].shearModulus_S,sp[0].bulkModulus_M,sp[0].rho_M,sp[0].rho_S,sp[0].frac_M,
	// 																										sp[0].bulk_S_wt[5]+sp[0].bulk_S_wt[6],0.0/0.0, sp[0].bulk_S_wt[0], 0.0/0.0, 0.0/0.0);
	// 	}
	// }
	// else{
	// 	fprintf(tot_min, "%i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n", gv.numPoint+1, z_b.P, z_b.T-273.15, 0.0/0.0, 0.0/0.0, 0.0/0.0, 0.0/0.0, sp[0].bulkModulus_M, 0.0/0.0, 0.0/0.0,sp[0].frac_M,
	// 																										0.0/0.0,sp[0].bulk_M_wt[5]+sp[0].bulk_M_wt[6],0.0/0.0, sp[0].bulk_M_wt[0], sp[0].bulk_M_wt[10]);
	// }
	
	// fclose(tot_min);

}


/**
 	output used for the matlab interface
*/
void output_matlab(				global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	FILE *loc_min;
	char out_lm[255];
	
	int rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int i,j,m,n,k;

	/* output active phase fraction*/
	if (numprocs==1){	sprintf(out_lm,	"%s_matlab_output.txt"		,gv.outpath); 	}
	else 			{	sprintf(out_lm,	"%s_matlab_output.%i.txt"		,gv.outpath, rank); 	}

	loc_min 	= fopen(out_lm, 	"a"); 
	fprintf(loc_min, "============================================================\n");
	
	for (i = 0; i < gv.len_cp; i++){
		if ( cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	"%4s ", cp[i].name);
		}
	}
	for (i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			fprintf(loc_min, 	"%4s ", gv.PP_list[i]);
		}
	}	
	fprintf(loc_min, " {%10.5f, %10.5f} kbar/°C\n\n",z_b.P,z_b.T-273.15);
	
	fprintf(loc_min, "\nEnd-members fractions[wt fr]\n");	

	for (m = 0; m < gv.n_cp_phase; m++){
		fprintf(loc_min, 	" %5s", "");
		for (j = 0; j < sp[0].SS[m].n_em; j++){
			fprintf(loc_min, 	"%10s ", sp[0].SS[m].emNames[j]);
		}
		for (k = j; k < gv.len_ox; k++){
			fprintf(loc_min, 	"%10s ", "-");
		}	
		fprintf(loc_min, "\n");

		fprintf(loc_min, 	" %5s", sp[0].ph[m]);
		for (j = 0; j < sp[0].SS[m].n_em; j++){
			fprintf(loc_min, 	"%10.5f ", sp[0].SS[m].emFrac_wt[j]);
		}
		for (k = j; k < gv.len_ox; k++){
			fprintf(loc_min, 	"%10s ", "-");
		}		
		fprintf(loc_min, "\n");	
	}

	fprintf(loc_min, "\n\nSite fractions:\n");		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	" %5s", cp[i].name);
			for (j = 0; j < (cp[i].n_sf); j++){
				fprintf(loc_min, 	"%8.5f ", cp[i].sf[j]); // *-1.0 because inequality are given as -x <= 0 in NLopt
			}
			for (k = j; k < 18; k++){
				fprintf(loc_min, 	"%8s ", "-");
			}		
			fprintf(loc_min, "\n");	
		}
	}
	
	fprintf(loc_min, "\n\nOxide compositions [wt fr]:\n");	
	fprintf(loc_min, "%5s"," ");
	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, " %10s", gv.ox[i]);
	}
	fprintf(loc_min, "\n %5s","SYS");
	for (int i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, "%10.5f ",sp[0].bulk_wt[i]);
	}
	fprintf(loc_min, "\n");	

	for (int m = 0; m < gv.n_cp_phase; m++){
		fprintf(loc_min, 	" %5s", sp[0].ph[m]);
		for (int j = 0; j < gv.len_ox; j++){
			fprintf(loc_min, "%10.5f ", sp[0].SS[m].Comp_wt[j]);
		}
		fprintf(loc_min, "\n");
	}
	n = 0;
	for (m = gv.n_cp_phase; m < gv.n_phase; m++){
		fprintf(loc_min, 	" %5s", sp[0].ph[m]);
		for (j = 0; j < gv.len_ox; j++){
			fprintf(loc_min, "%10.5f ", sp[0].PP[n].Comp_wt[j]);
		}
		fprintf(loc_min, "\n");
		n += 1;
	}

	// fprintf(loc_min, "\n\nEnd-members compositions[wt fr]\n");	
	// fprintf(loc_min, "%5s %5s", "SS", "EM");
	// for (i = 0; i < gv.len_ox; i++){
	// 	fprintf(loc_min, " %10s", gv.ox[i]);
	// }
	// fprintf(loc_min, "\n");
	// for (int m = 0; m < gv.n_cp_phase; m++){
	// 	for (j = 0; j < sp[0].SS[m].n_em; j++){
	// 		fprintf(loc_min, 	"%5s ", sp[0].ph[m]);
	// 		fprintf(loc_min, 	"%5s ", sp[0].SS[m].emNames[j]);
	// 		for (int k = 0; k < gv.len_ox; k++){
	// 			fprintf(loc_min, 	"%10.5f ", sp[0].SS[m].emComp_wt[j][k]);
	// 		}	
	// 		fprintf(loc_min, "\n");
	// 	}
	// 	fprintf(loc_min, "\n");
	// }

	double G;
	fprintf(loc_min, "\n\n");	
	fprintf(loc_min, "Stable mineral assemblage:\n");	
	fprintf(loc_min, "%6s%15s %13s %17s %17s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","phase","fraction[wt]","G[J]" ,"V_molar[cm3/mol]","V_partial[cm3]" ,"Cp[kJ/K]","Rho[kg/m3]","Alpha[1/K]","Entropy[J/K]","Enthalpy[J]","BulkMod[GPa]","ShearMod[GPa]","Vp[km/s]","Vs[km/s]");

	n = 0;		
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){

			G = 0.0;
			for (j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}

			fprintf(loc_min, "%6s", cp[i].name);
			fprintf(loc_min, "%+15.5f %+13.5f %+17.5f %+17.5f %+12.5f %+12.5f %+12.8f %+12.6f %+14.4f %+12.2f %+12.2f %+12.2f %+12.2f",
						sp[0].ph_frac_wt[n],
						G,
						cp[i].volume*10.,
						cp[i].volume*10.*cp[i].ss_n_mol*cp[i].factor,
						cp[i].phase_cp,
						cp[i].phase_density,
						cp[i].phase_expansivity,
						cp[i].phase_entropy,
						cp[i].phase_enthalpy,
						cp[i].phase_bulkModulus/10.,
						cp[i].phase_shearModulus/10.,
						sqrt((cp[i].phase_bulkModulus/10. +4.0/3.0*cp[i].phase_shearModulus/10.)/(cp[i].phase_density/1e3)),
						sqrt(cp[i].phase_shearModulus/10.0/(cp[i].phase_density/1e3)) 
					);
			fprintf(loc_min, "\n");
			n += 1;
		}
	}
	int n_ss = n;
	for (i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){ 
			fprintf(loc_min, "%6s", gv.PP_list[i]);
			fprintf(loc_min, "%+15.5f %+13.5f %+17.5f %+17.5f %+12.5f %+12.5f %+12.8f %+12.6f %+14.4f %+12.2f %+12.2f %+12.2f %+12.2f",
					sp[0].ph_frac_wt[n],
					PP_ref_db[i].gbase,
					PP_ref_db[i].volume*10.,
					PP_ref_db[i].volume*10.*gv.pp_n[i]*PP_ref_db[i].factor,
					PP_ref_db[i].phase_cp,
					PP_ref_db[i].phase_density,
					PP_ref_db[i].phase_expansivity,
					PP_ref_db[i].phase_entropy,
					PP_ref_db[i].phase_enthalpy,
					PP_ref_db[i].phase_bulkModulus/10.,
					PP_ref_db[i].phase_shearModulus/10.,
					sqrt((PP_ref_db[i].phase_bulkModulus/10. +4.0/3.0*PP_ref_db[i].phase_shearModulus/10.)/(PP_ref_db[i].phase_density/1e3)),
					sqrt(PP_ref_db[i].phase_shearModulus/10.0/(PP_ref_db[i].phase_density/1e3))
				);
			fprintf(loc_min, "\n");
			n+=1;
		}
	}
	
	G = 0.0;
	for (j = 0; j < gv.len_ox; j++){
		G += z_b.bulk_rock[j]*gv.gam_tot[j];
	}
	fprintf(loc_min, "%6s %14s %+13.5f %17s %+17.5f %+12.5f %+12.5f %12s %+12.6f %+14.4f %+12.5f %+12.5f %+12.5f %+12.5f\n",
			"SYS",
			" ",
			G,
			" ",
			gv.system_volume*10.,
			gv.system_cp,
			gv.system_density,
			" ",
			gv.system_entropy,
			gv.system_enthalpy,
			gv.system_bulkModulus,
			gv.system_shearModulus,
			gv.system_Vp,
			gv.system_Vs
	);
	
	/* output solution phase composition */ 	
	fprintf(loc_min, "\n\nGamma[J] (chemical potential of oxides):\n");
	for (i = 0; i < gv.len_ox; i++){
		fprintf(loc_min, "%6s %+12.5f\n", gv.ox[i], gv.gam_tot[i]);
	}

	fprintf(loc_min, "\n\nSystem fugacity:\n");
	fprintf(loc_min, 	"%6s %+10e\n", "log10(fO2)",log10(gv.system_fO2));
	fprintf(loc_min, 	"%6s %+10e\n", "log10(dQFM)",log10(gv.system_deltaQFM));

	fprintf(loc_min, "\n\nSystem activity:\n");
	fprintf(loc_min, 	"%6s %+10e\n", "aH2O",gv.system_aH2O);
	fprintf(loc_min, 	"%6s %+10e\n", "aSiO2",gv.system_aSiO2);
	fprintf(loc_min, 	"%6s %+10e\n", "aTiO2",gv.system_aTiO2);
	fprintf(loc_min, 	"%6s %+10e\n", "aAl2O3",gv.system_aAl2O3);
	fprintf(loc_min, 	"%6s %+10e\n", "aMgO",gv.system_aMgO);
	fprintf(loc_min, 	"%6s %+10e\n", "aFeO",gv.system_aFeO);

	fprintf(loc_min, "\n\nG-hyperplane distance[J]:\n");
	for (i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			fprintf(loc_min, 	"%5s %+10e\n", cp[i].name,cp[i].df);
		}
	}
	fprintf(loc_min, "\n\n");

	/* save initial guess for THERMOCALC */
	if (gv.output_matlab == 2){
		fprintf(loc_min, "Initial guess for THERMOCALC\n");
		fprintf(loc_min, "%% ----------------------------------------------------------\n");
		fprintf(loc_min, "%% at P = %10f, T = %10f, for: ",z_b.P,z_b.T-273.15);
		for (i = 0; i < gv.len_cp; i++){
			if ( cp[i].ss_flags[1] == 1){
				if (strlen(cp[i].name) == 1){
					fprintf(loc_min, 	"%1s ", cp[i].name);
				}
				else if (strlen(cp[i].name) == 2){
					fprintf(loc_min, 	"%2s ", cp[i].name);
				}
				else if (strlen(cp[i].name) == 3){
					fprintf(loc_min, 	"%3s ", cp[i].name);
				}
				else if (strlen(cp[i].name) == 4){
					fprintf(loc_min, 	"%4s ", cp[i].name);
				}
				else if (strlen(cp[i].name) == 5){
					fprintf(loc_min, 	"%5s ", cp[i].name);
				}
			}
		}
		for (i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){ 
				if (strlen(gv.PP_list[i]) == 1){
					fprintf(loc_min, "%1s ", gv.PP_list[i]);
				}
				else if (strlen(gv.PP_list[i]) == 2){
					fprintf(loc_min, "%2s ", gv.PP_list[i]);
				}
				else if (strlen(gv.PP_list[i]) == 3){
					fprintf(loc_min, "%3s ", gv.PP_list[i]);
				}
				else if (strlen(gv.PP_list[i]) == 4){
					fprintf(loc_min, "%4s ", gv.PP_list[i]);
				}
				else if (strlen(gv.PP_list[i]) == 5){
					fprintf(loc_min, "%5s ", gv.PP_list[i]);
				}
			}
		}

		fprintf(loc_min, "\n");
		fprintf(loc_min, "%% ----------------------------------------------------------\n");
		fprintf(loc_min, "ptguess %10f %10f\n",z_b.P,z_b.T-273.15);
		fprintf(loc_min, "%% ----------------------------------------------------------\n");
		n = 0;
		for (i = 0; i < gv.len_cp; i++){
			if ( cp[i].ss_flags[1] == 1){
				n +=1;
				for (j = 0; j < SS_ref_db[cp[i].id].n_xeos; j++){
					if (strlen(cp[i].name) == 1){
						fprintf(loc_min, 	"xyzguess %5s(%1s) %10f\n", SS_ref_db[cp[i].id].CV_list[j],cp[i].name, cp[i].xeos[j]);
					}
					else if (strlen(cp[i].name) == 2){
						fprintf(loc_min, 	"xyzguess %5s(%2s) %10f\n", SS_ref_db[cp[i].id].CV_list[j],cp[i].name, cp[i].xeos[j]);
					}
					else if (strlen(cp[i].name) == 3){
						fprintf(loc_min, 	"xyzguess %5s(%3s) %10f\n", SS_ref_db[cp[i].id].CV_list[j],cp[i].name, cp[i].xeos[j]);
					}
					else if (strlen(cp[i].name) == 4){
						fprintf(loc_min, 	"xyzguess %5s(%4s) %10f\n", SS_ref_db[cp[i].id].CV_list[j],cp[i].name, cp[i].xeos[j]);
					}
					else if (strlen(cp[i].name) == 5){
						fprintf(loc_min, 	"xyzguess %5s(%5s) %10f\n", SS_ref_db[cp[i].id].CV_list[j],cp[i].name, cp[i].xeos[j]);
					}
				}
				if (n < n_ss){
					fprintf(loc_min, "%% -----------------------------\n");
				}
				
			}
		}
		fprintf(loc_min, "%% —————————————————————————————\n");
	}


	fclose(loc_min);	
}



/**
  Save final result of minimization
*/
void save_results_function(		global_variable 	 gv,
								bulk_info 	 		 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
){
	
	FILE *loc_min;
	char out_lm[255];
	
	int i,j, rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (gv.output_matlab >= 1){
		output_matlab(					gv,											/** global variables (e.g. Gamma) 	*/
										z_b,										/** bulk-rock informations 			*/
										PP_ref_db,									/** pure phase database 			*/
										SS_ref_db,									/** solution phase database 		*/
										cp,
										sp						);
	}

	if (gv.verbose == 1 && gv.output_matlab == 0){
		output_thermocalc(				gv,											/** global variables (e.g. Gamma) 	*/
										z_b,										/** bulk-rock informations 			*/
										PP_ref_db,									/** pure phase database 			*/
										SS_ref_db,									/** solution phase database 		*/
										cp,
										sp						);
	}
	
	/** ----------------------------------------------------------------------------------------------- **/
	/** MATLAB GRID OUTPUT **/
	if (gv.verbose == 0){
		output_gui(						gv,											/** global variables (e.g. Gamma) 	*/
										z_b,										/** bulk-rock informations 			*/
										PP_ref_db,									/** pure phase database 			*/
										SS_ref_db,									/** solution phase database 		*/
										cp,
										sp						);
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


// 	char tot_out_lm[255];
// 	char tot_in_lm[255];

// 	if (numprocs == 1){ return; }

// 	sprintf(tot_out_lm,	"%s_wave_output.txt"		,gv.outpath);
//    	FILE *fp2a = fopen(tot_out_lm, "w"); 

// 	fprintf(fp2a, "Number P[kbar]\t T[C]\t Vs0[km/s]\t Vp0[km/s]\t Kb_S[GPa]\t Ks_S[GPa]\t Kb_L[GPa]\t rhoL[kg/m3]\t rhoS[kg/m3]\t frac_melt\t NaK_S[wt]\t NaK_M[wt]\t Si_S[wt]\t Si_M[wt]\t H_M[wt]\n");

// 	// Open file to be merged 
// 	for (i = 0; i < numprocs; i++){
// 		// open file
// 		sprintf(tot_in_lm,	"%s_wave_output.%i.txt"		,gv.outpath, i);
// 		FILE *fp1a = fopen(tot_in_lm, "r"); 
			
// 		// Copy contents of first file to file3.txt 
// 		while ((c = fgetc(fp1a)) != EOF){ 
// 			fputc(c, fp2a); 
// 		}
// 		fclose(fp1a); 
// 	}
//    fclose(fp2a);

}


/**
  Parallel file dump for phase diagrams
*/
void mergeParallel_matlab(global_variable gv){

	int i, rank, numprocs,  MAX_LINE_LENGTH=200;
	char out_lm[255];
	char in_lm[255];
	char c; 
	char buf[MAX_LINE_LENGTH];
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (numprocs == 1){ return; }

	sprintf(out_lm,	"%s_matlab_output.txt"		,gv.outpath);
   	FILE *fp2 = fopen(out_lm, "w"); 

	// Open file to be merged 
	for (i = 0; i < numprocs; i++){
		// open file
		sprintf(in_lm,	"%s_matlab_output.%i.txt"		,gv.outpath, i);
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

	fprintf(fp2, "// BULK-ROCK[len_ox]\tP[kbar]\tT[°C]\tGAMMA[G]\n");

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

