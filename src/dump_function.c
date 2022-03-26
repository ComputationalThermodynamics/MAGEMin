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
	if (gv.verbose != 2){
		/** MATLAB GRID OUTPUT **/
		if (numprocs==1){	sprintf(out_lm,	"%s_pseudosection_output.txt"		,gv.outpath); 		}
		else 			{	sprintf(out_lm,	"%s_pseudosection_output.%i.txt"	,gv.outpath, rank); }
		loc_min 	= fopen(out_lm, 	"w"); 
		fprintf(loc_min, "// NUMBER\t	STATUS[S,R1,R2,F]\tP[kbar]\tT[C]\tG_sys[G]\tbr_norm[wt]\tGAMMA[G]; PHASE[name]\tMODE[wt]\tRHO[kg.m-3]\tX-EOS\n");
		fclose(loc_min);	
			
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
			fprintf(loc_min, "// BULK-ROCK[len_ox]\tP[kbar]\tT[°C]\tGAMMA[G]\n");

			fclose(loc_min);	
		}
	}
}

/**
  Save final result of minimization
*/
void dump_results_function(		global_variable 	gv,
								struct bulk_info 	z_b,

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
		fprintf(loc_min, " {%10.5f, %10.5f} kbar/°C\n\n",z_b.P,z_b.T-273.15);
		
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
		
		fprintf(loc_min, "\nOxide compositions [mol%%] (normalized):\n");	
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
					fprintf(loc_min, "%10.5f ", cp[i].ss_comp[j]);
				}
				fprintf(loc_min, "\n");
			}
		}
		
		double G;
		fprintf(loc_min, "\n");	
		fprintf(loc_min, "Stable mineral assemblage:\n");	
		fprintf(loc_min, "%6s%12s %12s %12s %12s %12s %12s %12s %12s\n","phase","mode","f","G" ,"V" ,"Cp","rho","Thermal_Exp","ShearModu");
					
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
				fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.2f",cp[i].ss_n,cp[i].factor,G,cp[i].volume,cp[i].phase_cp,cp[i].phase_density,cp[i].phase_expansivity,cp[i].phase_shearModulus);
				fprintf(loc_min, "\n");
			}
		}
		
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){ 
				fprintf(loc_min, "%6s", gv.PP_list[i]);
				fprintf(loc_min, "%+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.5f %+12.8f %+12.2f",gv.pp_n[i],PP_ref_db[i].factor,PP_ref_db[i].gbase,PP_ref_db[i].volume,PP_ref_db[i].phase_cp,PP_ref_db[i].phase_density,PP_ref_db[i].phase_expansivity,PP_ref_db[i].phase_shearModulus);
				fprintf(loc_min, "\n");
			}
		}
		
		
			
		G = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			G += z_b.bulk_rock[j]*gv.gam_tot[j];
		}
		fprintf(loc_min, "%6s %24s %+12.5f\n","SYS"," ",G);

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
	if (gv.verbose != 2){
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
		fprintf(loc_min, "%i %i %.10f %.10f %.10f %.10f", gv.numPoint+1, gv.status, z_b.P,z_b.T-273.15,gv.G_system,gv.BR_norm);
		for (i = 0; i < gv.len_ox; i++){
			fprintf(loc_min, 	" \t %0.10f ", gv.gam_tot[i]);
		}
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

	fprintf(fp2, "// NUMBER\tSTATUS[S,R1,R2,F]\tP[kbar]\tT[C]\tG_sys[G]\tbr_norm[wt]\tGAMMA[G]; PHASE[name]\tMODE[wt]\tRHO[kg.m-3]\tX-EOS\n");

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

