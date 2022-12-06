/**
       Input/Output function, to read and save data        
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "MAGEMin.h"

/** 
  read in input data from file 
*/
void read_in_data(		global_variable 	 gv,
						io_data 			*input_data,												/** input data structure */
						int      			 n_points
){
	char line[1000];
	FILE* input_file = fopen(gv.File,"rt");	
	if (gv.File != NULL && input_file != NULL){						/** if input file is provided and exists */
		int k = 0;
		int l = 0;
		/* loop through all the lines of the input file making sure that the number of points is not exceeded */
		while( fgets(line, sizeof(line), input_file) != NULL && k < n_points){
			/* if this is the first line belonging to a PT point to take into account */
			if (l == 0){
				/* first allocate memory to fill gamma array */
				input_data[k].in_bulk      = malloc (gv.len_ox * sizeof (double) ); 
				for (int z = 0; z < gv.len_ox; z++){
					input_data[k].in_bulk[z] = 0.0; 
				}

				sscanf(line, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
					&input_data[k].n_phase, 
					&input_data[k].P, 
					&input_data[k].T,
					&input_data[k].in_bulk[0],
					&input_data[k].in_bulk[1],
					&input_data[k].in_bulk[2],
					&input_data[k].in_bulk[3],
					&input_data[k].in_bulk[4],
					&input_data[k].in_bulk[5],
					&input_data[k].in_bulk[6],
					&input_data[k].in_bulk[7],
					&input_data[k].in_bulk[8],
					&input_data[k].in_bulk[9],
					&input_data[k].in_bulk[10]	);
				
				/* allocate memory depending on the number of provided solution phases */
				input_data[k].phase_names = malloc(input_data[k].n_phase * sizeof(char*));
				for (int i = 0; i < input_data[k].n_phase; i++){
					input_data[k].phase_names[i] = malloc(20 * sizeof(char));
				}
				
				/* allocate memory for compositional variables */
				//input_data[k].sum_phase_xeos = malloc(input_data[k].n_phase * sizeof(double));
				input_data[k].phase_xeos 	 = malloc(input_data[k].n_phase * sizeof(double*));
				for (int i = 0; i < input_data[k].n_phase; i++){
					input_data[k].phase_xeos[i] = malloc((gv.len_ox) * sizeof(double));
				}
				/* initialize x-eos to zeros in case there is mistake in the input file */
				for (int i = 0; i < input_data[k].n_phase; i++){
					for (int j = 0; j < (gv.len_ox); j++){
						input_data[k].phase_xeos[i][j] = gv.bnd_val;
					}
				}
				
				/* allocate memory for endmember fractions */
				//input_data[k].sum_phase_emp = malloc(input_data[k].n_phase * sizeof(double));
				input_data[k].phase_emp 	= malloc(input_data[k].n_phase * sizeof(double*));
				for (int i = 0; i < input_data[k].n_phase; i++){
					input_data[k].phase_emp[i] = malloc((gv.len_ox+1) * sizeof(double));
				}
				/* initialize x-eos to zeros in case there is mistake in the input file */
				for (int i = 0; i < input_data[k].n_phase; i++){
					for (int j = 0; j < (gv.len_ox+1); j++){
						input_data[k].phase_emp[i][j] = 0.0;
					}
				}
			}
			
			/* Lines belonging to the provided x-eos for each solution phase listed */
			/* allocates memory only if the number of phases is not 0 				*/
			if(l > 0 && l < input_data[k].n_phase+1){
				sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
					 input_data[k].phase_names[l-1], 
					&input_data[k].phase_xeos[l-1][0],
					&input_data[k].phase_xeos[l-1][1],
					&input_data[k].phase_xeos[l-1][2],
					&input_data[k].phase_xeos[l-1][3],
					&input_data[k].phase_xeos[l-1][4],
					&input_data[k].phase_xeos[l-1][5],
					&input_data[k].phase_xeos[l-1][6],
					&input_data[k].phase_xeos[l-1][7],
					&input_data[k].phase_xeos[l-1][8],
					&input_data[k].phase_xeos[l-1][9],
					&input_data[k].phase_xeos[l-1][10],

					&input_data[k].phase_emp[l-1][0],
					&input_data[k].phase_emp[l-1][1],
					&input_data[k].phase_emp[l-1][2],
					&input_data[k].phase_emp[l-1][3],
					&input_data[k].phase_emp[l-1][4],
					&input_data[k].phase_emp[l-1][5],
					&input_data[k].phase_emp[l-1][6],
					&input_data[k].phase_emp[l-1][7],
					&input_data[k].phase_emp[l-1][8],
					&input_data[k].phase_emp[l-1][9],
					&input_data[k].phase_emp[l-1][10],
					&input_data[k].phase_emp[l-1][11]	);	
			}
			
			l++;
			/* reset line count when reading a phase is over and increase the phase count */
			if(l > input_data[k].n_phase){l = 0; k += 1;}
		}
		fclose(input_file);
	}
};


/* 
	This initializes the output structure
*/
out_data InitializeOutput(		global_variable gv,
								Databases 		DB)
{
	int 		max_num_EM, i, j;
	out_data 	output;

	max_num_EM 			=	15;
	output.max_num_EM 	= max_num_EM;
	output.Gamma 		= 	malloc ((gv.len_ox) * sizeof (double) ); 
	
	// Determine number of stable SS & PP and allocate structures
	int num_ss, num_pp, num_tot;
	
	num_ss	= 0; 
	num_pp	= 0;
	
	for (int i = 0; i < gv.len_cp; i++){ 
		if (DB.cp[i].ss_flags[1] == 1){
			num_ss += 1;
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			num_pp += 1;
		}
	}		
	num_tot 				= 	num_ss + num_pp;		// total number of phases
	output.n_SS 			= 	num_ss;
	output.n_PP 			= 	num_pp;

	// Allocate structures required later
	output.StableFractions 	= 	malloc((num_tot) * sizeof (double) 		);
	output.Phasedensity  	= 	malloc((num_tot) * sizeof (double) 		);
	output.n_em    			= 	malloc((num_ss)  * sizeof (int) 		);

	output.StableSolutions 	= 	malloc((num_tot) * sizeof(char*)		);
	for (int i = 0; i < num_tot; i++){
		output.StableSolutions[i] = malloc(20 * sizeof(char));
	}

	// Add info about xEOS and p_EM for each of the stable solution models
	output.max_num_EM = max_num_EM;			// ensures that julia knows the correct size
	output.xEOS  	  = malloc (num_ss 	* sizeof (double*)); 
	output.p_EM   	  = malloc (num_ss 	* sizeof (double*)); 
	for (i = 0; i < num_ss; i++){
		output.xEOS[i] 	= malloc ((max_num_EM-1)* sizeof (double) );
		output.p_EM[i] 	= malloc ((max_num_EM  )* sizeof (double) );
	}
	for (i = 0; i < num_ss; i++){
		for (j = 0; j < max_num_EM; j++){
			output.xEOS[i][j] = 0.0;
			output.p_EM[i][j] = 0.0;
		}
	}

	return output;
}

/* 
	This initializes the output structure
*/
void FreeOutput(out_data output)
{
	int max_num_EM, i, j, num_tot;

	max_num_EM 		=	15;
	num_tot 		= 	output.n_SS + output.n_PP;		// total number of phases

	// Allocate structures required later
	free(output.StableFractions);
	free(output.Phasedensity);
	free(output.n_em);
	
	for (int i = 0; i < num_tot; i++){
		free(output.StableSolutions[i]);
	}
	free(output.StableSolutions);

	// Add info about xEOS and p_EM for each of the stable solution models
	for (i = 0; i < output.n_SS; i++){
		free(output.xEOS[i]);
		free(output.p_EM[i]);
	}
	free(output.xEOS); 
	free(output.p_EM); 

}

/**
 * This adds the results to an output struct, which can be passed on to julia.
 * 
 * NOTE: we likely need a cleanup struct for this as well, as we do malloc and julia will likely not release 
**/
void AddResults_output_struct(		global_variable 	gv,
									bulk_info 			z_b,
									Databases 			DB, 
									out_data 			output		)
{
	int i,j,num;

	printf("\n ********* Outputting data: P=%f \n",z_b.P);
	output.status 	= 	gv.status;
	output.number   = 	gv.numPoint+1;
	output.P 		= 	z_b.P;
	output.T 		= 	z_b.T-273.15;
	output.G_system =	gv.G_system;
	output.BR_norm  =	gv.BR_norm;

	// Gamma
	for (i = 0; i < gv.len_ox; i++){
		output.Gamma[i] = gv.gam_tot[i];
	}

	// Copy the names & fractions of solid solutions
	num = 0; 
	for (int i = 0; i < gv.len_cp; i++){ 
		if (DB.cp[i].ss_flags[1] == 1){
	
		//	strcpy(output.StableSolutions[num],gv.SS_list[i]);
			output.StableFractions[num] = 	DB.cp[i].ss_n;
			output.Phasedensity[num] 	=	DB.cp[i].phase_density;
			output.n_em[num]			=	DB.cp[i].n_em;

			num += 1;
		}
	}

	// Copy the names & fractions of pure phases (if present)
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			strcpy(output.StableSolutions[num],gv.PP_list[i]);
			output.StableFractions[num] = 	gv.pp_n[i];
			output.Phasedensity[num] 	=	DB.PP_ref_db[i].phase_density;

			num += 1;
		}
	}

#if 1
	
	// Add xEOS and p_EM info to the output struct
	num = 0;
	for (int i = 0; i < gv.len_cp; i++){ 
		if (DB.cp[i].ss_flags[1] == 1){
			for (int j = 0; j < (output.n_em[num]-1); j++){
	//			output.xEOS[j][num] = DB.SS_ref_db[i].xeos[j];
			}
			for (int j = 0; j < (output.n_em[num]); j++){
	//			output.p_EM[num][j] = DB.SS_ref_db[i].p[j];
			}
			num += 1;
		}
	}


	printf("# of stable SS=%i PP=%i \n",output.n_SS, output.n_PP);

#endif

}
