/**
       Input/Output function, to read and save data        
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "MAGEMin.h"


#define n_ox 11
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
				input_data[k].in_bulk      = malloc (n_ox * sizeof (double) ); 
				for (int z = 0; z < n_ox; z++){
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
					input_data[k].phase_xeos[i] = malloc((n_ox) * sizeof(double));
				}
				/* initialize x-eos to zeros in case there is mistake in the input file */
				for (int i = 0; i < input_data[k].n_phase; i++){
					for (int j = 0; j < (n_ox); j++){
						input_data[k].phase_xeos[i][j] = gv.bnd_val;
					}
				}
				
				/* allocate memory for endmember fractions */
				//input_data[k].sum_phase_emp = malloc(input_data[k].n_phase * sizeof(double));
				input_data[k].phase_emp 	= malloc(input_data[k].n_phase * sizeof(double*));
				for (int i = 0; i < input_data[k].n_phase; i++){
					input_data[k].phase_emp[i] = malloc((n_ox+1) * sizeof(double));
				}
				/* initialize x-eos to zeros in case there is mistake in the input file */
				for (int i = 0; i < input_data[k].n_phase; i++){
					for (int j = 0; j < (n_ox+1); j++){
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


