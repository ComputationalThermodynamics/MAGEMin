/**
        Pure phase minimization function                    
-----------------------------------------------------------

This function simply update the driving forces of pure phase each time the G-hyperplane is tilted.
                  
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "MAGEMin.h"
#include "gem_function.h"


/**
  main pure phase minimization routine
*/
void pp_min_function(		global_variable 	 gv,
							bulk_info 	 		 z_b,
							PP_ref 				*PP_ref_db
){
	// update delta_G of pure phases using Gamma
	for (int i = 0; i < gv.len_pp; i++){
		/* if pure phase is active or on hold (PP cannot be removed from consideration */
		if (gv.pp_flags[i][0] == 1){
			PP_ref_db[i].gb_lvl = PP_ref_db[i].gbase;
			
			/* level the phase using chemical potential of oxides (gamma) */
			for (int j = 0; j < gv.len_ox; j++) {
				PP_ref_db[i].gb_lvl -= PP_ref_db[i].Comp[j]*gv.gam_tot[j];
			}
			
			gv.pp_xi[i] = exp(-PP_ref_db[i].gb_lvl/(z_b.R*z_b.T));
		}
	}
};


/** 
  initialize pure phase database */
global_variable init_em_db(		int 				EM_database,
								bulk_info 			z_b,
								global_variable 	gv,
								PP_ref 			   *PP_ref_db
){

		/* initialize endmember database */
		char state[] = "equilibrium";	
		int sum_zel;
		for (int i = 0; i < gv.len_pp; i++){

			if (strcmp( gv.PP_list[i], "qfm") == 0){

				PP_ref q 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"q", 
												state				);
				PP_ref fa 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"fa", 
												state				);

				PP_ref mt 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = -3.0 * fa.Comp[j] + 3.0*q.Comp[j] + 2.0*mt.Comp[j];
				}								
				PP_ref_db[i].gbase   =  -3.0 * fa.gbase + 3.0*q.gbase + 2.0*mt.gbase + z_b.T*0.019145*gv.QFM_n;
				PP_ref_db[i].factor  =  -3.0 * fa.factor + 3.0*q.factor + 2.0*mt.factor;
				PP_ref_db[i].phase_shearModulus  =  -3.0 * fa.phase_shearModulus + 3.0*q.phase_shearModulus + 2.0*mt.phase_shearModulus;

			}
			else{
				PP_ref_db[i] = G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												gv.PP_list[i], 
												state				);
			}

			sum_zel = 0;
			for (int j = 0; j < z_b.zEl_val; j++){
				
				/* If pure-phase contains an oxide absent in the bulk-rock then do not take it into account */
				if (PP_ref_db[i].Comp[z_b.zEl_array[j]] != 0.0){
					sum_zel += 1;
				}
				
				/* If pure-phase contains an oxide absent in the bulk-rock then do not take it into account */
				else{
					if (sum_zel != 0){
						gv.pp_flags[i][0] = 0;
						gv.pp_flags[i][1] = 0;
						gv.pp_flags[i][2] = 0;
						gv.pp_flags[i][3] = 1;
					}
					else{
						gv.pp_flags[i][0] = 1;
						gv.pp_flags[i][1] = 0;
						gv.pp_flags[i][2] = 1;
						gv.pp_flags[i][3] = 0;
					}
				}

			}

			/* If pure-phase contains an oxide absent in the bulk-rock then do not take it into account */
			if (gv.QFM_buffer == 0 && strcmp( gv.PP_list[i], "qfm") == 0 ){

				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}
			

			if (gv.verbose==1){
				printf(" %4s:  %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase);

				/* display molar composition */

				if (EM_database == 0){
					printf("\n S   A   C   M   F   K   N   T   O   M   H  \n");
				}
				else if (EM_database == 1){
					printf("\n S   A   C   M   F   K   N   T   O   H  \n");
				}
				else if (EM_database == 2 || EM_database == 6){
					printf("\n S   A   C   M   F   K   N   T   O   Cr  H  \n");
				}
				else if (EM_database == 4){
					printf("\n S   A   M   F   O   H   S\n");
				}
				else if (EM_database == 5){
					printf("\n S   A   C   M   F   K   N   O   H   C  \n");
				}

				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",PP_ref_db[i].Comp[j]);
				}
				printf("\n");

			}


		}
		if (gv.verbose==1){
			printf("\n");
		}		
		return gv;
};
