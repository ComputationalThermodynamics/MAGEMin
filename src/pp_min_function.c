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

			if 		(strcmp( gv.PP_list[i], "qif") 	== 0){

				PP_ref iron 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"iron", 
												state				);
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


				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = fa.Comp[j] -2.0*iron.Comp[j] - q.Comp[j];
				}		

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  fa.gbase -2.0*iron.gbase - q.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = fa.phase_shearModulus -2.0*iron.phase_shearModulus - q.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "nno") 	== 0){

				PP_ref Ni 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"Ni", 
												state				);
				PP_ref NiO 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"NiO", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = 2.*NiO.Comp[j] -2.*Ni.Comp[j];
				}		

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  2.*NiO.gbase - 2.*Ni.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 2.*NiO.phase_shearModulus - 2.*Ni.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "mw") 	== 0){

				PP_ref mt 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);
				PP_ref wu 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"wu", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = 2.0*mt.Comp[j] -6.0*wu.Comp[j];
				}		

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  2.0*mt.gbase -6.0*wu.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 2.0*mt.phase_shearModulus -6.0*wu.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "cco") 	== 0){

				PP_ref co2 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"CO2", 
												state				);
				PP_ref gph 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"gph", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = co2.Comp[j] - gph.Comp[j];
				}		

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  co2.gbase - gph.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = co2.phase_shearModulus - gph.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "hm") 	== 0){

				PP_ref mt 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);
				PP_ref hem 	= G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"hem", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = 6.0*hem.Comp[j] - 4.0*mt.Comp[j];
				}		

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  6.0*hem.gbase -4.0*mt.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 6.0*hem.phase_shearModulus -4.0*mt.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "qfm") 	== 0){

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

				/* Calculate the number of atoms in the bulk-rock composition */
				double fbc     = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					fbc += z_b.bulk_rock[j]*z_b.apo[j];
				}
				
				/* Calculate the number of atom in the solution */
				double ape = 0.0;
				for (int j = 0; j < gv.len_ox; j++){
					ape += PP_ref_db[i].Comp[j]*z_b.apo[j];
				}
				/* Calculate normalizing factor */
				double factor = fbc/ape;

				PP_ref_db[i].gbase   =  -3.0 * fa.gbase + 3.0*q.gbase + 2.0*mt.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  =  -3.0 * fa.phase_shearModulus + 3.0*q.phase_shearModulus + 2.0*mt.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "O2") 	== 0){
				PP_ref_db[i] = G_EM_function(	EM_database, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												0.001, 					//for computing oxygen fugacity pressure = 1bar, 0.001 kbar
												z_b.T, 
												gv.PP_list[i], 
												state				);
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

			}

			/* If pure-phase contains an oxide absent in the bulk-rock then do not take it into account */
			if (sum_zel != 0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}
			else{
				if (gv.pp_flags[i][0] != 0){ 			//here  we check if the pure phase is deactivated from the start (O2 for instance)
					gv.pp_flags[i][0] = 1;
					gv.pp_flags[i][1] = 0;
					gv.pp_flags[i][2] = 1;
					gv.pp_flags[i][3] = 0;
				}
			}

			/* If buffer not active then remove it */
			if (gv.pp_flags[i][4] == 1 && strcmp(gv.buffer,gv.PP_list[i]) != 0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}		

			if (gv.verbose==1){
				printf(" %4s:  %+10f %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase, PP_ref_db[i].factor);

				/* display molar composition */

				if (EM_database == 0){
					printf("\n S   A   C   M   F   K   N   T   O   Mn  H  \n");
				}
				else if (EM_database == 1){
					printf("\n S   A   C   M   F   K   N   T   O   H  \n");
				}
				else if (EM_database == 2 || EM_database == 3 || EM_database == 6 || EM_database == 5){
					printf("\n S   A   C   M   F   K   N   T   O   Cr  H  \n");
				}
				else if (EM_database == 4){
					printf("\n S   A   M   F   O   H   S\n");
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
