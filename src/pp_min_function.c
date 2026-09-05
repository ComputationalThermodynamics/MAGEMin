/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
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
#include "GH_database/GH_gem_function.h"


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
		double buffer_n;
		char state[] = "equilibrium";	
		int sum_zel;
		for (int i = 0; i < gv.len_pp; i++){
			if 		(strcmp( gv.PP_list[i], "qif") 	== 0){

				PP_ref iron 	= G_EM_function(gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"iron", 
												state				);
				PP_ref q 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"q", 
												state				);
				PP_ref fa 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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

				PP_ref Ni 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"Ni", 
												state				);
				PP_ref NiO 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = 2.0*NiO.Comp[j] -2.0*Ni.Comp[j];
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

				PP_ref_db[i].gbase   =  2.0*NiO.gbase - 2.0*Ni.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 2.0*NiO.phase_shearModulus - 2.0*Ni.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
 			else if (strcmp( gv.PP_list[i], "aH2O") == 0){

				PP_ref H2O 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"H2O", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = H2O.Comp[j];
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
				if (gv.buffer_n < 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + H2O.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = H2O.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
 			else if (strcmp( gv.PP_list[i], "aO2") == 0){

				PP_ref O2 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"O2", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = O2.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + O2.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = O2.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}			
 			else if (strcmp( gv.PP_list[i], "aMgO") == 0){


				PP_ref MgO 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"per", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = MgO.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + MgO.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = MgO.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
 			else if (strcmp( gv.PP_list[i], "aFeO") == 0){

				PP_ref FeO 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"fper", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = FeO.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + FeO.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = FeO.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "aAl2O3") == 0){

				PP_ref Al2O3 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"cor", 
													state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = Al2O3.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + Al2O3.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = Al2O3.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "aTiO2") == 0){

				PP_ref TiO2 	= G_EM_function(	gv.research_group,
                                                	gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"ru", 
													state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = TiO2.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}

				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + TiO2.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = TiO2.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "mw") 	== 0){

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);
				PP_ref wu 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
			else if (strcmp( gv.PP_list[i], "iw") 	== 0){

				PP_ref iron 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"iron", 
													state				);
													
				PP_ref wu 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = 2.0*wu.Comp[j] -2.0*iron.Comp[j];
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

				PP_ref_db[i].gbase   =  2.0*wu.gbase -2.0*iron.gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 2.0*wu.phase_shearModulus -2.0*iron.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "cco") 	== 0){

				PP_ref co2 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"CO2", 
												state				);
				PP_ref gph 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);
				PP_ref hem 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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

				PP_ref q 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"q", 
												state				);
				PP_ref fa 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"fa", 
												state				);

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
				PP_ref_db[i] = G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												0.001, 					//for computing oxygen fugacity pressure = 1bar, 0.001 kbar
												z_b.T, 
												gv.PP_list[i], 
												state				);
			}
			else if (i >= gv.len_pp - gv.n_mu_fix){
				/* native mu-fix buffer phase (buffer= style, single positive fictive phase): see notes/mu-mu-native-legendre-transform.md */
				int k = i - (gv.len_pp - gv.n_mu_fix);
				int valid_idx = (gv.mu_fix_idx[k] >= 0 && gv.mu_fix_idx[k] < gv.len_ox);

				if (!valid_idx){
					/* n_mu_fix/mu_fix_idx are inconsistent (e.g. fewer indices
					   provided than n_mu_fix, leaving this slot at its -1
					   default) - fall back to an inert (all-zero composition,
					   zero gbase) phase instead of reading z_b.apo[] out of
					   bounds below. Zero Comp + zero gbase gives an exactly-
					   zero driving force, so this can never be favorably
					   swapped into the basis regardless of pp_flags. */
					printf(" WARNING: mu_fix_idx[%d]=%d is out of range [0,%d) - disabling fictive phase '%s'\n", k, gv.mu_fix_idx[k], gv.len_ox, gv.PP_list[i]);
				}

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = (valid_idx && j == gv.mu_fix_idx[k]) ? 1.0 : 0.0;
				}

				double factor = 1.0;
				if (valid_idx){
					double fbc = 0.0;
					for (int j = 0; j < gv.len_ox; j++){
						fbc += z_b.bulk_rock[j]*z_b.apo[j];
					}
					factor = fbc/z_b.apo[gv.mu_fix_idx[k]];
				}

				PP_ref_db[i].gbase              = valid_idx ? gv.mu_fix_val[k] : 0.0;
				PP_ref_db[i].factor             = factor;
				PP_ref_db[i].factor_norm        = 1.0;
				PP_ref_db[i].phase_density      = 0.0;
				PP_ref_db[i].phase_shearModulus = 0.0;
				PP_ref_db[i].phase_bulkModulus  = 0.0;
				PP_ref_db[i].phase_cp           = 0.0;
				PP_ref_db[i].volume             = 0.0;
				PP_ref_db[i].mass               = 0.0;
			}
			else{
				PP_ref_db[i] = G_EM_function(	gv.research_group,
                                                gv.EM_dataset,
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
			if ( (gv.pp_flags[i][4] == 1 && strcmp(gv.buffer,gv.PP_list[i]) != 0) || (gv.act_PP[i] == 0) ){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}

			if (gv.verbose==1){
				printf("\n %4s:  %+10f %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase, PP_ref_db[i].factor);

				/* display molar composition */
				if (EM_database == 0){
					printf(" S   A   C   M   F   K   N   T   O   Mn  H\n");
				}
				else if (EM_database == 1){
					printf(" S   A   C   M   F   K   N   T   O   H\n");
				}
				else if (EM_database == 11){
					printf(" S   A   C   M   F   K   N   T   O   H\n");
				}
				else if (EM_database == 2){
					printf(" S   A   C   M   F   K   N   T   O   Cr  H\n");
				}
				else if (EM_database == 22){
					printf(" S   A   C   M   F   K   N   T   O   Cr\n");
				}
				else if (EM_database == 3){
					printf(" S   A   C   M   F   K   N   T   O   Cr\n");
				}
				else if (EM_database == 4){
					printf(" S   A   M   F   O   H   S\n");
				}
				else if (EM_database == 5){
					printf(" S   A   M   F   O   H   S   C   N   Cr  CO\n");
				}
				else if (EM_database == 6){
					printf(" S   A   M   F   O   H   S   C   N \n");
				}
				else if (EM_database == 7){
					printf(" S   A   C   M   F   K   N   T   O   Mn  H   CO2 S\n");
				}
				else if (EM_database == 8){
					printf(" S   A   C   M   F   K   N   T   O   Mn  Cr  H   CO2 S\n");
				}
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",PP_ref_db[i].Comp[j]);
				}
				printf("\n");

				if (gv.pp_flags[i][3] == 1){
					printf(" Pure phase is not considered in the calculation\n");
				}
			}

		}
		if (gv.verbose==1){
			printf("\n");
		}
		return gv;
};



/** 
  initialize pure phase database */
global_variable init_em_db_sb(	int 				EM_database,
								bulk_info 			z_b,
								global_variable 	gv,
								PP_ref 			   *PP_ref_db
){

		/* initialize endmember database */
		double buffer_n;
		char state[] = "equilibrium";	
		int sum_zel;
		for (int i = 0; i < gv.len_pp; i++){
		
 			if (strcmp( gv.PP_list[i], "aMgO") == 0){

				PP_ref MgO 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"pe", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = MgO.Comp[j];
					if (gv.EM_database == 1 || gv.EM_database == 2){
						PP_ref_db[i].Comp[j] /= 4.0; // divide by 4 for sb21
					}
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				double G0 = MgO.gbase;
				if (gv.EM_database == 1 || gv.EM_database == 2){
					G0 /= 4.0; // divide by 4 for sb21
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + G0;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = MgO.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
 			else if (strcmp( gv.PP_list[i], "aFeO") == 0){

				PP_ref FeO 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = FeO.Comp[j];
					if (gv.EM_database == 1 || gv.EM_database == 2 ){
						PP_ref_db[i].Comp[j] /= 4.0; // divide by 4 for sb21
					}
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				double G0 = FeO.gbase;
				if (gv.EM_database == 1 || gv.EM_database == 2){
					G0 /= 4.0; // divide by 4 for sb21
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + G0;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = FeO.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "aAl2O3") == 0){

				PP_ref Al2O3 	= G_EM_function(	gv.research_group,
                                               		gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"co", 
													state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = Al2O3.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}

				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + Al2O3.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = Al2O3.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
 			else if (strcmp( gv.PP_list[i], "aO2") == 0){

				PP_ref O2 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"O2", 
												state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = O2.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + O2.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = O2.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "qif") 	== 0){

				PP_ref fea 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"fea", 
													state				);
				PP_ref fee 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"fee", 
													state				);
				PP_ref feg 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"feg", 
													state				);		
													
				PP_ref q 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"qtz", 
												state				);
				PP_ref coe	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"coe", 
												state				);
				PP_ref st 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"st", 
												state				);
				PP_ref fa 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = fa.Comp[j] -2.0*fea.Comp[j] - q.Comp[j];
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
				double min_Fe_gbase = fea.gbase;
				if (fee.gbase < min_Fe_gbase) { min_Fe_gbase = fee.gbase; }
				if (feg.gbase < min_Fe_gbase) { min_Fe_gbase = feg.gbase; }

				double min_SiO2_gbase = q.gbase;
				if (coe.gbase < min_SiO2_gbase) { min_SiO2_gbase = coe.gbase; }

				if (st.gbase  < min_SiO2_gbase) { min_SiO2_gbase = st.gbase;  }
				PP_ref_db[i].gbase   =  fa.gbase -2.0*min_Fe_gbase - min_SiO2_gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = fa.phase_shearModulus -2.0*fea.phase_shearModulus - q.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "mw") 	== 0){


				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mag", 
												state				);

				PP_ref smt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"smag", 
												state				);

				PP_ref hmt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"hmag", 
												state				);
												
				PP_ref wu 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = 2.0*mt.Comp[j] -3.0*wu.Comp[j]/2.0;
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

				double min_mag_gbase = mt.gbase;
				if (smt.gbase < min_mag_gbase) { min_mag_gbase = smt.gbase; }
				if (hmt.gbase < min_mag_gbase) { min_mag_gbase = hmt.gbase; }

				PP_ref_db[i].gbase   =  2.0*min_mag_gbase -3.0*wu.gbase/2.0 + z_b.T*0.019145*gv.buffer_n; // divided by 4 Fe4O4 -> FeO
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 2.0*mt.phase_shearModulus -3.0*wu.phase_shearModulus/2.0;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "iw") 	== 0){

				PP_ref fea 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"fea", 
													state				);
				PP_ref fee 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"fee", 
													state				);
				PP_ref feg 	= G_EM_function(	gv.research_group,
													gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"feg", 
													state				);		
													
													
				PP_ref wu 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
					PP_ref_db[i].Comp[j] = wu.Comp[j]/2.0 -2.0*fea.Comp[j];
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

				double min_Fe_gbase = fea.gbase;
				if (fee.gbase < min_Fe_gbase) { min_Fe_gbase = fee.gbase; }
				if (feg.gbase < min_Fe_gbase) { min_Fe_gbase = feg.gbase; }

				PP_ref_db[i].gbase   =  wu.gbase/2.0 -2.0*min_Fe_gbase+ z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = wu.phase_shearModulus/2.0 -2.0*fea.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "hm") 	== 0){

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mag", 
												state				);

				PP_ref smt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"smag", 
												state				);

				PP_ref hmt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"hmag", 
												state				);

				PP_ref hem 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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

				double min_mag_gbase = mt.gbase;
				if (smt.gbase < min_mag_gbase) { min_mag_gbase = smt.gbase; }
				if (hmt.gbase < min_mag_gbase) { min_mag_gbase = hmt.gbase; }

				PP_ref_db[i].gbase   =  6.0*hem.gbase -4.0*min_mag_gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = 6.0*hem.phase_shearModulus -4.0*mt.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "qfm") 	== 0){

				PP_ref q 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"qtz", 
												state				);
				PP_ref coe	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"coe", 
												state				);
				PP_ref st 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"st", 
												state				);

				PP_ref fa 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"fa", 
												state				);

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mag", 
												state				);

				PP_ref smt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"smag", 
												state				);

				PP_ref hmt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"hmag", 
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

				double min_SiO2_gbase = q.gbase;
				if (coe.gbase < min_SiO2_gbase) { min_SiO2_gbase = coe.gbase; }
				if (st.gbase  < min_SiO2_gbase) { min_SiO2_gbase = st.gbase;  }

				double min_mag_gbase = mt.gbase;
				if (smt.gbase < min_mag_gbase) { min_mag_gbase = smt.gbase; }
				if (hmt.gbase < min_mag_gbase) { min_mag_gbase = hmt.gbase; }

				PP_ref_db[i].gbase   =  -3.0 * fa.gbase + 3.0*min_SiO2_gbase + 2.0*min_mag_gbase + z_b.T*0.019145*gv.buffer_n;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  =  -3.0 * fa.phase_shearModulus + 3.0*q.phase_shearModulus + 2.0*mt.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}
			else if (strcmp( gv.PP_list[i], "O2") 	== 0){
				PP_ref_db[i] = G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												0.001, 					//for computing oxygen fugacity pressure = 1bar, 0.001 kbar
												z_b.T, 
												gv.PP_list[i], 
												state				);
			}
			else if (i >= gv.len_pp - gv.n_mu_fix){
				/* native mu-fix buffer phase (buffer= style, single positive fictive phase): see notes/mu-mu-native-legendre-transform.md */
				int k = i - (gv.len_pp - gv.n_mu_fix);
				int valid_idx = (gv.mu_fix_idx[k] >= 0 && gv.mu_fix_idx[k] < gv.len_ox);

				if (!valid_idx){
					/* n_mu_fix/mu_fix_idx are inconsistent (e.g. fewer indices
					   provided than n_mu_fix, leaving this slot at its -1
					   default) - fall back to an inert (all-zero composition,
					   zero gbase) phase instead of reading z_b.apo[] out of
					   bounds below. Zero Comp + zero gbase gives an exactly-
					   zero driving force, so this can never be favorably
					   swapped into the basis regardless of pp_flags. */
					printf(" WARNING: mu_fix_idx[%d]=%d is out of range [0,%d) - disabling fictive phase '%s'\n", k, gv.mu_fix_idx[k], gv.len_ox, gv.PP_list[i]);
				}

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = (valid_idx && j == gv.mu_fix_idx[k]) ? 1.0 : 0.0;
				}

				double factor = 1.0;
				if (valid_idx){
					double fbc = 0.0;
					for (int j = 0; j < gv.len_ox; j++){
						fbc += z_b.bulk_rock[j]*z_b.apo[j];
					}
					factor = fbc/z_b.apo[gv.mu_fix_idx[k]];
				}

				PP_ref_db[i].gbase              = valid_idx ? gv.mu_fix_val[k] : 0.0;
				PP_ref_db[i].factor             = factor;
				PP_ref_db[i].factor_norm        = 1.0;
				PP_ref_db[i].phase_density      = 0.0;
				PP_ref_db[i].phase_shearModulus = 0.0;
				PP_ref_db[i].phase_bulkModulus  = 0.0;
				PP_ref_db[i].phase_cp           = 0.0;
				PP_ref_db[i].volume             = 0.0;
				PP_ref_db[i].mass               = 0.0;
			}
			else{
				PP_ref_db[i] = G_EM_function(	gv.research_group,
                                                gv.EM_dataset,
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
			if ( (gv.pp_flags[i][4] == 1 && strcmp(gv.buffer,gv.PP_list[i]) != 0) ){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}

			if (gv.verbose==1){
				printf("\n %4s:  %+10f %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase, PP_ref_db[i].factor);

				/* display molar composition */

				if (EM_database == 0){
					printf(" S   C   A   F   M   N\n");
				}
				else if (EM_database == 1){
					printf(" S   C   A   F   M   N\n");
				}
				else if (EM_database == 2){
					printf(" S   C   A   M   N   O   Cr  Fe\n");
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

/**
    Initialize pure-phase database for the "gh" (Ghiorso/MELTS) research
    group: a handful of common rock-forming pure phases ported from xMELTS'
    own solid-phase database (quartz, cristobalite, tridymite, corundum,
    sillimanite, rutile, sphene) plus O2 and H2O, all resolved directly by
    GH_G_EM_function (see GH_PP_endmembers.h). O2 is evaluated at a fixed
    1 bar reference pressure regardless of system pressure, mirroring
    init_em_db/init_em_db_sb's convention for the same fO2-buffer phase.
*/
global_variable init_em_db_gh(	int 				EM_database,
								bulk_info 			z_b,
								global_variable 	gv,
								PP_ref 			   *PP_ref_db
){
		/* runs before GH_SS_objective_init_function (MAGEMin.c's own call
		   order), so GH_actual_EM_database must ALSO be set here - not
		   just there - otherwise the very first gbase computation for gh
		   (this function, and init_ss_db_gh's own equivalent) would still
		   see the stale default. See GH_gem_function.c's header comment
		   and [[gh-multicalibration-xmelts-rmelts-pmelts]]. */
		double buffer_n;

		GH_actual_EM_database = gv.EM_database;
		char state[] = "equilibrium";
		int sum_zel;
		for (int i = 0; i < gv.len_pp; i++){
			if (strcmp( gv.PP_list[i], "O2") == 0){
				PP_ref_db[i] = G_EM_function(	gv.research_group,
												gv.EM_dataset,
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock,
												z_b.apo,
												0.001,				//for computing oxygen fugacity pressure = 1bar, 0.001 kbar
												z_b.T,
												gv.PP_list[i],
												state					);
			}
			else if (strcmp( gv.PP_list[i], "aAl2O3") == 0){

				PP_ref Al2O3 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"cor", 
													state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = Al2O3.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}
				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + Al2O3.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = Al2O3.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "aTiO2") == 0){

				PP_ref TiO2 	= G_EM_function(	gv.research_group,
                                                	gv.EM_dataset, 
													gv.len_ox,
													z_b.id,
													z_b.bulk_rock, 
													z_b.apo, 
													z_b.P, 
													z_b.T, 
													"ru", 
													state				);

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = TiO2.Comp[j];
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
				if (gv.buffer_n <= 1e-8){
					buffer_n = 1e-8;
				}
				else if (gv.buffer_n >= 1.0){
					buffer_n = 1.0-1e-8;
				}
				else{
					buffer_n = gv.buffer_n;
				}

				PP_ref_db[i].gbase   =  z_b.R * z_b.T*log(buffer_n) + TiO2.gbase;
				PP_ref_db[i].factor  =  factor;
				PP_ref_db[i].phase_shearModulus  = TiO2.phase_shearModulus;
				gv.pp_flags[i][4] 	= 1;
			}	
			else if (strcmp( gv.PP_list[i], "qfm") 	== 0){

				PP_ref q 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"q", 
												state				);
				PP_ref fa 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"fa", 
												state				);

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
			else if (strcmp( gv.PP_list[i], "hm") 	== 0){

				PP_ref mt 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock, 
												z_b.apo, 
												z_b.P, 
												z_b.T, 
												"mt", 
												state				);
				PP_ref hem 	= G_EM_function(	gv.research_group,
                                                gv.EM_dataset, 
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
			else if (i >= gv.len_pp - gv.n_mu_fix){
				/* native mu-fix buffer phase (buffer= style, single positive fictive phase): see notes/mu-mu-native-legendre-transform.md */
				int k = i - (gv.len_pp - gv.n_mu_fix);
				int valid_idx = (gv.mu_fix_idx[k] >= 0 && gv.mu_fix_idx[k] < gv.len_ox);

				if (!valid_idx){
					/* n_mu_fix/mu_fix_idx are inconsistent (e.g. fewer indices
					   provided than n_mu_fix, leaving this slot at its -1
					   default) - fall back to an inert (all-zero composition,
					   zero gbase) phase instead of reading z_b.apo[] out of
					   bounds below. Zero Comp + zero gbase gives an exactly-
					   zero driving force, so this can never be favorably
					   swapped into the basis regardless of pp_flags. */
					printf(" WARNING: mu_fix_idx[%d]=%d is out of range [0,%d) - disabling fictive phase '%s'\n", k, gv.mu_fix_idx[k], gv.len_ox, gv.PP_list[i]);
				}

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = (valid_idx && j == gv.mu_fix_idx[k]) ? 1.0 : 0.0;
				}

				double factor = 1.0;
				if (valid_idx){
					double fbc = 0.0;
					for (int j = 0; j < gv.len_ox; j++){
						fbc += z_b.bulk_rock[j]*z_b.apo[j];
					}
					factor = fbc/z_b.apo[gv.mu_fix_idx[k]];
				}

				PP_ref_db[i].gbase              = valid_idx ? gv.mu_fix_val[k] : 0.0;
				PP_ref_db[i].factor             = factor;
				PP_ref_db[i].factor_norm        = 1.0;
				PP_ref_db[i].phase_density      = 0.0;
				PP_ref_db[i].phase_shearModulus = 0.0;
				PP_ref_db[i].phase_bulkModulus  = 0.0;
				PP_ref_db[i].phase_cp           = 0.0;
				PP_ref_db[i].volume             = 0.0;
				PP_ref_db[i].mass               = 0.0;
			}
			else{
				PP_ref_db[i] = G_EM_function(	gv.research_group,
												gv.EM_dataset,
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock,
												z_b.apo,
												z_b.P,
												z_b.T,
												gv.PP_list[i],
												state					);
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
			if ( (gv.pp_flags[i][4] == 1 && strcmp(gv.buffer,gv.PP_list[i]) != 0) || (gv.act_PP[i] == 0) ){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}

			if (gv.verbose==1){
				printf("\n %4s:  %+10f %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase, PP_ref_db[i].factor);
				if (GH_actual_EM_database == 2){
					printf(" S   A   C   M   F   K   N   T   O   Mn  Cr  H\n");
				}
				else {
					printf(" S   A   C   M   F   K   N   T   O   Mn  Cr  H   CO2\n");
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

/* Source: Pourteau et al. (2014), Contrib Mineral Petrol; Berman (1988)
   EOS formalism via the THERIAK JUN92.bs lineage. */
global_variable init_em_db_br(	int 				EM_database,
								bulk_info 			z_b,
								global_variable 	gv,
								PP_ref 			   *PP_ref_db
){
		char state[] = "equilibrium";
		int sum_zel;
		for (int i = 0; i < gv.len_pp; i++){
			if (i >= gv.len_pp - gv.n_mu_fix){
				/* native mu-fix buffer phase (buffer= style, single positive fictive phase): see notes/mu-mu-native-legendre-transform.md */
				int k = i - (gv.len_pp - gv.n_mu_fix);
				int valid_idx = (gv.mu_fix_idx[k] >= 0 && gv.mu_fix_idx[k] < gv.len_ox);

				if (!valid_idx){
					/* n_mu_fix/mu_fix_idx are inconsistent (e.g. fewer indices
					   provided than n_mu_fix, leaving this slot at its -1
					   default) - fall back to an inert (all-zero composition,
					   zero gbase) phase instead of reading z_b.apo[] out of
					   bounds below. Zero Comp + zero gbase gives an exactly-
					   zero driving force, so this can never be favorably
					   swapped into the basis regardless of pp_flags. */
					printf(" WARNING: mu_fix_idx[%d]=%d is out of range [0,%d) - disabling fictive phase '%s'\n", k, gv.mu_fix_idx[k], gv.len_ox, gv.PP_list[i]);
				}

				strcpy(PP_ref_db[i].Name, gv.PP_list[i]);
				for (int j = 0; j < gv.len_ox; j++){
					PP_ref_db[i].Comp[j] = (valid_idx && j == gv.mu_fix_idx[k]) ? 1.0 : 0.0;
				}

				double factor = 1.0;
				if (valid_idx){
					double fbc = 0.0;
					for (int j = 0; j < gv.len_ox; j++){
						fbc += z_b.bulk_rock[j]*z_b.apo[j];
					}
					factor = fbc/z_b.apo[gv.mu_fix_idx[k]];
				}

				PP_ref_db[i].gbase              = valid_idx ? gv.mu_fix_val[k] : 0.0;
				PP_ref_db[i].factor             = factor;
				PP_ref_db[i].factor_norm        = 1.0;
				PP_ref_db[i].phase_density      = 0.0;
				PP_ref_db[i].phase_shearModulus = 0.0;
				PP_ref_db[i].phase_bulkModulus  = 0.0;
				PP_ref_db[i].phase_cp           = 0.0;
				PP_ref_db[i].volume             = 0.0;
				PP_ref_db[i].mass               = 0.0;
			}
			else{
				PP_ref_db[i] = G_EM_function(	gv.research_group,
												gv.EM_dataset,
												gv.len_ox,
												z_b.id,
												z_b.bulk_rock,
												z_b.apo,
												z_b.P,
												z_b.T,
												gv.PP_list[i],
												state					);
			}

			sum_zel = 0;
			for (int j = 0; j < z_b.zEl_val; j++){
				if (PP_ref_db[i].Comp[z_b.zEl_array[j]] != 0.0){
					sum_zel += 1;
				}
			}

			if (sum_zel != 0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}
			else{
				if (gv.pp_flags[i][0] != 0){
					gv.pp_flags[i][0] = 1;
					gv.pp_flags[i][1] = 0;
					gv.pp_flags[i][2] = 1;
					gv.pp_flags[i][3] = 0;
				}
			}

			if (gv.act_PP[i] == 0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}

			if (gv.verbose==1){
				printf("\n %4s:  %+10f %+10f\n",gv.PP_list[i],PP_ref_db[i].gbase, PP_ref_db[i].factor);
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
