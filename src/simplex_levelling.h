/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __run_levelling_function_H_
#define __run_levelling_function_H_

#include "MAGEMin.h"
#include "all_solution_phases.h"

/* simplex levelling declaration part */
void update_dG(						simplex_data 	    *splx_data			);


void update_local_gamma(			double 				*A1, 
									double 				*g0_A, 
									double 				*gam, 
									int 				 n					);

void update_global_gamma( 			bulk_info 			 z_b,
									simplex_data 	    *splx_data			);

void update_global_gamma_LU( 		bulk_info 			 z_b,
									simplex_data 	    *splx_data			);


void swap_pure_phases(				bulk_info 	 		 z_b,
									simplex_data 		*splx_data,
									global_variable 	 gv,
									
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db			);

void swap_pure_endmembers(			bulk_info 	 		 z_b,
									simplex_data 		*splx_data,
									global_variable 	 gv,
										
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db			);


void swap_pseudocompounds(			bulk_info 	 		 z_b,
									simplex_data 		*splx_data,
									global_variable 	 gv,
										
									PP_ref 				*PP_ref_db,
									SS_ref 				*SS_ref_db			);

void swap_PGE_pseudocompounds(		bulk_info 	 		 z_b,
									simplex_data 	    *splx_data,
									global_variable 	 gv,
									
									PP_ref 			    *PP_ref_db,
									SS_ref 			    *SS_ref_db			);
		
global_variable Levelling(			bulk_info 			z_b,
									global_variable 	gv,

									PC_type            *PC_read,
									P2X_type		   *P2X_read,
									obj_type 		   *SS_objective,
									simplex_data	   *splx_data,
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db,
									csd_phase_set  	   *cp					);

global_variable Initial_guess(		bulk_info 			z_b,
									global_variable 	gv,

									PC_type            *PC_read,
									P2X_type 		   *P2X_read,
									simplex_data	   *splx_data,
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db,
									csd_phase_set  	   *cp					);

global_variable Metastable_calc(	bulk_info 			z_b,
									global_variable 	gv,

									PC_type            *PC_read,
									P2X_type 		   *P2X_read,
									simplex_data	   *splx_data,
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db,
									csd_phase_set  	   *cp					);	
															
void destroy_simplex_A(				simplex_data 	   *splx_data			);

void destroy_simplex_B(				simplex_data 	   *splx_data			);

void print_levelling(				bulk_info 			z_b,
									global_variable 	gv,
									
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db			);


#endif
