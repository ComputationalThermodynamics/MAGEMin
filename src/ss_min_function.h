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
#ifndef __SS_MIN_FUNCTION_H_
#define __SS_MIN_FUNCTION_H_

#include "MAGEMin.h"
#include "all_solution_phases.h"

SS_ref SS_UPDATE_function(				global_variable 	 gv,
										SS_ref 				 SS_ref_db, 
										bulk_info 	 		 z_b,
										char    			*name			);		
								
csd_phase_set CP_UPDATE_function(		global_variable 	 gv,
										SS_ref 				 SS_ref_db,
										csd_phase_set  		 cp, 
										bulk_info 	 		 z_b			);		

global_variable split_cp(				global_variable 	 gv,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

void init_PGE_from_LP(					global_variable 	 gv,
										PC_type 			*PC_read,
										obj_type 			*SS_objective,
										bulk_info 	 		 z_b,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

void ss_min_PGE(						global_variable 	 gv,
										PC_type				*PC_read,
										obj_type 			*SS_objective,
										NLopt_type			*NLopt_opt,
										bulk_info 	 		 z_b,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp				);
										
void ss_min_LP(							global_variable 	 gv,
										PC_type				*PC_read,
										
										obj_type 			*SS_objective,
										NLopt_type			*NLopt_opt,
										bulk_info 	 		 z_b,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp				);

void copy_to_Ppc(						int 				 pc_check,
										int 				 add,
										int 				 ph_id,
										global_variable 	 gv,

										obj_type 			*SS_objective,
										SS_ref 			    *SS_ref_db		);

void copy_to_cp(						int 				 i, 
										int 				 ph_id,
										global_variable 	 gv,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

global_variable init_ss_db(				int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										SS_ref 				*SS_ref_db		);
#endif
