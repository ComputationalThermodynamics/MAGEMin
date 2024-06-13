/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __PGE_FUNCTION_H_
#define __PGE_FUNCTION_H_

#include "MAGEMin.h"
#include "all_solution_phases.h"

global_variable PGE(			bulk_info 			z_b,
								global_variable 	gv,
								
								PC_type 			*PC_read,
								obj_type 			*SS_objective,
								NLopt_type 			*NLopt_opt,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					);
								
global_variable init_LP(							bulk_info 	 		 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													PC_type 			*PC_read,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db,
													csd_phase_set  		*cp	
);

global_variable run_LP(								bulk_info 			 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db
);

global_variable run_LP_ig(							bulk_info 			 z_b,
													simplex_data 		*splx_data,
													global_variable 	 gv,
													
													PP_ref 				*PP_ref_db,
													SS_ref 				*SS_ref_db
);

global_variable LP(				bulk_info 			z_b,
								global_variable 	gv,
								PC_type 			*PC_read,
								
								obj_type 			*SS_objective,
								NLopt_type			*NLopt_opt,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					);

global_variable LP2(			bulk_info 			z_b,
								global_variable 	gv,
								
								obj_type 			*SS_objective,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					);

double norm_vector(double *array ,int n);

#endif
