#ifndef __PGE_FUNCTION_H_
#define __PGE_FUNCTION_H_

#include "MAGEMin.h"

global_variable PGE(			bulk_info 			z_b,
								global_variable 	gv,
								
								obj_type 			*SS_objective,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					);

global_variable LP_PGE(			bulk_info 			z_b,
								global_variable 	gv,
								
								obj_type 			*SS_objective,
								simplex_data	    *splx_data,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp					);

double norm_vector(double *array ,int n);

#endif
