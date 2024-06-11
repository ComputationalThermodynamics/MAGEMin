#ifndef __NLOPT_OPT_FUNCTION_H_
#define __NLOPT_OPT_FUNCTION_H_

#include "../MAGEMin.h"

typedef void (*sf_type) (		unsigned  		 m,
								double 			*result,
								unsigned 		 n,
								const double 	*x,
								double 			*grad,
								void 			*data		);

	
global_variable NLopt_global_opt_function(	bulk_info 			 z_b,
											global_variable 	 gv, 
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
	 										csd_phase_set  		*cp
);

SS_ref NLopt_opt_function(		global_variable 	gv, 
								SS_ref 				SS_ref_db,  
								int     			index				);

#endif

