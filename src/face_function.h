#ifndef __FACE_FUNCTION_H_
#define __FACE_FUNCTION_H_

#include "MAGEMin.h"

global_variable face_detect(			bulk_info 	 		 z_b,
										global_variable 	 gv,
										PP_ref  			*PP_ref_db,
										SS_ref  			*SS_ref_db		);

double face_deviation(					global_variable 	 gv,
										const double 		*comp			);

void face_project_gamma(				global_variable 	 gv,
										double 				*gam			);

#endif
