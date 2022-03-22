#ifndef __PHASE_UPDATE_FUNCTION_H_
#define __PHASE_UPDATE_FUNCTION_H_

#include "MAGEMin.h"

global_variable phase_update_function(		struct bulk_info 	z_b,
											global_variable 	gv,
											
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp 
);

global_variable phase_merge_function(		struct bulk_info 	z_b,
											global_variable 	gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp 
);


struct str
{
	double value;
	int    index;
};

/* compare double function */
int cmp_dbl(const void *a, const void *b);
	
/* compare int function */
int cmp_int(const void *a, const void *b);

#endif
