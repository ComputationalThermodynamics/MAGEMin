#ifndef __DUMP_FUNCTION_H_
#define __DUMP_FUNCTION_H_

#include "MAGEMin.h"
void dump_init(global_variable gv);

void dump_results_function(		global_variable 	gv,
								struct bulk_info 	z_b,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp
);

void mergeParallelFiles(global_variable gv);

void mergeParallel_LocalMinima_Files(global_variable gv);

void mergeParallel_LevellingGamma_Files(global_variable gv);

#endif
