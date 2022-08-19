#ifndef __DUMP_FUNCTION_H_
#define __DUMP_FUNCTION_H_

#include "MAGEMin.h"
void dump_init(global_variable gv);

void fill_output_struct(		global_variable 	gv,
								bulk_info 			z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
);


void dump_results_function(		global_variable 	gv,
								bulk_info 			z_b,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
);

void mergeParallelFiles(global_variable gv);

void mergeParallel_residual_Files(global_variable gv);

void mergeParallel_LocalMinima_Files(global_variable gv);

void mergeParallel_LevellingGamma_Files(global_variable gv);

#endif
