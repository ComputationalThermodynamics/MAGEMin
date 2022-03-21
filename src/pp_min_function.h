#ifndef __PP_MIN_FUNCTION_H_
#define __PP_MIN_FUNCTION_H_

#include "MAGEMin.h"

void pp_min_function(		global_variable gv,
							struct bulk_info z_b,
							PP_ref *PP_ref_db				);

/* initialize pure phase database using P-T conditions */
global_variable init_em_db(	int EM_database,
							struct bulk_info z_b,
							
							global_variable gv,
							PP_ref *PP_ref_db				);

#endif
