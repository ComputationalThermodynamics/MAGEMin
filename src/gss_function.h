#ifndef __GSS_FUNCTION_H_
#define __GSS_FUNCTION_H_

#include "MAGEMin.h"
SS_ref G_SS_mp_EM_function(	global_variable  gv, 
							SS_ref 			 SS_ref_db,
							int 			 EM_database,
							bulk_info 		 z_b,
							char 			*name					);

SS_ref G_SS_mb_EM_function(	global_variable  gv, 
							SS_ref 			 SS_ref_db,
							int 			 EM_database,
							bulk_info 		 z_b,
							char 			*name					);

SS_ref G_SS_ig_EM_function(	global_variable  gv, 
							SS_ref 			 SS_ref_db,
							int 			 EM_database,
							bulk_info 		 z_b,
							char 			*name					);

SS_ref G_SS_um_EM_function(	global_variable  gv, 
							SS_ref 			 SS_ref_db,
							int 			 EM_database,
							bulk_info 		 z_b,
							char 			*name					);
#endif
