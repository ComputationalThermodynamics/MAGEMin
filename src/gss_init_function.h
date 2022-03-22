#ifndef __GSS_INIT_FUNCTION_H_
#define __GSS_INIT_FUNCTION_H_

#include "MAGEMin.h"
<<<<<<< HEAD

SS_ref G_SS_INIT_EM_function(SS_ref SS_ref_db, int EM_database, char *name, global_variable gv);
=======
>>>>>>> dev

SS_ref G_SS_INIT_EM_function(		SS_ref 			SS_ref_db, 
									int 			EM_database, 
									char 		   *name, 
									global_variable gv				);

csd_phase_set CP_INIT_function(		csd_phase_set 	cp, 
									global_variable gv				);
									
stb_system SP_INIT_function(		stb_system 		sp,
									global_variable gv				);

#endif
