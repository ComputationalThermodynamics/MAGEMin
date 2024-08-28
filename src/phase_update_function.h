/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __PHASE_UPDATE_FUNCTION_H_
#define __PHASE_UPDATE_FUNCTION_H_

#include "MAGEMin.h"

global_variable check_PC(					bulk_info 	 		 z_b,
											global_variable  	 gv,
											PC_type				*PC_read,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp				
);

global_variable check_PC_driving_force(		bulk_info 	 z_b,
											global_variable  	 gv,

											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp			
);

global_variable phase_update_function(		bulk_info 			z_b,
											global_variable 	gv,
											
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
											csd_phase_set  		*cp 
);

global_variable phase_merge_function(		bulk_info 			z_b,
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
