/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __NS_OPT_FUNCTION_H_
#define __NS_OPT_FUNCTION_H_

#include "../MAGEMin.h"

void TC_NS_init(					global_variable 	 gv,
									SS_ref 				*SS_ref_db		);
void SB_NS_init(					global_variable 	 gv,
									SS_ref 				*SS_ref_db		);

void NS_free(						SS_ref 				*SS_ref_db		);

void ns_reduce_system(				global_variable 	 gv,
									bulk_info 	 		 z_b,
									SS_ref 				*SS_ref_db		);

double ns_x_of_sf(					global_variable 	 gv,
									SS_ref 				*SS_ref_db,
									const double 		*sf				);

double ns_eval(						global_variable 	 gv,
									SS_ref 				*SS_ref_db,
									const double 		*sf,
									double 				*grad			);

SS_ref NS_opt_function(				global_variable 	 gv,
									SS_ref 				 SS_ref_db		);
int ns_project_x(					global_variable 	 gv,
									SS_ref 				*SS_ref_db,
									double 				*x				);
int ns_pc_mode(						global_variable 	 gv,
									SS_ref 				*SS_ref_db		);

#endif
