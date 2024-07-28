/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus, Jamison Assunção
 **   Contributors : Dominguez, H., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __PP_MIN_FUNCTION_H_
#define __PP_MIN_FUNCTION_H_

#include "MAGEMin.h"

void pp_min_function(		global_variable  gv,
							bulk_info 		 z_b,
							PP_ref 			*PP_ref_db				);

/* initialize pure phase database using P-T conditions */
global_variable init_em_db(	int EM_database,
							bulk_info 		 z_b,
							
							global_variable  gv,
							PP_ref 			*PP_ref_db				);

#endif
