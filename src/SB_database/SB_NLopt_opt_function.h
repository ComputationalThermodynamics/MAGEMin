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
#ifndef __SB_NLOPT_OPT_FUNCTION_H_
#define __SB_NLOPT_OPT_FUNCTION_H_

#include "../MAGEMin.h"

typedef SS_ref (*NLopt_type) (		global_variable 	 gv,
									SS_ref 				 SS_ref_db		);

void SB_NLopt_opt_init(	        	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
#endif

