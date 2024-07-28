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
#ifndef __DUMP_FUNCTION_H_
#define __DUMP_FUNCTION_H_

#include "MAGEMin.h"
void dump_init(global_variable gv);

void fill_output_struct(		global_variable 	 gv,
								simplex_data	    *splx_data,
								bulk_info 			 z_b,

								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
);


void save_results_function(		global_variable 	gv,
								bulk_info 			z_b,
								PP_ref 				*PP_ref_db,
								SS_ref 				*SS_ref_db,
								csd_phase_set  		*cp,
								stb_system  		*sp
);

void mergeParallelFiles(global_variable gv);

void mergeParallel_matlab(global_variable gv);

void mergeParallel_residual_Files(global_variable gv);

void mergeParallel_LocalMinima_Files(global_variable gv);

void mergeParallel_LevellingGamma_Files(global_variable gv);

#endif
