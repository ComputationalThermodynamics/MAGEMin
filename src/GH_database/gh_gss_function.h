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
#ifndef __GH_GSS_FUNCTION_H_
#define __GH_GSS_FUNCTION_H_

#include "../MAGEMin.h"

SS_ref G_SS_gh_EM_function(    global_variable  gv,
                                SS_ref           SS_ref_db,
                                int              EM_dataset,
                                bulk_info        z_b,
                                char            *name                  );

#endif
