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
#ifndef __GH_GEM_FUNCTION_H_
#define __GH_GEM_FUNCTION_H_

#include "../gem_function.h"


extern int GH_actual_EM_database;
extern int GH_H2O_liquid_context;

PP_ref GH_G_EM_function(   int          EM_database,
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock,
                            double      *apo,
                            double       P,
                            double       T,
                            char        *name,
                            char        *state          );

#endif
