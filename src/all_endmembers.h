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
#ifndef __TC_all_em_H_
#define __TC_all_em_H_

    /* This includes the activity model from THERMOCALC */
    #include "./TC_database/TC_endmembers.h"
    #include "./TC_database/TC_gem_function.h"
    /* This include Stixrude & Lithgow-Bertelloni solution phase models */
    #include "./SB_database/SB_endmembers.h"
    #include "./SB_database/SB_gem_function.h"
    //#include "STIX_solution_phases.h"

    /* This includes the Ghiorso/MELTS liquid model (research group "gh") */
    #include "./GH_database/GH_endmembers.h"
    #include "./GH_database/GH_PP_endmembers.h"
    #include "./GH_database/GH_gem_function.h"

#endif