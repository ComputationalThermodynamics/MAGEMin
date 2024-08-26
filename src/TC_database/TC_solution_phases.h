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
#ifndef __TC_sol_H_
#define __TC_sol_H_

    #include "gss_init_function.h"
    #include "gss_function.h"
    #include "objective_functions.h"
    #include "NLopt_opt_function.h"

    /* include pseudocompounds */
    #include "SS_xeos_PC_mp.h" 				//mp is first, it contains the structure definition
    #include "SS_xeos_PC_mb.h" 
    #include "SS_xeos_PC_ig.h"
    #include "SS_xeos_PC_um.h"

#endif