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