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
/**
  Function to calculate chemical potential of endmembers/pure phases  
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "MAGEMin.h"
#include "all_endmembers.h"
#include "gem_function.h"

PP_ref G_EM_function(	    char        *research_group,
                            int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
    PP_ref PP_ref_db;

    /**
        Here we select the right EOS for each research group 
    */
	if 	(strcmp(research_group, "tc") 	== 0){
        PP_ref_db = TC_G_EM_function(	EM_dataset,
                                        len_ox,
                                        id,
                                        bulk_rock,
                                        apo,
                                        P, T,
                                        name,
                                        state);
    }
    else{
        PP_ref_db = TC_G_EM_function(	EM_dataset,
                                        len_ox,
                                        id,
                                        bulk_rock,
                                        apo,
                                        P, T,
                                        name,
                                        state);
    }


    return PP_ref_db;
}