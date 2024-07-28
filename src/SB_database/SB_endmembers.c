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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 

#include "SB_endmembers.h"


/* Select required thermodynamic database */
EM_db_sb Access_SB_EM_DB(int id, int EM_dataset) {
	EM_db_sb Entry_EM;

	// if (EM_dataset == 2011){	
	//  	Entry_EM = arr_em_db_sb_2011[id]; 
	// }
	// else if (EM_dataset == 2024){		
	// 	Entry_EM = arr_em_db_sb_2024[id]; 
	// }

	// else{
	// 	printf(" Wrong endmember dataset, values should be 2011 or 2024\n");
	// 	printf(" -> using default 2011 to avoid ugly crash\n");
	// 	Entry_EM = arr_em_db_sb_2011[id]; 
	// }
	
	return Entry_EM;
}
