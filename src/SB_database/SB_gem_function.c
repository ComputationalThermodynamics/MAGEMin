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
/**
  Function to calculate chemical potential of endmembers/pure phases for Stixrude database
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "SB_gem_function.h"

#define eps 1e-8


PP_ref SB_G_EM_function(	int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
	/* Get thermodynamic data */
	EM_db_sb EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_SB_EM_DB(p_id, EM_dataset);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[id[i]];
	}
	
	double kbar2bar = 1e3;
	double RTlnf 	= 0.0;
	double t0 		= 298.15;
	double p0 		= 0.001;
	double R  		= 0.0083144; 
 	
	/* fill structure to send back to main */
	PP_ref PP_ref_db;
	
	/* Calculate normalizing factor using bulk-rock composition */
	double factor  = 0.0;
	
	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (i = 0; i < len_ox; i++){
		fbc += bulk_rock[i]*apo[i];
	}
	
	/* Calculate the number of atom in the solution */
	double ape = 0.0;
	for (i = 0; i < len_ox; i++){
		ape += composition[i]*apo[i];
	}
	
	/* Calculate normalizing factor */
	factor = fbc/ape;

	strcpy(PP_ref_db.Name, name);
	for (i = 0; i < len_ox; i++){
		PP_ref_db.Comp[i] = composition[i];
	}
	// PP_ref_db.gbase   =  gbase;
	PP_ref_db.factor  =  factor;
	PP_ref_db.phase_shearModulus  =  (EM_return.input_4[0]*kbar2bar + (P - p0)*(EM_return.input_4[1])*kbar2bar + (T - t0)*(EM_return.input_4[2]))/kbar2bar;

	// printf(" %4s %+10f\n",name,PP_ref_db.gbase);
	// for (i = 0; i < len_ox; i++){
	// 	printf("%+10f",PP_ref_db.Comp[i]*PP_ref_db.factor); 
	// }
	// printf("\n");

	return (PP_ref_db);
}
