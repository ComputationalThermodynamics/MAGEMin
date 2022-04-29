#ifndef __IO_FUNCTION_H_
#define __IO_FUNCTION_H_

#include "MAGEMin.h"

void read_in_data(		global_variable 				gv,
						io_data 					   *input_data,
						char    					   *file_name,
						int      			 			n_points	);

/**
	Output structure (send back data) 
*/ 
void AddResults_output_struct(	global_variable 		gv,
								struct bulk_info 		z_b,
								double P, double T,
								Databases 				DB,
								out_data 				output		);


out_data InitializeOutput(		global_variable 		gv,
								Databases 				DB			);

void FreeOutput(				out_data 				output);

#endif
