#ifndef __run_levelling_function_H_
#define __run_levelling_function_H_

#include "MAGEMin.h"

void inverseMatrix(					double 			   *A1, 
									int 				n			);

void VecMatMul(						double 			   *B1,
									double 			   *A1,
									double 			   *B,
									int 				n			);

void MatVecMul(						double 			   *A1,
									double 			   *br,
									double 			   *n_vec,
									int 				n			);

global_variable Levelling(			struct bulk_info 	z_b,
									global_variable 	gv,
									simplex_data	   *splx_data,
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db,
									csd_phase_set  	   *cp			);

void destroy_simplex_A(				simplex_data *splx_data			);

void destroy_simplex_B(				simplex_data *splx_data			);

void print_levelling(				struct bulk_info 	z_b,
									global_variable 	gv,
									
									PP_ref 			   *PP_ref_db,
									SS_ref 			   *SS_ref_db	);

#endif
