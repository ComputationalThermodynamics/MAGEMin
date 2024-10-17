#ifndef __TOOLKIT_H_
#define __TOOLKIT_H_

#include "MAGEMin.h"

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
/* set of matrix operations */
void 	print_help(global_variable gv);

bulk_info retrieve_bulk_PT(				global_variable      gv,
										io_data 		    *input_data,
										int					 sgleP,
										bulk_info 			 z_b			);

void convert_system_comp(				global_variable      gv,
										char 				*sys_in,
										bulk_info 			 z_b			);
										
void 	get_act_sf_id(int *result, double *A, int n);
void 	inverseMatrix(int *ipiv, double *A1, int n, double *work, int lwork);
void 	MatMatMul( double **A, int nrowA, double **B, int ncolB, int common, double **C);
void 	VecMatMul(double *B1, double *A1, double *B, int n);
void 	MatVecMul(double *A1, double *br, double *n_vec, int n);
void 	pseudo_inverse(	double *matrix,
						double *B,
						int m,
						int n				);

int 	get_act_sf(double *A, int n);
int 	get_active_em(double *array, int n);
int 	EndsWithTail(char *name, char* tail);	
int    	RootBracketed(double x1,double x2);

void 	print_1D_double_array(double nx, double *array, char *title);
void 	print_2D_double_array(double nx, double ny, double **array, char *title);
void 	print_1D_int_array(double nx, int *array, char *title);

double  rnd(double a);
double  SUPCRT_to_HSC(double *ElH, double *comp, int size);
double  HSC_to_SUPCRT(double *ElH, double *comp, int size);
double* norm_array(double *array, int size);
double  sum_norm_xipi(double *xi, double *pi, int size);
double  sign(double x);
double  AFunction(double x, double *data);
double  Minimum(double x1,double x2);
double  Maximum(double x1,double x2);
double 	norm_vector(double *array ,int n);
double 	euclidean_distance(double *array1 ,double *array2 ,int n);
double 	partial_euclidean_distance(double *array1 ,double *array2 ,int n);
double 	VecVecMul(double *B0, double *B1, int n);
double 	BrentRoots(  double 	x1, 
					double 	x2,
					double *data, 
					double 	Tolerance,
					
					int 	mode,
					int		maxIterations,
					double *valueAtRoot,
					
					int    *niter, 
					int    *error 			);

/* printing function (verbose) */
void print_cp(					global_variable gv,
								csd_phase_set  *cp
															);

void print_SS_informations(		global_variable gv,
								SS_ref SS_ref_db,
								int		iss					);
								
/* functon related to hyperplane manipulation */
SS_ref rotate_hyperplane(		global_variable gv,
								SS_ref SS_ref_db			);
								
SS_ref non_rot_hyperplane(	global_variable gv,
							SS_ref 			SS_ref_db		);
							
SS_ref raw_hyperplane(		global_variable  gv,
							SS_ref 			 SS_ref_db,
							double 			*gb				);


SS_ref restrict_SS_HyperVolume(	global_variable gv, 
								SS_ref SS_ref_db,
								double box_size				);		

SS_ref check_SS_bounds(			global_variable gv, 
								SS_ref SS_ref_db			);																	

/* function related to update on global variables structure */
global_variable get_pp_id(			global_variable  	 gv								);

global_variable get_ss_id(			global_variable  	 gv,
									csd_phase_set  		*cp								);

/* Melt-fraction correction for P-wave and S-wave velocities */		
global_variable wave_melt_correction( 	global_variable     gv,
										bulk_info 			z_b,	
										double  			aspectRatio					);
										
/* This routine computes a correction of P-wave and S-wave velocities using melt fraction reduction.  */
double anelastic_correction( int 	 water,
							 double  Vs0,
							 double  P,
							 double  T 												);
																																	
global_variable compute_phase_mol_fraction(	global_variable 	 gv,
											bulk_info 	 		 z_b,
											PP_ref  			*PP_ref_db,
											SS_ref  			*SS_ref_db,
											csd_phase_set  		*cp					);

global_variable compute_activities(			int					 EM_database,	
											global_variable 	 gv,
											PP_ref  			*PP_ref_db,
											bulk_info 			 z_b				);

global_variable compute_density_volume_modulus(		int 				 EM_database,
													bulk_info 	 		 z_b,
													global_variable 	 gv,
													PP_ref  			*PP_ref_db,
													SS_ref  			*SS_ref_db,
													csd_phase_set  		*cp					);									
#endif
