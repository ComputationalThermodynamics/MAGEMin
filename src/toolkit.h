#ifndef __TOOLKIT_H_
#define __TOOLKIT_H_

#include "MAGEMin.h"

/* set of matrix operations */
void 	print_help(global_variable gv);
void 	_DCDCT_fct(int *id, double *result, double **A, int n_act_sf, int n_xeos);
void 	_DC_Null_fct(int *id, double *result, double **A, double *B, int n_xeos, int n_act_sf);
void 	_Epsilon_C_fct(int *id, double *result, double *A, double *b, int n_xeos, int n_sf);
void 	_Epsilon_J_fct(double *result, double *A, double *b, int n_xeos);
void 	_I_DC_Null_fct(int *id, double *result, double *A, double **B, double *eye, int n_act_sf, int n_xeos);
void 	_FillEyeMatrix(double *A, int n);
void 	get_act_sf_id(int *result, double *A, int n);
void 	inverseMatrix(double *A1, int n);
void 	MatMatMul( double **A, int nrowA, double **B, int ncolB, int common, double **C);
void 	VecMatMul(double *B1, double *A1, double *B, int n);
void 	MatVecMul(double *A1, double *br, double *n_vec, int n);
void 	pseudo_inverse(	double *matrix,
						double *B,
						int m,
						int n				);

int 	get_max_n_pc(int tot_pc, int n_pc);
int 	get_act_sf(double *A, int n);
int 	get_active_em(double *array, int n);
int 	EndsWithTail(char *name, char* tail);	
int    	RootBracketed(double x1,double x2);

double* norm_array(double *array, int size);
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

/*	Reduce row echelon form function (should be deleted eventually) */
typedef struct TMatrix {
	double **m;
	int nRows; int nCols;
} TMATRIX;

TMATRIX createMatrix(int nRows, int nCols);
TMATRIX rref(TMATRIX stoeMat, int *pivot, double tolerance);

void 	freeMatrix(TMATRIX oMatrix);
void 	cleanUpMatrix(TMATRIX stoeMat,  double tolerance);

/* function related to update on global variables structure */
global_variable get_pp_id(			global_variable  	 gv								);

global_variable get_ss_id(			global_variable  	 gv,
									csd_phase_set  		*cp								);

/* Melt-fraction correction for P-wave and S-wave velocities */		
void wave_melt_correction( 	 double  Kb_L,
							 double  Kb_S,
							 double  Ks_S,
							 double  rhoL,
							 double  rhoS,
							 double  Vp0,
							 double  Vs0,
							 double  meltFrac,
							 double  aspectRatio,
							 double *V_cor												);

/* This routine computes a correction of P-wave and S-wave velocities using melt fraction reduction.  */
double anelastic_correction( int 	 water,
							 double  Vs0,
							 double  P,
							 double  T 													);
																																	
global_variable get_sol_phase_infos( 	io_data 			 input_data,
										bulk_info 	 		 z_b,
										global_variable 	 gv,

										PP_ref  			*PP_ref_db,
										SS_ref  			*SS_ref_db,
										csd_phase_set  		*cp							);

#endif
