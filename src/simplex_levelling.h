#ifndef __run_levelling_function_H_
#define __run_levelling_function_H_

#include "MAGEMin.h"

void inverseMatrix(double *A1, int n);
void VecMatMul(double *B1, double *A1, double *B, int n);
void MatVecMul(double *A1, double *br, double *n_vec, int n);

global_variable Levelling(
	struct bulk_info z_b,
	global_variable gv,
	
	PP_ref *PP_ref_db,
	SS_ref *SS_ref_db,
	csd_phase_set  *cp
);

typedef struct simplex_datas
{
	/* global variables */
	double  *gamma_ps;		/** chemical potential of oxides (pure species round) 	*/
	double  *gamma_ss;		/** chemical potential of oxides (solution phase round) */
	double  *gamma_tot;		/** update global gamma									*/
	double  *gamma_delta;	/** delta gamma between two levelling rounds			*/
	
	double   min_F;		/** min F */
	int      ph2swp;	/** index of phase to swap */
	int      n_swp;     /** number of phase added to the reference assemblage */
	int      swp;       /** swap occured? */
	int     *pivot;		/** pivot point when doing RREF toget the rational basis of the null space */
	
	/* Reference assemblage */
	double  *A;			/** stoechiometry matrix */
	double  *A1;		/** inverse of stoechiometry matrix */
	int    **ph_id_A;	/** id of phases */
	
	double  *g0_A;		/** save reference gibbs energy of pseudocompound */
	double  *dG_A;		/** driving force matrix */
	double  *n_vec;		/** phase fractions */
	int    	 n_Ox;		/** number of active oxides */
	int		 len_ox;	/** number of total oxides */
	
	/* Potential candidates */
	int      n_pp;		/** number of pure phases */
	int      n_em_ss;	/** number of endmembers in solutions phases */
	
	double  *B;			/** stoechiometry matrix */
	double  *B1;		/** inverse of stoechiometry matrix entry to be added*/
	int     *ph_id_B;	/** id of phases */
	
	double   g0_B;		/** save reference gibbs energy of pseudocompound */
	double   dG_B;		/** driving force matrix */
	int      n_B;		/** number of pseudocompounds */
	
	int n_local_min;
	int n_filter;
	
} simplex_data;

void print_levelling(
	struct bulk_info z_b,
	
	global_variable gv,
	
	PP_ref *PP_ref_db,
	SS_ref *SS_ref_db
);

#endif
