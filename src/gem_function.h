#ifndef __GEM_FUNCTION_H_
#define __GEM_FUNCTION_H_

double sum_array(double *array, int size);

int check_sign(double v1, double v2);

/*  Store pure phases composition, gbase, and bulk-rock factor */
typedef struct PP_refs {
	char   	Name[20];			/** Name */
    double 	Comp[11];        	/** composition [0-10]*/
    double 	gbase; 
    double 	gb_lvl;         	/**driving force, delta_G with G-hyperplane */
    double 	factor;
    double  phase_density;		/** molar density of the phase */
    double  phase_shearModulus;		/** molar density of the phase */
    double  phase_cp;			/** molar cp of the phase */
	double  phase_expansivity;
	double  phase_bulkModulus;
    double  volume;				/** molar volume of the phase */
    double  mass;				/** molar mass of the phase */
    
} PP_ref;

PP_ref G_EM_function(int EM_database, double *bulk_rock, double P, double T, char *name, char *state);

#endif
