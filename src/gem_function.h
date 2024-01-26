#ifndef __GEM_FUNCTION_H_
#define __GEM_FUNCTION_H_

double sum_array(double *array, int size);

int check_sign(double v1, double v2);

/*  Store pure phases composition, gbase, and bulk-rock factor */
typedef struct PP_refs {
	char   	Name[20];			    /** Name                                    */
    double 	Comp[11];        	    /** composition [0-10]                      */
    double 	gbase; 
    double 	gb_lvl;         	    /**driving force, delta_G with G-hyperplane */
    double 	factor;
    double  phase_density;		    /** molar density of the phase              */
    double  phase_shearModulus;		/** molar density of the phase              */
    double  phase_shearModulus_v;
    double  phase_cp;			    /** molar cp of the phase                   */
	double  phase_expansivity;
	double  phase_isoTbulkModulus;
	double  volume_P0;
	double  thetaExp;
	double  phase_entropy;
	double  phase_enthalpy;
	double  phase_bulkModulus;
    double  volume;				    /** molar volume of the phase               */
    double  mass;				    /** molar mass of the phase                 */
    double  charge;
    
} PP_ref;

PP_ref G_EM_function(   int          EM_database, 
                        int          len_ox,
                        int         *id,
                        double      *bulk_rock, 
                        double      *apo,
                        double       P, 
                        double       T, 
                        char        *name, 
                        char        *state          );

/**
  compute the Gibbs Free energy from the thermodynamic database
*/
typedef struct solvent_properties {
    double g;
    double density;
    double epsilon;
    double Z;
} solvent_prop;

/* initialize properties of the solvent, i.e., water in this case */
PP_ref G_FS_function(   int              len_ox,
                        solvent_prop    *wat,
                        int             *id,
                        double          *bulk_rock, 
                        double          *ElH,
                        double          *apo,
                        double           P, 
                        double           T, 
                        char            *name, 
                        char            *state          );

void propSolvent_JN91_calc(     solvent_prop    *wat,
                                double           TK       );

void propSolvent_FE97_calc(     solvent_prop    *wat,
                                double           Pbar,
                                double           TK       );

void propSolvent_SV14_calc(     solvent_prop    *wat,
                                double           Pbar,
                                double           TK       );

void rho_wat_calc(              solvent_prop    *wat,
                                double           Pbar,
                                double           TK,
                                char            *opt      );                                
#endif
