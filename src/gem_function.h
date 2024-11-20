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
#ifndef __GEM_FUNCTION_H_
#define __GEM_FUNCTION_H_

/*  Store pure phases composition, gbase, and bulk-rock factor */
typedef struct PP_refs {
	char   	Name[20];			    /** Name                                    */
    double 	Comp[15];        	    /** composition [0-10]                      */
    double 	Comp_mol[15];        	/** composition [0-10]                      */
    double 	Comp_wt[15];        	/** composition [0-10]                      */
    double 	gbase; 
    double 	gb_lvl;         	    /**driving force, delta_G with G-hyperplane */
    double 	factor;
    double  factor_norm;
    double  phase_density;		    /** molar density of the phase              */
    double  phase_shearModulus;		/** molar density of the phase              */
    double  phase_bulkModulus;
    double  phase_cp;			    /** molar cp of the phase                   */
	double  phase_expansivity;
	double  phase_isoTbulkModulus;
	// double  volume_P0;
	double  thetaExp;
	double  phase_entropy;
	double  phase_enthalpy;
    double  volume;				    /** molar volume of the phase               */
    double  mass;				    /** molar mass of the phase                 */
    double  charge;
    
} PP_ref;


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
);

#endif
