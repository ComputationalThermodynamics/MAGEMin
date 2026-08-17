/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __TC_GEM_FUNCTION_H_
#define __TC_GEM_FUNCTION_H_

#include "../gem_function.h"

double sum_array(double *array, int size);

int check_sign(double v1, double v2);

PP_ref TC_G_EM_function(    int          EM_database, 
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock, 
                            double      *apo,
                            double       P, 
                            double       T, 
                            char        *name, 
                            char        *state          );

/**
    Return type for DEW2019 aqueous species (mirrors MAGEMin.jl's AF_ref, DEW_19_gbase.jl).
    Distinct from PP_ref since fl_DEW's solver (aq_min_iterative port) needs each
    species' formation-reaction stoichiometry against the oxide components + H+ (MuComp),
    not just its mass-balance composition (Comp) - no other TC_gem_function.c EOS needs this.
**/
typedef struct AQ_refs {
    char   Name[20];
    double Comp[17];        /** oxide-basis composition (O-corrected), sliced to the active oxide set [0..len_ox-1] */
    double MuComp[18];      /** formation-reaction stoichiometry vs active oxide components [0..len_ox-1] + H+ [len_ox] */
    double z;                /** electrical charge */
    double gbase_supcrt;    /** kJ/mol, SUPCRT-convention G (diagnostic, matches DEW_19_gbase.jl's gbase_supcrt) */
    double gbase;           /** kJ/mol, HSC-convention G - the value the rest of MAGEMin should use (matches gbase_hsc) */
    double factor;          /** always 1.0 - kept for structural parity with PP_ref/AF_ref, not computed from bulk-rock here */
} AQ_ref;

double DEW_psat(         double T_celsius                                          );
double DEW_gSolvent(     double T, double P, double rho_w                          );
double DEW_JN_epsilon(   double T, double rho_w                                    );
double DEW_DM_epsilon(   double T, double rho_w                                    );
double DEW_epsilon(      double T, double P, double rho_w                          );
double DEW_born_B(       double T, double P, double rho_w                          );
double DEW_Agamma(       double T, double P, double rho_w                          );
double DEW_Bgamma(       double T, double P, double rho_w                          );

AQ_ref G_DEW_function(   int              len_ox,
                        int             *id,
                        double          *ElEntropy,
                        double           P,
                        double           T,
                        double           rho_w,
                        char            *name          );

#endif
