/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
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
