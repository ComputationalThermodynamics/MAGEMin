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
#ifndef __SB_GEM_FUNCTION_H_
#define __SB_GEM_FUNCTION_H_

#include "../gem_function.h"

/* 0 = legacy (Perple_X-style) SLB EOS solver, 1 = burnman-style (Brent volume solve + 3rd order shear),
   2 = same as 1 but with HeFESTo's analytic vibrational/spinodal volume bounds */
void SB_set_eos_formulation(int mode);

/* 0 = legacy compute_G0() behavior unchanged (default), 1 = apply the
   Perple_X-comparison fix: destabilize (NAN) on non-convergence instead of
   silently using the unconverged volume, and tighten the v/v0 sanity bound */
void SB_set_eos_correction(int mode);

PP_ref SB_G_EM_function(    int          EM_database, 
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock, 
                            double      *apo,
                            double       P, 
                            double       T, 
                            char        *name, 
                            char        *state          );
                           
#endif
