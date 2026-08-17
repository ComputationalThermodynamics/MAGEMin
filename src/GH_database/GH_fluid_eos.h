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
#ifndef __GH_FLUID_EOS_H_
#define __GH_FLUID_EOS_H_

double GH_pitzer_sterner_G(int is_H2O, double T, double Pbar);
double GH_pitzer_sterner_H2O_hTs_G(double T, double Pbar);
double GH_pitzer_sterner_mix_G(double x_h2o, double T, double Pbar, double *dGdx_h2o);
double GH_haar_H2O_G(double T, double P);
double GH_wdh78_G(double T, double P);
double GH_duan_CO2_G(double T);

/* Duan & Zhang (2006) H2O-CO2 fluid - the ACTUAL rhyolite-MELTS "fluid"
   solution-phase model (fluidPhase.c's gmixFlu/actFlu engine + the
   "h2oduan"/"co2duan" component standard states propertiesOfPureH2O/CO2).
   Ported 2026-07-17, replacing the Pitzer-Sterner mixture previously
   used by "fl" by mistake (that model is only real MELTS' LIQUID H2O/CO2
   standard-state helper, never its fluid phase). T in K, p in BARS,
   result in J/mol. */
double GH_duan_pure_G(int is_h2o, double t, double p);
double GH_duan_mix_G(double x_h2o, double t, double p, double *dGdx_h2o);
void   GH_duan_mix_muGex(double x_h2o, double t, double p, double *muW, double *muC);

#endif
